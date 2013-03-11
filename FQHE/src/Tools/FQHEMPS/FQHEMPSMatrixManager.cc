////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of manager for FQHE MPS matrices                  //
//                                                                            //
//                        last modification : 31/10/2012                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"

#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"

#include "Options/Options.h"

#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSLaughlinMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RQuasiholeSectorMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3Matrix.h"
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3QuasiholeSectorMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSLaughlinQuasiholeMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSBlockMatrix.h"

#include "Matrix/SparseRealMatrix.h"


#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;


// default constructor 
//

FQHEMPSMatrixManager::FQHEMPSMatrixManager()
{
  this->Options = 0;
}

// destructor
//

FQHEMPSMatrixManager::~FQHEMPSMatrixManager()
{
}
  
// add an option group containing all options related to the MPS matrix construction 
//
// manager = pointer to the option manager
// comment = additional comment that is displayed in the behind each option group

void FQHEMPSMatrixManager::AddOptionGroup(OptionManager* manager, const char* comment)
{
  this->Options = manager;
  char* TmpC = 0;
  if (comment == 0)
    TmpC = new char[255];
  else
    TmpC = new char[255 + strlen(comment)];
  if (comment == 0)
    sprintf(TmpC, "system options");
  else
    sprintf(TmpC, "system options (%s)", comment);  
  OptionGroup* SystemGroup  = new OptionGroup (TmpC);
  if (comment == 0)
    sprintf(TmpC, "precalculation options");
  else
    sprintf(TmpC, "precalculation options (%s)", comment);    
  OptionGroup* PrecalculationGroup = new OptionGroup (TmpC);
  if (comment == 0)
    sprintf(TmpC, "output options");
  else
    sprintf(TmpC, "output options (%s)", comment);    
  OptionGroup* OutputGroup = new OptionGroup (TmpC);
  delete[] TmpC;
  (*(this->Options)) += SystemGroup;
  (*(this->Options)) += PrecalculationGroup;
  (*(this->Options)) += OutputGroup;

  (*SystemGroup) += new SingleIntegerOption  ('\n', "p-truncation", "truncation level", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-index", "index of the Laughlin state to generate", 3);
  (*SystemGroup) += new BooleanOption  ('\n', "k-2", "consider the (k=2,r) series of clustered states");
  (*SystemGroup) += new BooleanOption  ('\n', "rr-3", "consider the k= 3 Read-Rezayi state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "r-index", "r index of the (k,r) clustered state", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new SingleStringOption  ('\n', "cft", "use a file that described the CFT to be used");
  (*SystemGroup) += new BooleanOption  ('\n', "quasihole-sector", "look at the quasihole sector for the (k=2,r) series of clustered states");
  (*SystemGroup) += new SingleStringOption  ('\n', "with-quasiholes", "state has to be built with quasihole whose location is given in a text file");
  (*PrecalculationGroup) += new SingleStringOption('\n', "import-bmatrices", "import the B matrices from a given binary file instead of computing them");
  (*PrecalculationGroup) += new BooleanOption ('\n', "export-bmatrices", "export the B matrices in a binary file");
  (*PrecalculationGroup) += new SingleStringOption('\n', "export-bmatrixname", "use a custom output file name to export the B matrices instead of the default one");
  (*PrecalculationGroup) += new BooleanOption ('\n', "only-export", "only export the B matrices in a binary file and exit from the program");
  (*OutputGroup) += new BooleanOption ('c', "normalize-cylinder", "express the MPS in the normalized cylinder basis");
  (*OutputGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the cylinder", 1);
  (*OutputGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
}


// get the MPS matrice class defined by the running options
//
// nbrFluxQuanta = number of flux quanta
// return value = pointer to the MPS matrice class

AbstractFQHEMPSMatrix* FQHEMPSMatrixManager::GetMPSMatrices(int nbrFluxQuanta)
{
  bool CylinderFlag = this->Options->GetBoolean("normalize-cylinder");
  double AspectRatio = this->Options->GetDouble("aspect-ratio");
  double Kappa = 0.0;
  double Perimeter = 0.0;
  if (CylinderFlag)
    {
      if (this->Options->GetDouble("cylinder-perimeter") > 0.0)
	{
	  Perimeter = this->Options->GetDouble("cylinder-perimeter");
	}
      else
	{
	  Perimeter = sqrt(2.0 * M_PI * (nbrFluxQuanta + 1) * AspectRatio);
	}
      Kappa = (2.0 * M_PI) / Perimeter;
      cout<<"Cylinder geometry, perimeter = " << Perimeter << " , kappa= " << Kappa << endl;
    }

  AbstractFQHEMPSMatrix* MPSMatrix = 0; 
  char* ExportFileName = 0;
  int NbrBMatrices = 2;
  if (this->Options->GetString("with-quasiholes") == 0)
    {
      if (this->Options->GetBoolean("k-2") == true)
	{
	  if (this->Options->GetBoolean("quasihole-sector") == false)
	    {
	      if (this->Options->GetString("import-bmatrices") != 0)
		{
		  MPSMatrix = new FQHEMPSClustered2RMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), 
							   this->Options->GetString("import-bmatrices"), CylinderFlag, Kappa);
		}
	      else
		{
		  if (this->Options->GetString("cft") != 0)
		    {
		      MPSMatrix = new FQHEMPSClustered2RMatrix(this->Options->GetInteger("p-truncation"), NbrBMatrices, this->Options->GetString("cft"),
							       CylinderFlag, Kappa);
		    }
		  else
		    {
		      MPSMatrix = new FQHEMPSClustered2RMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), NbrBMatrices,
							       CylinderFlag, Kappa);
		    }
		}
	    }
	  else
	    {
	      if (this->Options->GetString("import-bmatrices") != 0)
		{
		  MPSMatrix = new FQHEMPSClustered2RQuasiholeSectorMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), 
									      this->Options->GetString("import-bmatrices"), CylinderFlag, Kappa);
		}
	      else
		{
		  if (this->Options->GetString("cft") != 0)
		    {
		      MPSMatrix = new FQHEMPSClustered2RQuasiholeSectorMatrix(this->Options->GetInteger("p-truncation"), NbrBMatrices, this->Options->GetString("cft"), 
									      CylinderFlag, Kappa);
		    }
		  else
		    {
		      MPSMatrix = new FQHEMPSClustered2RQuasiholeSectorMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), NbrBMatrices,
									      CylinderFlag, Kappa);
		    }
// 		  AbstractFQHEMPSMatrix* MPSMatrix1 = new FQHEMPSClustered2RMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), NbrBMatrices,
// 										   CylinderFlag, Kappa);
// 		  AbstractFQHEMPSMatrix* MPSMatrix2 = new FQHEMPSClustered2RQuasiholeSectorMatrix(this->Options->GetInteger("r-index"), 2, this->Options->GetInteger("p-truncation"), NbrBMatrices,
// 												  CylinderFlag, Kappa);
// 		  MPSMatrix = new FQHEMPSBlockMatrix(MPSMatrix1, MPSMatrix2);
// 		  delete MPSMatrix1;
// 		  delete MPSMatrix2;
		}
	    }	  
	  ExportFileName = new char [512];
	  if (CylinderFlag == false)
	    {
	      sprintf(ExportFileName, "fqhemps_bmatrices_unnormalized_clustered_k_2_r_%ld_q_2_p_%ld_n_%d.dat", this->Options->GetInteger("r-index"), this->Options->GetInteger("p-truncation"), NbrBMatrices);
	    }
	  else
	    {
	      sprintf(ExportFileName, "fqhemps_bmatrices_cylinder_clustered_k_2_r_%ld_q_2_p_%ld_n_%d_perimeter_%.6f.dat", this->Options->GetInteger("r-index"), this->Options->GetInteger("p-truncation"), NbrBMatrices, Perimeter);
	    }
	}
      else
	{
	  if (this->Options->GetBoolean("rr-3") == true)
	    {
	      if (this->Options->GetBoolean("quasihole-sector") == false)
		{
		  if (this->Options->GetString("import-bmatrices") != 0)
		    {
		      MPSMatrix = new FQHEMPSReadRezayi3Matrix(2, this->Options->GetInteger("p-truncation"), this->Options->GetString("import-bmatrices"), 
							       CylinderFlag, Kappa);
		    }
		  else
		    {
		      MPSMatrix = new FQHEMPSReadRezayi3Matrix(2, this->Options->GetInteger("p-truncation"), NbrBMatrices,
							       CylinderFlag, Kappa);
		    }
		}
	      else
		{
		  if (this->Options->GetString("import-bmatrices") != 0)
		    {
		      MPSMatrix = new FQHEMPSReadRezayi3QuasiholeSectorMatrix(2, this->Options->GetInteger("p-truncation"), this->Options->GetString("import-bmatrices"), 
									      CylinderFlag, Kappa);
		    }
		  else
		    {
		      MPSMatrix = new FQHEMPSReadRezayi3QuasiholeSectorMatrix(2, this->Options->GetInteger("p-truncation"), NbrBMatrices,
									      CylinderFlag, Kappa);
		    }
		}
	      ExportFileName = new char [512];
	      if (CylinderFlag == false)
		{
		  sprintf(ExportFileName, "fqhemps_bmatrices_unnormalized_readrezayi3_q_2_p_%ld_n_%d.dat", this->Options->GetInteger("p-truncation"), NbrBMatrices);
		}
	      else
		{
		  sprintf(ExportFileName, "fqhemps_bmatrices_cylinder_readrezayi3_q_2_p_%ld_n_%d_perimeter_%.6f.dat", this->Options->GetInteger("p-truncation"), NbrBMatrices, Perimeter);
		}
	    }
	  else
	    {
	      if (this->Options->GetString("import-bmatrices") != 0)
		{
		  MPSMatrix = new FQHEMPSLaughlinMatrix(this->Options->GetInteger("laughlin-index"), this->Options->GetInteger("p-truncation"), this->Options->GetString("import-bmatrices"),
							CylinderFlag, Kappa);
		}
	      else
		{
		  MPSMatrix = new FQHEMPSLaughlinMatrix(this->Options->GetInteger("laughlin-index"), this->Options->GetInteger("p-truncation"), NbrBMatrices,
							CylinderFlag, Kappa);
		}
	      ExportFileName = new char [512];
	      if (CylinderFlag == false)
		{
		  sprintf(ExportFileName, "fqhemps_bmatrices_unnormalized_laughlin_q_%ld_p_%ld_n_%d.dat", this->Options->GetInteger("laughlin-index"), 
			  this->Options->GetInteger("p-truncation"), NbrBMatrices);
		}
	      else
		{
		  sprintf(ExportFileName, "fqhemps_bmatrices_cylinder_laughlin_q_%ld_p_%ld_n_%d_perimeter_%.6f.dat", this->Options->GetInteger("laughlin-index"), 
			  this->Options->GetInteger("p-truncation"), NbrBMatrices, Perimeter);
		}
	    }
	}
    }
  else
    {
      if (this->Options->GetBoolean("k-2") == true)
	{
	  MPSMatrix = 0;
	}
      else
	{
	  if (this->Options->GetBoolean("rr-3") == true)
	    {
	      MPSMatrix = 0;
	    }
	  else
	    {
	      if (this->Options->GetString("import-bmatrices") != 0)
		{
		  MPSMatrix = new FQHEMPSLaughlinQuasiholeMatrix(this->Options->GetInteger("laughlin-index"), this->Options->GetInteger("p-truncation"), this->Options->GetString("import-bmatrices"),
							CylinderFlag, Kappa);
		}
	      else
		{
		  MPSMatrix = new FQHEMPSLaughlinQuasiholeMatrix(this->Options->GetInteger("laughlin-index"), this->Options->GetInteger("p-truncation"), NbrBMatrices,
							CylinderFlag, Kappa);
		}
	    }
	}
    }
  if (this->Options->GetBoolean("export-bmatrices"))
    {
      if (this->Options->GetString("export-bmatrixname") != 0)
	MPSMatrix->SaveMatrices(this->Options->GetString("export-bmatrixname"));
      else
	MPSMatrix->SaveMatrices(ExportFileName);
      cout << "number of non-zero matrix elements:" << endl;
      for (int i = 0; i < NbrBMatrices; ++i)
	cout << "B[" << i << "] = " << MPSMatrix->GetMatrices()[i].ComputeNbrNonZeroMatrixElements() << endl;
    }

  delete[] ExportFileName;
  return MPSMatrix;
}


// get the cylinder perimeter (in magnetic length unit) if the cylinder geometry if used
//
// nbrFluxQuanta = number of flux quanta
// return value = cylinder perimeter (negative if another geometry is used)

double FQHEMPSMatrixManager::GetCylinderPerimeter(int nbrFluxQuanta)
{
  if (this->Options->GetBoolean("normalize-cylinder"))
    {
      if (this->Options->GetDouble("cylinder-perimeter") > 0.0)
	return this->Options->GetDouble("cylinder-perimeter");
      else
	return (sqrt(2.0 * M_PI * (nbrFluxQuanta + 1) * this->Options->GetDouble("aspect-ratio")));
    }
  return -1.0;
}
