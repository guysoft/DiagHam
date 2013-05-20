////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the clustered (k=2,r) states            //
//                         in their quasihole sector                          //
//                                                                            //
//                        last modification : 11/02/2013                      //
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
#include "Tools/FQHEMPS/FQHEMPSClustered2RQuasiholeSectorMatrix.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/LongRationalMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "Architecture/ArchitectureOperation/FQHEMPSEvaluateCFTOperation.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor 
//

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix()
{
}

// constructor 
//
// rindex = r index (i.e. clustered (k=2,r) states) 
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, bool cylinderFlag, double kappa, 
										 AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->CentralCharge = LongRational ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SelfDualFlag = ((this->RIndex & 1) == 0);
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
  this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));
  this->MatrixElementNormalization = 1.0 / M_SQRT2;
  this->SquareMatrixElementNormalization = LongRational(1, 2);
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  if (this->SelfDualFlag == true)
    this->TransferMatrixDegeneracy /= 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
  this->CreateBMatrices(0, architecture);
}

// constructor 
//
// rindex = r index (i.e. clustered (k=2,r) states) 
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// cftDirectory = path to the directory where all the pure CFT matrices are stored
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, 
										 bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->CentralCharge = LongRational ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SelfDualFlag = ((this->RIndex & 1) == 0);
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
  this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));
  this->MatrixElementNormalization = 1.0 / M_SQRT2;
  this->SquareMatrixElementNormalization = LongRational(1, 2);
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
  this->CreateBMatrices(cftDirectory, architecture);
}

// constructor from a file describing the state
//
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// fileName = name of the file that contains the state description
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix(int pLevel, int nbrBMatrices, char* fileName, bool cylinderFlag, double kappa, 
										 AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->PLevel = pLevel;
  
  ConfigurationParser StateDefinition;
  if (StateDefinition.Parse(fileName) == false)
    {
      StateDefinition.DumpErrors(cout) << endl;
    }
  else
    {
      bool ErrorFlag = true;
      ErrorFlag = StateDefinition.GetAsSingleInteger("RIndex", this->RIndex);
      ErrorFlag = StateDefinition.GetAsSingleInteger("LaughlinIndex", this->LaughlinIndex);
      ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightIdentity", this->WeightIdentity);
      ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightPsi", this->WeightPsi);
      ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightSigma", this->WeightSigma);
      ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightPhi", this->WeightPhi);
      ErrorFlag = StateDefinition.GetAsBoolean("SelfDual", this->SelfDualFlag);
      ErrorFlag = StateDefinition.GetAsSingleLongRational("CentralCharge", this->CentralCharge);
      if (StateDefinition["PsiSquareMatrixElement"] != 0)
	{
	  ErrorFlag = StateDefinition.GetAsSingleLongRational("PsiSquareMatrixElement", this->SquareMatrixElementNormalization);
	}
      else
	{
	  this->SquareMatrixElementNormalization = LongRational(1, 2);
	}
      this->MatrixElementNormalization = sqrt(fabs(this->SquareMatrixElementNormalization.GetNumericalValue()));
      if (StateDefinition["EMatrixDegeneracy"] != 0)
	{
	  ErrorFlag = StateDefinition.GetAsSingleInteger("EMatrixDegeneracy", this->TransferMatrixDegeneracy);	  
	}
      else
	{
	  switch (this->RIndex)
	    {
	    case 2:
	      this->TransferMatrixDegeneracy = 2;
	      break;
	    case 3:
	      this->TransferMatrixDegeneracy = 5;
	      break;
	    case 6:
	      this->TransferMatrixDegeneracy = 5;
	      break;
	    }
	}
      if (StateDefinition["Name"] != 0)
	{
	  this->BMatrixOutputName = new char[strlen(StateDefinition["Name"]) + 1]; 
	  strcpy(this->BMatrixOutputName, StateDefinition["Name"]);
	}
      else
	{
	  this->BMatrixOutputName = new char[256]; 
	  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
	}
      if (ErrorFlag == true)
	{
	  this->WeightPrimaryFieldMatrixElement = this->WeightPsi;
	  this->CreateBMatrices(StateDefinition["CFTMatrixDirectory"], architecture);
	}
    }
}


// constructor from stored B matrices
//
// rindex = r index (i.e. clustered (k=2,r) states) 
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag, double kappa)
{
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->LoadMatrices(fileName);
  this->CentralCharge = LongRational ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SelfDualFlag = ((this->RIndex & 1) == 0);
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
  this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));
  this->MatrixElementNormalization = 1.0 / M_SQRT2;
  this->SquareMatrixElementNormalization = LongRational(1, 2);
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
}

// destructor
//

FQHEMPSClustered2RQuasiholeSectorMatrix::~FQHEMPSClustered2RQuasiholeSectorMatrix()
{
}
  
// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSClustered2RQuasiholeSectorMatrix::CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture)
{
  LongRational CentralCharge12 (this->CentralCharge);
  cout << "central charge = " << CentralCharge12 << endl;
  CentralCharge12 /= 12l;
  double WeightSigmaNumerical = this->WeightSigma.GetNumericalValue();
  double WeightPhiNumerical = this->WeightPhi.GetNumericalValue();
  double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();
  double WeightPrimaryFieldMatrixElementNumerical = this->WeightPrimaryFieldMatrixElement.GetNumericalValue();
  long* Partition = new long[2 * (this->PLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevel + 2];
  this->TemporaryOccupationNumber = new unsigned long [this->PLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductSigma = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductPhi = new RealSymmetricMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductSigma = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductPhi = new LongRationalMatrix[this->PLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi01 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi10 = new LongRationalMatrix*[this->PLevel + 1];
  RealMatrix* OrthogonalBasisSigmaLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPhiLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisSigmaRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPhiRight = new RealMatrix[this->PLevel + 1];
  LongRational** RationalMultiplicityFactor = new LongRational*[this->PLevel + 1];
  double** MultiplicityFactor = new double*[this->PLevel + 1];

  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, this->PLevel + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevel + 1];
      RationalMatrixPsi01[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi10[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMultiplicityFactor[i] = new LongRational[U1BosonBasis[i]->GetHilbertSpaceDimension()];
      MultiplicityFactor[i] = new double[U1BosonBasis[i]->GetHilbertSpaceDimension()];
      for (int j = 0; j < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++j)
	{
	  U1BosonBasis[i]->GetOccupationNumber(j, TmpPartition);	    
	  RationalMultiplicityFactor[i][j] = 1l;
	  MultiplicityFactor[i][j] = 1.0;
 	  for (int k = 1; k <= i; ++k)
	    {
	      if (TmpPartition[k] > 1ul)
		{
		  RationalMultiplicityFactor[i][j].FactorialDivide(TmpPartition[k]);
		  double Tmp = 1.0;
		  for (unsigned long l = 2l; l <= TmpPartition[k]; ++l)
		    Tmp *=  (double) l;
		  MultiplicityFactor[i][j] /= Tmp;
		}
	    }
	}
    }

  cout << "weight: " <<   this->WeightSigma << " " << this->WeightPhi << endl;
  char* TmpScalarProductSigmaFileName = 0; 
  char* TmpScalarProductPhiFileName = 0;
  if (cftDirectory != 0)
    {
      TmpScalarProductSigmaFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductPhiFileName = new char[512 + strlen(cftDirectory)];
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      if (cftDirectory != 0)
	{
	  if (this->UseRationalFlag == true)
	    {
	      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_scalarproducts_sigma_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	      if (this->SelfDualFlag == false)
		sprintf (TmpScalarProductPhiFileName, "%s/cft_%s_scalarproducts_phi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	    }
	  else
	    {
	      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_num_scalarproducts_sigma_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	      if (this->SelfDualFlag == false)
		sprintf (TmpScalarProductPhiFileName, "%s/cft_%s_num_scalarproducts_phi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	    }
	}
      if ((cftDirectory != 0) && (IsFile(TmpScalarProductSigmaFileName)))
	{		
	  if (this->UseRationalFlag == true)
	    {
	      RationalScalarProductSigma[i].ReadMatrix(TmpScalarProductSigmaFileName);
	    }
	  else
	    {
	      ScalarProductSigma[i].ReadMatrix(TmpScalarProductSigmaFileName);
	    }
	}
      else
	{
	  if (this->UseRationalFlag == true)
	    {
	      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12, 
						     this->WeightSigma,
						     RationalScalarProductSigma,  i- 1);
	      Operation1.ApplyOperation(architecture);
	      RationalScalarProductSigma[i] = Operation1.GetRationalMatrixElements();
	      if  ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  RationalScalarProductSigma[i].WriteMatrix(TmpScalarProductSigmaFileName);
		}
	    }
	  else
	    {
	      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12Numerical, 
						     WeightSigmaNumerical,
						     ScalarProductSigma,  i- 1);
	      Operation1.ApplyOperation(architecture);
	      ScalarProductSigma[i] = Operation1.GetMatrixElements();
	      if  ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  ScalarProductSigma[i].WriteMatrix(TmpScalarProductSigmaFileName);
		}
	    }
	}
      if (this->SelfDualFlag == false)
	{
	  if ((cftDirectory != 0) && (IsFile(TmpScalarProductPhiFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalScalarProductPhi[i].ReadMatrix(TmpScalarProductPhiFileName);
		}
	      else
		{
		  ScalarProductPhi[i].ReadMatrix(TmpScalarProductPhiFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12, 
							 this->WeightPhi,
							 RationalScalarProductPhi,  i - 1);
		  Operation2.ApplyOperation(architecture);
		  RationalScalarProductPhi[i] = Operation2.GetRationalMatrixElements();
		  if  ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalScalarProductPhi[i].WriteMatrix(TmpScalarProductPhiFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12Numerical, 
							 WeightPhiNumerical,
							 ScalarProductPhi,  i - 1);
		  Operation2.ApplyOperation(architecture);
		  ScalarProductPhi[i] = Operation2.GetRationalMatrixElements();
		  if  ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      ScalarProductPhi[i].WriteMatrix(TmpScalarProductPhiFileName);
		    }
		}
	    }
	}
      
      RealSymmetricMatrix TmpMatrix;
      if (this->UseRationalFlag == true)
 	{
 	  LongRationalMatrix TmpRationalMatrix(RationalScalarProductSigma[i].GetNbrRow(), RationalScalarProductSigma[i].GetNbrColumn());
   	  for (int k = 0; k < RationalScalarProductSigma[i].GetNbrRow(); ++k)
   	    for (int l = 0; l < RationalScalarProductSigma[i].GetNbrColumn(); ++l)
   	      {
   		TmpRationalMatrix[l][k] = RationalScalarProductSigma[i][l][k] * (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
   	      }
 	  TmpMatrix = TmpRationalMatrix;
  	}
      else
	{
	  TmpMatrix = RealSymmetricMatrix (ScalarProductSigma[i].GetNbrRow(), ScalarProductSigma[i].GetNbrColumn());
	  for (int k = 0; k < ScalarProductSigma[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductSigma[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductSigma[i].GetMatrixElement(k, l, Tmp);
		Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		TmpMatrix.SetMatrixElement(k, l, Tmp);
	      }
	}
      
      RealMatrix TmpBasis(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension());
      TmpBasis.SetToIdentity();
      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
      TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
      TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
      double Error = 0.0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) > Error)
	  Error = fabs(TmpDiag(n, n));
      Error *= 1e-14;
      if (Error < 1e-14)
	Error = 1e-14;
      int Count  = 0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	{
	  if (fabs(TmpDiag(n, n)) < Error)
	    ++Count;
	}
      cout << "nbr of null vectors identity sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisSigmaLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  OrthogonalBasisSigmaRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisSigmaLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisSigmaRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisSigmaLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisSigmaRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisSigmaLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisSigmaRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisSigmaLeft[i] = RealMatrix();
	  OrthogonalBasisSigmaRight[i] = RealMatrix();
	}

      if (this->SelfDualFlag == false)
	{
	  if (this->UseRationalFlag == true)
	    {
	      LongRationalMatrix TmpRationalMatrix2(RationalScalarProductPhi[i].GetNbrRow(), RationalScalarProductPhi[i].GetNbrColumn());
	      for (int k = 0; k < RationalScalarProductPhi[i].GetNbrRow(); ++k)
		for (int l = 0; l < RationalScalarProductPhi[i].GetNbrColumn(); ++l)
		  {
		    TmpRationalMatrix2[l][k] = RationalScalarProductPhi[i][l][k] * (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
		  }
	      TmpMatrix = TmpRationalMatrix2;
	    }
	  else
	    {
	      TmpMatrix = RealSymmetricMatrix (ScalarProductPhi[i].GetNbrRow(), ScalarProductPhi[i].GetNbrColumn());
	      for (int k = 0; k < ScalarProductPhi[i].GetNbrRow(); ++k)
		for (int l = k; l < ScalarProductPhi[i].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    ScalarProductPhi[i].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		    TmpMatrix.SetMatrixElement(k, l, Tmp);
		  }
	    }
	  TmpBasis.SetToIdentity();
#ifdef __LAPACK__
	  TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
	  TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
	  Error = 0.0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      Error = fabs(TmpDiag(n, n));
	  Error *= 1e-14;
	  if (Error < 1e-14)
	    Error = 1e-14;
	  Count  = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) < Error)
	      ++Count;
	  cout << "nbr of null vectors Phi sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
	  if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	    {
	      OrthogonalBasisPhiLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	      OrthogonalBasisPhiRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	      Count = 0;
	      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		if (fabs(TmpDiag(n, n)) > Error)
		  {
		    OrthogonalBasisPhiLeft[i][Count].Copy(TmpBasis[n]);
		    OrthogonalBasisPhiRight[i][Count].Copy(TmpBasis[n]);
		    if (TmpDiag(n, n) > 0)
		      {
			OrthogonalBasisPhiLeft[i][Count] /=  sqrt(TmpDiag(n, n));
			OrthogonalBasisPhiRight[i][Count] /=  sqrt(TmpDiag(n, n));
		      }
		    else
		      {
			OrthogonalBasisPhiLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
			OrthogonalBasisPhiRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		      }
		    ++Count;
		  }
	    }
	  else
	    {
	      OrthogonalBasisPhiLeft[i] = RealMatrix();
	      OrthogonalBasisPhiRight[i] = RealMatrix();
	    }
	}
      cout << "---------------------------------" << endl;
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      if (this->UseRationalFlag == true)
 	{
    	  for (int k = 0; k < RationalScalarProductSigma[i].GetNbrRow(); ++k)
   	    for (int l = 0; l < RationalScalarProductSigma[i].GetNbrColumn(); ++l)
   	      {
   		RationalScalarProductSigma[i][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
   	      }
 	  ScalarProductSigma[i] = RationalScalarProductSigma[i];
  	}
      else
	{
	  for (int k = 0; k < ScalarProductSigma[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductSigma[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductSigma[i].GetMatrixElement(k, l, Tmp);
		Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		ScalarProductSigma[i].SetMatrixElement(k, l, Tmp);
	      }
	}
      if (this->SelfDualFlag == false)
	{
	  if (this->UseRationalFlag == true)
	    {
	      for (int k = 0; k < RationalScalarProductPhi[i].GetNbrRow(); ++k)
		for (int l = 0; l < RationalScalarProductPhi[i].GetNbrColumn(); ++l)
		  {
		    RationalScalarProductPhi[i][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
		  }
	      ScalarProductPhi[i] = RationalScalarProductPhi[i];
	    }
	  else
	    {
	      for (int k = 0; k < ScalarProductPhi[i].GetNbrRow(); ++k)
		for (int l = k; l < ScalarProductPhi[i].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    ScalarProductPhi[i].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		    ScalarProductPhi[i].SetMatrixElement(k, l, Tmp);
		  }
	    }
	}
   }
  this->StartingIndexPerPLevel = new int* [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->IdentityBasisDimension = new int [this->PLevel + 1];
  this->PsiBasisDimension = new int [this->PLevel + 1];
  this->U1BasisDimension = new int [this->PLevel + 1];	
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->IdentityBasisDimension[i] = OrthogonalBasisSigmaLeft[i].GetNbrColumn();
      if (this->SelfDualFlag == false)
	{
	  this->PsiBasisDimension[i] = OrthogonalBasisPhiLeft[i].GetNbrColumn();
	}
      else
	{
	  this->PsiBasisDimension[i] = 0.0;
	}
      this->U1BasisDimension[i] = U1BosonBasis[i]->GetHilbertSpaceDimension();
    }
  this->TotalStartingIndexPerPLevel[0] = 0;
  this->StartingIndexPerPLevel[0] = new int [1];
  this->StartingIndexPerPLevel[0][0] = 0;

  
  int NValueShift;
  int QValue;
  int QValueDenominator;
  double ExtraCylinderFactor = 1.0;
  if ((this->RIndex & 1) == 0)
    {
      QValue = 1 + (this->RIndex / 2);
      this->NbrNValue = ((2 * this->PLevel) + QValue) + this->RIndex / 2 + 1;
      NValueShift = 2 * this->PLevel - 1;
      --this->NbrNValue;
      ++NValueShift;
      QValueDenominator = 1;
    }
  else
    {
      QValue = 2 + this->RIndex;
      this->NbrNValue = ((4 * this->PLevel) + QValue) + this->RIndex + 1;
      NValueShift = 4 * this->PLevel - 2;
      this->NbrNValue -= 2;
      NValueShift += 2;
      QValueDenominator = 2;
      ExtraCylinderFactor = 4.0;
    }

  this->NbrNValuesPerPLevel = new int [this->PLevel + 1];
  this->NInitialValuePerPLevel = new int [this->PLevel + 1];
  this->NLastValuePerPLevel = new int [this->PLevel + 1];     
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->ComputeChargeIndexRange(i, this->NInitialValuePerPLevel[i], this->NLastValuePerPLevel[i]);
      this->NbrNValuesPerPLevel[i] =  this->NLastValuePerPLevel[i] - this->NInitialValuePerPLevel[i] + 1;
    }
     
  if (this->SelfDualFlag == false)
    {
      this->NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * (OrthogonalBasisSigmaLeft[0].GetNbrColumn() + OrthogonalBasisPhiLeft[0].GetNbrColumn())) * this->NbrNValue;
    }
  else
    {
      this->NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * OrthogonalBasisSigmaLeft[0].GetNbrColumn()) * this->NbrNValue;
    }
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      this->StartingIndexPerPLevel[i] = new int [i + 1];      
      this->StartingIndexPerPLevel[i][0] = this->TotalStartingIndexPerPLevel[i];
      int Tmp = 0;
      int Tmp2;
      if (this->SelfDualFlag == false)
	{
	  for (int j = 0; j < i; ++j)
	    {
	      Tmp2 = U1BosonBasis[j]->GetHilbertSpaceDimension() * (OrthogonalBasisSigmaLeft[i - j].GetNbrColumn() + OrthogonalBasisPhiLeft[i - j].GetNbrColumn()) * this->NbrNValue;
	      this->StartingIndexPerPLevel[i][j + 1] = Tmp2 + this->StartingIndexPerPLevel[i][j];
	      Tmp += Tmp2;
	    }
	  Tmp += U1BosonBasis[i]->GetHilbertSpaceDimension() * (OrthogonalBasisSigmaLeft[0].GetNbrColumn() + OrthogonalBasisPhiLeft[0].GetNbrColumn()) * this->NbrNValue;
	}
      else
	{
	  for (int j = 0; j < i; ++j)
	    {
	      Tmp2 = U1BosonBasis[j]->GetHilbertSpaceDimension() * OrthogonalBasisSigmaLeft[i - j].GetNbrColumn() * this->NbrNValue;
	      this->StartingIndexPerPLevel[i][j + 1] = Tmp2 + this->StartingIndexPerPLevel[i][j];
	      Tmp += Tmp2;
	    }
	  Tmp += U1BosonBasis[i]->GetHilbertSpaceDimension() * OrthogonalBasisSigmaLeft[0].GetNbrColumn() * this->NbrNValue;
	}

      this->NbrIndicesPerPLevel[i] =  Tmp;
    }
  
  int MatrixSize = this->NbrIndicesPerPLevel[this->PLevel] + this->TotalStartingIndexPerPLevel[this->PLevel];
  cout << "B matrix size = " << MatrixSize << endl;

  cout << "computing Psi matrix elements" << endl;
  LongRational Weight (this->WeightPsi);
  for (int j = 0; j <= this->PLevel; ++j)
    {
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  cout << "Levels = " <<  i << " " << j << endl;
	  if (cftDirectory != 0)
	    {
	      if (this->UseRationalFlag == true)
		{
		  if (this->SelfDualFlag == false)
		    {
		      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_matrixelement_sigmaphi_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		      sprintf (TmpScalarProductPhiFileName, "%s/cft_%s_matrixelement_phisigma_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		    }
		  else
		    {
		      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_matrixelement_sigmasigma_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		    }
		}
	      else
		{
		  if (this->SelfDualFlag == false)
		    {
		      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_num_matrixelement_sigmaphi_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		      sprintf (TmpScalarProductPhiFileName, "%s/cft_%s_num_matrixelement_phisigma_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		    }
		  else
		    {
		      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_num_matrixelement_sigmasigma_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		    }
		}
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpScalarProductSigmaFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi01[i][j].ReadMatrix(TmpScalarProductSigmaFileName);
		}
	      else
		{
		  MatrixPsi01[i][j].ReadMatrix(TmpScalarProductSigmaFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, j, CentralCharge12, 
							 this->WeightSigma, this->WeightPhi, this->WeightPrimaryFieldMatrixElement,
							 RationalMatrixPsi01,  i - 1, j);
		  Operation1.ApplyOperation(architecture);
		  RationalMatrixPsi01[i][j] = Operation1.GetRationalMatrixElements();
		  if  ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi01[i][j].WriteMatrix(TmpScalarProductSigmaFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							 WeightSigmaNumerical, WeightPhiNumerical, WeightPrimaryFieldMatrixElementNumerical,
							 MatrixPsi01,  i - 1, j);
		  Operation1.ApplyOperation(architecture);
		  MatrixPsi01[i][j] = Operation1.GetMatrixElements();
		  if  ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi01[i][j].WriteMatrix(TmpScalarProductSigmaFileName);
		    }
		}
	    }
	  if (this->SelfDualFlag == false)
	    {
	      if ((cftDirectory != 0) && (IsFile(TmpScalarProductPhiFileName)))
		{
		  if (this->UseRationalFlag == true)
		    {
		      RationalMatrixPsi10[i][j].ReadMatrix(TmpScalarProductPhiFileName);
		    }
		  else
		    {
		      MatrixPsi10[i][j].ReadMatrix(TmpScalarProductPhiFileName);
		    }
		}
	      else
		{
		  if (this->UseRationalFlag == true)
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							     this->WeightPhi, this->WeightSigma, this->WeightPrimaryFieldMatrixElement,
							     RationalMatrixPsi10,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      RationalMatrixPsi10[i][j] = Operation2.GetRationalMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalMatrixPsi10[i][j].WriteMatrix(TmpScalarProductPhiFileName);
			}
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							     WeightPhiNumerical, WeightSigmaNumerical, WeightPrimaryFieldMatrixElementNumerical,
							     MatrixPsi10,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      MatrixPsi10[i][j] = Operation2.GetMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  MatrixPsi10[i][j].WriteMatrix(TmpScalarProductPhiFileName);
			}
		    }
		}
	    }
	}
    }
  for (int j = 0; j <= this->PLevel; ++j)
    {
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  if (this->UseRationalFlag == true)
	    {
	      for (int k = 0; k < RationalMatrixPsi01[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < RationalMatrixPsi01[i][j].GetNbrColumn(); ++l)
		  {
		    RationalMatrixPsi01[i][j][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[j][l]);
		  }
	      MatrixPsi01[i][j] = RationalMatrixPsi01[i][j];
	      if (this->SelfDualFlag == false)
		{
		  for (int k = 0; k < RationalMatrixPsi10[i][j].GetNbrRow(); ++k)
		    for (int l = 0; l < RationalMatrixPsi10[i][j].GetNbrColumn(); ++l)
		      {
			RationalMatrixPsi10[i][j][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[j][l]);
		      }
		  MatrixPsi10[i][j] = RationalMatrixPsi10[i][j];
		}
	    }
	  else
	    {
	      for (int k = 0; k < MatrixPsi01[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < MatrixPsi01[i][j].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    MatrixPsi01[i][j].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[j][l]);
		    MatrixPsi01[i][j].SetMatrixElement(k, l, Tmp);
		  }
	      if (this->SelfDualFlag == false)
		{
		  for (int k = 0; k < MatrixPsi10[i][j].GetNbrRow(); ++k)
		    for (int l = 0; l < MatrixPsi10[i][j].GetNbrColumn(); ++l)
		      {
			double Tmp;
			MatrixPsi10[i][j].GetMatrixElement(k, l, Tmp);
			Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[j][l]);
			MatrixPsi10[i][j].SetMatrixElement(k, l, Tmp);
		      }
		}
	    }
	  MatrixPsi01[i][j] *= this->MatrixElementNormalization;
	  if (this->SelfDualFlag == false)
	    {
	      MatrixPsi10[i][j] *= this->MatrixElementNormalization;
	    }
	}
    }
  cout << "building B matrices" << endl;


  SparseRealMatrix* BMatrices = new SparseRealMatrix[this->NbrBMatrices];
  int* TmpNbrElementPerRow = new int[MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;

  // B^[0]  matrix evaluation

  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaRight = OrthogonalBasisSigmaRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductSigma = ScalarProductSigma[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      if ((this->RIndex & 1) == 0)
		{
		  for (int j = 1; j < this->NbrNValue; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndex(j - 1, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
			    }
			}
		      if (this->SelfDualFlag == false)
			{
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex2)
				{
				  ++TmpNbrElementPerRow[this->Get2RMatrixIndex(j - 1, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
				}
			    }
			}
		    }
		}
	      else
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  for (int j = 2; j < this->NbrNValue; ++j)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndex(j - 2, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
			    }
			}
		    }
		  if (this->SelfDualFlag == false)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      for (int j = 2; j < this->NbrNValue; ++j)
				{
				  ++TmpNbrElementPerRow[this->Get2RMatrixIndex(j - 2, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
				}
			    }
			}
		    }
		}	      
	    }
	}
    }
  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaRight = OrthogonalBasisSigmaRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductSigma = ScalarProductSigma[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      if ((this->RIndex & 1) == 0)
		{
		  for (int j = 1; j < this->NbrNValue; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductSigma(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigmaRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisSigmaLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      if (this->CylinderFlag)
				Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
									 + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 1, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
		      if (this->SelfDualFlag == false)
			{
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpScalarProductPhi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPhiRight(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisPhiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  if (this->CylinderFlag)
				    Tmp *= exp(-this->Kappa * this->Kappa * (WeightPhiNumerical +  ((double) i)
									     + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									     + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 1, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
								this->Get2RMatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), Tmp);
				}
			    }
			}
		    }
		}
	      else
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  for (int j = 2; j < this->NbrNValue; ++j)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductSigma(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigmaRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisSigmaLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      if (this->CylinderFlag)
				Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
									 + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 2, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
		    }
		  if (this->SelfDualFlag == false)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      for (int j = 2; j < this->NbrNValue; ++j)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpScalarProductPhi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPhiRight(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisPhiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  if (this->CylinderFlag)
				    Tmp *= exp(-this->Kappa * this->Kappa * (WeightPhiNumerical +  ((double) i)
									     + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									     + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
				  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 2, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
								this->Get2RMatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), Tmp);
				}
			    }
			}
		    }
		}	      
	    }
	}
    }
  
   // B^[1]  matrix evaluation
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;
  FactorialCoefficient Coef;
  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[i - p][j - q];
		  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[i - p][j - q];	
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  if ((this->RIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + this->RIndex + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (4 * (j - i) + 2 * this->RIndex + NValueShift) / 2;
			      N1 = N2 + QValue - 2;
			    }	
			  if (this->SelfDualFlag == false)
			    {		  
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhi2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      ++TmpNbrElementPerRow[this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisSigma1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
				    } 
				}
			    }
			  else
			    {		  
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      ++TmpNbrElementPerRow[this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisSigma1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
				    } 
				}
			    }

			  if ((this->RIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + 2 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (4 * (j - i) + 4 + NValueShift) / 2;
			      N1 = N2 + QValue - 2;
			    }
			  if (this->SelfDualFlag == false)
			    {		  
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhi1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      ++TmpNbrElementPerRow[this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisSigma1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
				    }
				}
			    }
			}
		    }
		}	      
	    }
	}
    }

  BMatrices[1] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[i - p][j - q];
		  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[i - p][j - q];	
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  if ((this->RIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + this->RIndex + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (4 * (j - i) + 2 * this->RIndex + NValueShift) / 2;
			      N1 = N2 + QValue - 2;
			    }	
			  if (this->SelfDualFlag == false)
			    {		  
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhi2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      double Tmp = 0.0;
				      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					{
					  double Tmp1 = 0.0;			      
					  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					    {
					      Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPhi2(NeutralIndex4, NeutralIndex2);				  
					    }
					  Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
					}
				      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				      if (this->CylinderFlag)
					Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightSigmaNumerical + WeightPhiNumerical + ((double) (i + j))
										       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))
										       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift)) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))));
				      BMatrices[1].SetMatrixElement(this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisSigma1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
								    this->Get2RMatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisSigma2.GetNbrColumn(), this->StartingIndexPerPLevel[j][q]), Tmp);
				    } 
				}
			    }
			  else
			    {		  
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      double Tmp = 0.0;
				      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					{
					  double Tmp1 = 0.0;			      
					  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					    {
					      Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
					    }
					  Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
					}
				      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				      if (this->CylinderFlag)
					Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightSigmaNumerical + WeightPhiNumerical + ((double) (i + j))
										       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))
										       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift)) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))));
				      BMatrices[1].SetMatrixElement(this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisSigma1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
								    this->Get2RMatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisSigma2.GetNbrColumn(), this->StartingIndexPerPLevel[j][q]), Tmp);
				    } 
				}
			    }

			  if ((this->RIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + 2 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (4 * (j - i) + 4 + NValueShift) / 2;
			      N1 = N2 + QValue - 2;
			    }
			  if (this->SelfDualFlag == false)
			    {		  
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhi1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      double Tmp = 0.0;
				      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					{
					  double Tmp1 = 0.0;			      
					  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					    {
					      Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
					    }
					  Tmp += TmpOrthogonalBasisPhi1(NeutralIndex3, NeutralIndex1) * Tmp1;
					}
				      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				      if (this->CylinderFlag)
					Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightSigmaNumerical + WeightPhiNumerical + ((double) (i + j))
										       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))
										       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift)) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))));
				      BMatrices[1].SetMatrixElement(this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisSigma1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
								    this->Get2RMatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisSigma2.GetNbrColumn(), this->StartingIndexPerPLevel[j][q]), Tmp);
				    }
				}
			    }
			}
		    }
		}	      
	    }
	}
    }
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      this->RealBMatrices[i] = BMatrices[i];
    }
  delete[] BMatrices;
  delete[] ScalarProductSigma;
  delete[] ScalarProductPhi;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      delete[] MatrixPsi01[i];
      delete[] MatrixPsi10[i];
      delete U1BosonBasis[i];
      delete[] RationalMultiplicityFactor[i];
      delete[] MultiplicityFactor[i];
    }
  delete[] TmpNbrElementPerRow;
  delete[] U1BosonBasis;
  delete[] MatrixPsi01;
  delete[] MatrixPsi10;
  delete[] OrthogonalBasisSigmaLeft;
  delete[] OrthogonalBasisPhiLeft;
  delete[] OrthogonalBasisSigmaRight;
  delete[] OrthogonalBasisPhiRight;
  delete[] RationalMultiplicityFactor;
  delete[] MultiplicityFactor;
}


// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSClustered2RQuasiholeSectorMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int MinQ;
  int MaxQ;
  this->GetChargeIndexRange(0, MinQ, MaxQ);
  if (padding == true)
    {
      if ((this->RIndex & 1) == 0)
	rowIndex = this->PLevel + (this->RIndex / 2) - MinQ;
      else
	rowIndex = 2 * this->PLevel + this->RIndex - 1 - MinQ;      
      columnIndex = rowIndex;
    }
  else
    {
      if ((this->RIndex & 1) == 0)
	{
	  rowIndex = this->PLevel + this->RIndex - MinQ;
	  columnIndex = this->PLevel + this->RIndex - 1 - MinQ ;
	}
      else
	{
	  rowIndex = 2 * (this->PLevel + this->RIndex) - MinQ;
	  columnIndex = 2 * (this->PLevel + this->RIndex - 2) - MinQ;
	}
    }
}

