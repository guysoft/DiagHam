////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the clustered (k=2,r) states            //
//                                                                            //
//                        last modification : 30/10/2012                      //
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
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
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

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix()
{
  this->UseRationalFlag = true;
}

// constructor 
//
// rindex = r index (i.e. clustered (k=2,r) states) 
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// useRational = use arbitrary precision numbers for all the CFT calculations
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, bool useRational,
						   bool trimChargeIndices, bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->UseRationalFlag = useRational;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->CentralCharge = LongRational((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SquareMatrixElementNormalization = LongRational(1, 1);
  this->MatrixElementNormalization = 1.0;
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
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
// useRational = use arbitrary precision numbers for all the CFT calculations
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation
  
FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool useRational, 
						   bool trimChargeIndices, bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->UseRationalFlag = useRational;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->Kappa = kappa;
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->CentralCharge = LongRational((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SquareMatrixElementNormalization = LongRational(1, 1);
  this->MatrixElementNormalization = 1.0;
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
  this->CreateBMatrices(cftDirectory, architecture);
}


// constructor from a file describing the state
//
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// fileName = name of the file that contains the state description
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int pLevel, int nbrBMatrices, char* fileName, 
						   bool trimChargeIndices, bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->UniformChargeIndexRange = !trimChargeIndices;
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
      ErrorFlag = StateDefinition.GetAsSingleLongRational("CentralCharge", this->CentralCharge);
      if (StateDefinition["UseRational"] != 0)	
	this->UseRationalFlag = false;
      else
	ErrorFlag = StateDefinition.GetAsBoolean("UseRational", this->UseRationalFlag);
      if (StateDefinition["EMatrixDegeneracy"] != 0)
	{
	  ErrorFlag = StateDefinition.GetAsSingleInteger("EMatrixDegeneracy", this->TransferMatrixDegeneracy);	  
	}
      else
	{
	  this->TransferMatrixDegeneracy = this->RIndex + 2;
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
      if (StateDefinition["PsiSquareMatrixElement"] != 0)
	{
	  ErrorFlag = StateDefinition.GetAsSingleLongRational("PsiSquareMatrixElement", this->SquareMatrixElementNormalization);
	}
      else
	{
	  this->SquareMatrixElementNormalization = LongRational(1, 1);
	}
      this->MatrixElementNormalization = sqrt(fabs(this->SquareMatrixElementNormalization.GetNumericalValue()));
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
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, 
						   bool trimChargeIndices, bool cylinderFlag, double kappa)
{
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->LoadMatrices(fileName);
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->CentralCharge = LongRational((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SquareMatrixElementNormalization = LongRational(1, 1);
  this->MatrixElementNormalization = 1.0;
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
}

// destructor
//

FQHEMPSClustered2RMatrix::~FQHEMPSClustered2RMatrix()
{
  delete[] this->BMatrixOutputName;
}
  
// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

void FQHEMPSClustered2RMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  numerator = 2;
  denominator = 2 + this->RIndex;
}

// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSClustered2RMatrix::CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture)
{
  LongRational CentralCharge12 (this->CentralCharge);
  cout << "central charge = " << CentralCharge12 << endl;
  CentralCharge12 /= 12l;

  double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();
  double WeightPrimaryFieldMatrixElementNumerical = this->WeightPrimaryFieldMatrixElement.GetNumericalValue();

  double WeightIdentityNumerical = this->WeightIdentity.GetNumericalValue();
  double WeightPsiNumerical = this->WeightPsi.GetNumericalValue();
  long* Partition = new long[2 * (this->PLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevel + 2];
  this->TemporaryOccupationNumber = new unsigned long [this->PLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductIdentity = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductPsi = new LongRationalMatrix[this->PLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi01 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi10 = new LongRationalMatrix*[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[this->PLevel + 1];

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

  cout << "weight: " <<   this->WeightIdentity << " " << this->WeightPsi << endl;
  char* TmpScalarProductIdentityFileName = 0; 
  char* TmpScalarProductPsiFileName = 0;
  if (cftDirectory != 0)
    {
      TmpScalarProductIdentityFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductPsiFileName = new char[512 + strlen(cftDirectory)];
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      if (cftDirectory != 0)
	{
	  if (this->UseRationalFlag == true)
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_scalarproducts_identity_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	      sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_scalarproducts_psi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	    }
	  else
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_num_scalarproducts_identity_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	      sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_num_scalarproducts_psi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	    }
	}
      if ((cftDirectory != 0) && (IsFile(TmpScalarProductIdentityFileName)))
	{		
	  if (this->UseRationalFlag == true)
	    {
	      RationalScalarProductIdentity[i].ReadMatrix(TmpScalarProductIdentityFileName);
	    }
	  else
	    {
	      ScalarProductIdentity[i].ReadMatrix(TmpScalarProductIdentityFileName);
	    }
	}
      else
	{
	  if (this->UseRationalFlag == true)
	    {
	      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12, 
						     this->WeightIdentity,
						     RationalScalarProductIdentity,  i- 1);
	      Operation1.ApplyOperation(architecture);
	      RationalScalarProductIdentity[i] = Operation1.GetRationalMatrixElements();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  RationalScalarProductIdentity[i].WriteMatrix(TmpScalarProductIdentityFileName);
		}
	    }
	  else
	    {
	      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12Numerical, 
						     WeightIdentityNumerical,
						     ScalarProductIdentity,  i- 1);
	      Operation1.ApplyOperation(architecture);
	      ScalarProductIdentity[i] = Operation1.GetOverlapMatrix();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  ScalarProductIdentity[i].WriteMatrix(TmpScalarProductIdentityFileName);
		}
	    }
	}
      if ((cftDirectory != 0) && (IsFile(TmpScalarProductPsiFileName)))
	{
	  if (this->UseRationalFlag == true)
	    {
	      RationalScalarProductPsi[i].ReadMatrix(TmpScalarProductPsiFileName);
	    }
	  else
	    {
	      ScalarProductPsi[i].ReadMatrix(TmpScalarProductPsiFileName);
	    }
	}
      else
	{
	  if (this->UseRationalFlag == true)
	    {
	      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12, 
						     this->WeightPsi,
						     RationalScalarProductPsi,  i - 1);
	      Operation2.ApplyOperation(architecture);
	      RationalScalarProductPsi[i] = Operation2.GetRationalMatrixElements();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  RationalScalarProductPsi[i].WriteMatrix(TmpScalarProductPsiFileName);
		}
	    }
	  else
	    {
	      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12Numerical, 
						     WeightPsiNumerical,
						     ScalarProductPsi,  i - 1);
	      Operation2.ApplyOperation(architecture);
	      ScalarProductPsi[i] = Operation2.GetOverlapMatrix();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  ScalarProductPsi[i].WriteMatrix(TmpScalarProductPsiFileName);
		}
	    }
	}
      RealSymmetricMatrix TmpMatrix;
      if (this->UseRationalFlag == true)
 	{
 	  LongRationalMatrix TmpRationalMatrix(RationalScalarProductIdentity[i].GetNbrRow(), RationalScalarProductIdentity[i].GetNbrColumn());
 	  for (int k = 0; k < RationalScalarProductIdentity[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductIdentity[i].GetNbrColumn(); ++l)
	      {
		TmpRationalMatrix[l][k] = RationalScalarProductIdentity[i][l][k] * (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  TmpMatrix = TmpRationalMatrix;
 	}
       else
	{
	  TmpMatrix = RealSymmetricMatrix (ScalarProductIdentity[i].GetNbrRow(), ScalarProductIdentity[i].GetNbrColumn());
	  for (int k = 0; k < ScalarProductIdentity[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductIdentity[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductIdentity[i].GetMatrixElement(k, l, Tmp);
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
      cout << endl;
      cout << "nbr of null vectors identity sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisIdentityLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  OrthogonalBasisIdentityRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisIdentityLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisIdentityRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisIdentityLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisIdentityRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisIdentityLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisIdentityRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisIdentityLeft[i] = RealMatrix();
	  OrthogonalBasisIdentityRight[i] = RealMatrix();
	}

       if (this->UseRationalFlag == true)
 	{
	  LongRationalMatrix TmpRationalMatrix(RationalScalarProductPsi[i].GetNbrRow(), RationalScalarProductPsi[i].GetNbrColumn());
	  for (int k = 0; k < RationalScalarProductPsi[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductPsi[i].GetNbrColumn(); ++l)
	      {
		TmpRationalMatrix[l][k] = RationalScalarProductPsi[i][l][k] * (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  TmpMatrix = TmpRationalMatrix;
 	}
       else
	{
	  TmpMatrix = RealSymmetricMatrix (ScalarProductPsi[i].GetNbrRow(), ScalarProductPsi[i].GetNbrColumn());
	  for (int k = 0; k < ScalarProductPsi[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductPsi[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductPsi[i].GetMatrixElement(k, l, Tmp);
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
	{
	  if (fabs(TmpDiag(n, n)) < Error)
	    ++Count;
	}
      cout << "nbr of null vectors Psi sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;

      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisPsiLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  OrthogonalBasisPsiRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisPsiLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisPsiRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisPsiLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisPsiRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisPsiLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisPsiRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisPsiLeft[i] = RealMatrix();
	  OrthogonalBasisPsiRight[i] = RealMatrix();
	}
      cout << "---------------------------------" << endl;
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      if (this->UseRationalFlag == true)
 	{
 	  for (int k = 0; k < RationalScalarProductIdentity[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductIdentity[i].GetNbrColumn(); ++l)
	      {
		RationalScalarProductIdentity[i][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  ScalarProductIdentity[i] = RationalScalarProductIdentity[i];
 	}
       else
	{
	  for (int k = 0; k < ScalarProductIdentity[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductIdentity[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductIdentity[i].GetMatrixElement(k, l, Tmp);
		Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		ScalarProductIdentity[i].SetMatrixElement(k, l, Tmp);
	      }
	}
      if (this->UseRationalFlag == true)
 	{
	  for (int k = 0; k < RationalScalarProductPsi[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductPsi[i].GetNbrColumn(); ++l)
	      {
		RationalScalarProductPsi[i][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  ScalarProductPsi[i] = RationalScalarProductPsi[i];
 	}
       else
	{
	  for (int k = 0; k < ScalarProductPsi[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductPsi[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductPsi[i].GetMatrixElement(k, l, Tmp);
		Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		ScalarProductPsi[i].SetMatrixElement(k, l, Tmp);
	      }
	}
    }

  this->StartingIndexPerPLevel = new int* [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->IdentityBasisDimension = new int [this->PLevel + 1];
  this->PsiBasisDimension = new int [this->PLevel + 1];
  this->U1BasisDimension = new int [this->PLevel + 1];	
  this->NeutralSectorDimension = new int* [2];
  this->NeutralSectorDimension[0] = new int [this->PLevel + 1];
  this->NeutralSectorDimension[1] = new int [this->PLevel + 1];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->IdentityBasisDimension[i] = OrthogonalBasisIdentityLeft[i].GetNbrColumn();
      this->PsiBasisDimension[i] = OrthogonalBasisPsiLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[0][i] = OrthogonalBasisIdentityLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[1][i] = OrthogonalBasisPsiLeft[i].GetNbrColumn();
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
      QValueDenominator = 1;
    }
  else
    {
      QValue = 2 + this->RIndex;
      this->NbrNValue = ((4 * this->PLevel) + QValue) + this->RIndex + 1;
      NValueShift = 4 * this->PLevel - 2;
      QValueDenominator = 2;
      ExtraCylinderFactor = 4.0;
    }

  int MatrixSize = this->ComputeLinearizedIndexArrays();

//   this->NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() + OrthogonalBasisPsiLeft[0].GetNbrColumn())) * this->NbrNValuesPerPLevel[0];
//   for (int i = 1; i <= this->PLevel; ++i)
//     {
//       this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
//       this->StartingIndexPerPLevel[i] = new int [i + 1];      
//       this->StartingIndexPerPLevel[i][0] = this->TotalStartingIndexPerPLevel[i];
//       int Tmp = 0;
//       int Tmp2;
//       for (int j = 0; j < i; ++j)
// 	{
// 	  Tmp2 = U1BosonBasis[j]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[i - j].GetNbrColumn() + OrthogonalBasisPsiLeft[i - j].GetNbrColumn()) * this->NbrNValuesPerPLevel[i];
// 	  this->StartingIndexPerPLevel[i][j + 1] = Tmp2 + this->StartingIndexPerPLevel[i][j];
// 	  Tmp += Tmp2;
// 	}
//       Tmp += U1BosonBasis[i]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() + OrthogonalBasisPsiLeft[0].GetNbrColumn()) * this->NbrNValuesPerPLevel[i];
//       this->NbrIndicesPerPLevel[i] =  Tmp;
//     }
  
//  int MatrixSize = this->NbrIndicesPerPLevel[this->PLevel] + this->TotalStartingIndexPerPLevel[this->PLevel];
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
		  sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_matrixelement_identitypsi_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		  sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_matrixelement_psiidentity_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		}
	      else
		{
		  sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_num_matrixelement_identitypsi_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		  sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_num_matrixelement_psiidentity_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		}
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpScalarProductIdentityFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi01[i][j].ReadMatrix(TmpScalarProductIdentityFileName);
		}
	      else
		{
		  MatrixPsi01[i][j].ReadMatrix(TmpScalarProductIdentityFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, j, CentralCharge12, 
							 this->WeightIdentity, this->WeightPsi, this->WeightPrimaryFieldMatrixElement,
							 RationalMatrixPsi01,  i - 1, j);
		  Operation1.ApplyOperation(architecture);
		  RationalMatrixPsi01[i][j] = Operation1.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi01[i][j].WriteMatrix(TmpScalarProductIdentityFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							 WeightIdentityNumerical, WeightPsiNumerical, WeightPrimaryFieldMatrixElementNumerical,
							 MatrixPsi01,  i - 1, j);
		  Operation1.ApplyOperation(architecture);
		  MatrixPsi01[i][j] = Operation1.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi01[i][j].WriteMatrix(TmpScalarProductIdentityFileName);
		    }
		}
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpScalarProductPsiFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi10[i][j].ReadMatrix(TmpScalarProductPsiFileName);
		}
	      else
		{
		  MatrixPsi10[i][j].ReadMatrix(TmpScalarProductPsiFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							 this->WeightPsi, this->WeightIdentity, this->WeightPrimaryFieldMatrixElement,
							 RationalMatrixPsi10,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  RationalMatrixPsi10[i][j] = Operation2.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi10[i][j].WriteMatrix(TmpScalarProductPsiFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							 WeightPsiNumerical, WeightIdentityNumerical, WeightPrimaryFieldMatrixElementNumerical,
							 MatrixPsi10,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  MatrixPsi10[i][j] = Operation2.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi10[i][j].WriteMatrix(TmpScalarProductPsiFileName);
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
	      for (int k = 0; k < RationalMatrixPsi10[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < RationalMatrixPsi10[i][j].GetNbrColumn(); ++l)
		  {
		    RationalMatrixPsi10[i][j][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[j][l]);
		  }
	      MatrixPsi01[i][j] = RationalMatrixPsi01[i][j];
	      MatrixPsi10[i][j] = RationalMatrixPsi10[i][j];
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
	      for (int k = 0; k < MatrixPsi10[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < MatrixPsi10[i][j].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    MatrixPsi10[i][j].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[j][l]);
		    MatrixPsi10[i][j].SetMatrixElement(k, l, Tmp);
		  }
	    }
	  MatrixPsi01[i][j] *= MatrixElementNormalization;
	  MatrixPsi10[i][j] *= MatrixElementNormalization;
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
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[i - p];
	  RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductPsi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      if ((this->RIndex & 1) == 0)
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 1; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 1; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		}
	      else
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 2; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, j - 2, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 2; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, j - 2, p, ChargedIndex, NeutralIndex1)];
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
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[i - p];
	  RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductPsi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      if ((this->RIndex & 1) == 0)
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 1; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductIdentity(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentityRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisIdentityLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      if (this->CylinderFlag)
				Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
									 + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1),
 							    this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 1; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductPsi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsiRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisPsiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      if (this->CylinderFlag)
				Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
									 + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1),
 							    this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
			    }
			}
		    }
		}
	      else
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 2; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductIdentity(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentityRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisIdentityLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      if (this->CylinderFlag)
				Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
									 + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 2, p, ChargedIndex, NeutralIndex1),
 							    this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 2; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductPsi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsiRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisPsiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      if (this->CylinderFlag)
				Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
									 + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 2, p, ChargedIndex, NeutralIndex1),
 							    this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
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
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
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
			      N2 = (2 * (j - i) + this->RIndex + 1 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (4 * (j - i) + 2 * this->RIndex + 2 + NValueShift) / 2;
			      N1 = N2 + QValue - 2;
			    }			  
			  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
			      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
			    { 
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1)];
				    }
				}

			    }
			  if ((this->RIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + 1 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (4 * (j - i) + 2 + NValueShift) / 2;
			      N1 = N2 + QValue - 2;
			    }
			  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
			      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
			    { 
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1)];
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
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
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
			      N2 = (2 * (j - i) + this->RIndex + 1 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (4 * (j - i) + 2 * this->RIndex + 2 + NValueShift) / 2;
			      N1 = N2 + QValue - 2;
			    }			  
			  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
			      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
			    { 
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      double Tmp = 0.0;
				      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					{
					  double Tmp1 = 0.0;			      
					  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					    {
					  Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
					    }
					  Tmp += TmpOrthogonalBasisIdentity1(NeutralIndex3, NeutralIndex1) * Tmp1;
					}
				      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				      if (this->CylinderFlag)
					Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
										       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))
										       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift)) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))));
				      BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
								    this->Get2RMatrixIndexV2(j, 1, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
				    }
				  
				}
			    }
			  if ((this->RIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + 1 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (4 * (j - i) + 2 + NValueShift) / 2;
			      N1 = N2 + QValue - 2;
			    }
			  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
			      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
			    { 
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      double Tmp = 0.0;
				      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					{
					  double Tmp1 = 0.0;			      
					  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					    {
					      Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentity2(NeutralIndex4, NeutralIndex2);				  
					    }
					  Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
					}
				      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				      if (this->CylinderFlag)
					Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
										       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))
										       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift)) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))));
				      BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1),
								    this->Get2RMatrixIndexV2(j, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
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
  delete[] ScalarProductIdentity;
  delete[] ScalarProductPsi;
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
  delete[] OrthogonalBasisIdentityLeft;
  delete[] OrthogonalBasisPsiLeft;
  delete[] OrthogonalBasisIdentityRight;
  delete[] OrthogonalBasisPsiRight;
  delete[] RationalMultiplicityFactor;
  delete[] MultiplicityFactor;
}

// compute the scalar product matrices of the Virasoro descendant
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight = weight of the primary field that is considered
// return value = scalar product

LongRational FQHEMPSClustered2RMatrix::ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
									       LongRational& centralCharge12, LongRational& weight)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if ((position == partitionLength) || (position == 0))
    return 0l;
  if (partitionLength == 2)
    {
      if (partition[0] != -partition[1])
	{
	  return 0l;
	}
      else
	{
	  LongRational Tmp1 (centralCharge12);
	  Tmp1 *= partition[0] * (partition[0] * partition[0] - 1l);
	  LongRational Tmp2 (weight);
	  Tmp2 *= 2l * partition[0];
	  Tmp1 += Tmp2;
	  return Tmp1;
	}
    }
  LongRational Tmp(0l);
  if ((partition[position - 1] + partition[position]) == 0)
    {
      long TmpLength = 0l;
      long Store = partition[position - 1];
      for (int i = position + 1; i < partitionLength; ++i)
	TmpLength += partition[i];
      for (int i = position + 1; i < partitionLength; ++i)
	partition[i - 2] = partition[i];
      Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
	       + (2l * Store) * (weight - TmpLength)) * 
	      this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, weight));
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 2];
      partition[position - 1] = Store;
      partition[position] = -Store;
    }
  else
    {
      long Store1 = partition[position - 1];
      long Store2 = partition[position];
      partition[position - 1] += partition[position];
      for (int i = position + 1; i < partitionLength; ++i)
	partition[i - 1] = partition[i];
      if ((Store1 + Store2) > 0)
	{
	  Tmp += ((Store1 - Store2) 
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position, centralCharge12, weight));
	}
      else
	{
	  Tmp += ((Store1 - Store2) 
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position - 1, centralCharge12, weight));
	}
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 1];
      partition[position] = Store2;
      partition[position - 1] = Store1;
    }

  long Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  Tmp += this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, weight);
  Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  return Tmp;
}

// compute the scalar product matrices of the Virasoro descendant, using information from previous levels
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight = weight of the primary field that is considered
// precomputedScalarProduct = matrices where scalar product matrix elements computed for previous levels are stored
// precomputedScalarProductMaxPLevel = maxixum P level that can be accessed through precomputedScalarProduct
// basis = basis that related the partitions to their index
// temporaryOccupationNumber = local temporary to store the occupation numbers 
// return value = scalar product

LongRational FQHEMPSClustered2RMatrix::ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
									       LongRational& centralCharge12, LongRational& weight,
									       LongRationalMatrix* precomputedScalarProduct, int precomputedScalarProductMaxPLevel, 
									       BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if ((position == partitionLength) || (position == 0))
    return 0l;
  if (partitionLength == 2)
    {
      if (partition[0] != -partition[1])
	{
	  return 0l;
	}
      else
	{
	  LongRational Tmp1 (centralCharge12);
	  Tmp1 *= partition[0] * (partition[0] * partition[0] - 1l);
	  LongRational Tmp2 (weight);
	  Tmp2 *= 2l * partition[0];
	  Tmp1 += Tmp2;
	  return Tmp1;
	}
    }
  int TmpPosition = 0;
  while ((TmpPosition < position) && (partition[TmpPosition] >= 0))
    ++TmpPosition;
  if (TmpPosition == position)
    {
      TmpPosition = 1;
      int TmpPLevel1 = partition[0];
      bool FlagSorted = true;
      while ((TmpPosition < position) && (FlagSorted == true))
	{
	  TmpPLevel1 += partition[TmpPosition];
	  if (partition[TmpPosition - 1] < partition[TmpPosition])
	    FlagSorted = false;
	  ++TmpPosition;
	}
      if ((TmpPLevel1 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
	{
	  int TmpPLevel2 = -partition[TmpPosition];	  
	  FlagSorted = true;
	  ++TmpPosition;
	  while (TmpPosition < partitionLength)
	    {
	      TmpPLevel2 -= partition[TmpPosition];
	      if (partition[TmpPosition - 1] < partition[TmpPosition])
		FlagSorted = false;
	      ++TmpPosition;
	    }
	  if ((TmpPLevel2 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
	    {
	      for (int k = 0; k <= (this->PLevel + 1); ++k)
		temporaryOccupationNumber[k] = 0x0ul;	  
	      for (TmpPosition = 0; TmpPosition < position; ++TmpPosition)
		{
		  temporaryOccupationNumber[partition[TmpPosition]]++;	      
		}
	      temporaryOccupationNumber[0] = TmpPLevel1 - position;
	      int TmpIndex1 = basis[TmpPLevel1]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
	      for (int k = 0; k <= (this->PLevel + 1); ++k)
		temporaryOccupationNumber[k] = 0x0ul;	  
	      for (TmpPosition = position; TmpPosition < partitionLength; ++TmpPosition)
		{
		  temporaryOccupationNumber[-partition[TmpPosition]]++;	      
		}
	      temporaryOccupationNumber[0] = TmpPLevel2 - partitionLength + position;
	      int TmpIndex2 = basis[TmpPLevel2]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
	      LongRational Tmp;
	      precomputedScalarProduct[TmpPLevel1].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
	      return Tmp;
	    }
	}
    }
  LongRational Tmp(0l);
  if ((partition[position - 1] + partition[position]) == 0)
    {
      long TmpLength = 0l;
      long Store = partition[position - 1];
      for (int i = position + 1; i < partitionLength; ++i)
	TmpLength += partition[i];
      for (int i = position + 1; i < partitionLength; ++i)
	partition[i - 2] = partition[i];
      Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
	       + (2l * Store) * (weight - TmpLength)) * 
	      this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, weight,
							   precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 2];
      partition[position - 1] = Store;
      partition[position] = -Store;
    }
  else
    {
      long Store1 = partition[position - 1];
      long Store2 = partition[position];
      partition[position - 1] += partition[position];
      for (int i = position + 1; i < partitionLength; ++i)
	partition[i - 1] = partition[i];
      if ((Store1 + Store2) > 0)
	{
	  Tmp += ((Store1 - Store2) 
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position, centralCharge12, weight,
								 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	}
      else
	{
	  Tmp += ((Store1 - Store2) 
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position - 1, centralCharge12, weight,
								 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	}
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 1];
      partition[position] = Store2;
      partition[position - 1] = Store1;
    }

  long Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  Tmp += this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, weight,
						      precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber);
  Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  return Tmp;
}

// compute the scalar product matrices of the Virasoro descendant, using information from previous levels
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight = weight of the primary field that is considered
// precomputedScalarProduct = matrices where scalar product matrix elements computed for previous levels are stored
// precomputedScalarProductMaxPLevel = maxixum P level that can be accessed through precomputedScalarProduct
// basis = basis that related the partitions to their index
// temporaryOccupationNumber = local temporary to store the occupation numbers 
// return value = scalar product

double FQHEMPSClustered2RMatrix::ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
									 double& centralCharge12, double& weight,
									 RealSymmetricMatrix* precomputedScalarProduct, int precomputedScalarProductMaxPLevel, 
									 BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if ((position == partitionLength) || (position == 0))
    return 0.0;
  if (partitionLength == 2)
    {
      if (partition[0] != -partition[1])
	{
	  return 0.0;
	}
      else
	{
	  double Tmp1 = centralCharge12;
	  Tmp1 *= partition[0] * (partition[0] * partition[0] - 1l);
	  double Tmp2 = weight;
	  Tmp2 *= 2l * partition[0];
	  Tmp1 += Tmp2;
	  return Tmp1;
	}
    }
  int TmpPosition = 0;
  while ((TmpPosition < position) && (partition[TmpPosition] >= 0))
    ++TmpPosition;
  if (TmpPosition == position)
    {
      TmpPosition = 1;
      int TmpPLevel1 = partition[0];
      bool FlagSorted = true;
      while ((TmpPosition < position) && (FlagSorted == true))
	{
	  TmpPLevel1 += partition[TmpPosition];
	  if (partition[TmpPosition - 1] < partition[TmpPosition])
	    FlagSorted = false;
	  ++TmpPosition;
	}
      if ((TmpPLevel1 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
	{
	  int TmpPLevel2 = -partition[TmpPosition];	  
	  FlagSorted = true;
	  ++TmpPosition;
	  while (TmpPosition < partitionLength)
	    {
	      TmpPLevel2 -= partition[TmpPosition];
	      if (partition[TmpPosition - 1] < partition[TmpPosition])
		FlagSorted = false;
	      ++TmpPosition;
	    }
	  if ((TmpPLevel2 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
	    {
	      for (int k = 0; k <= (this->PLevel + 1); ++k)
		temporaryOccupationNumber[k] = 0x0ul;	  
	      for (TmpPosition = 0; TmpPosition < position; ++TmpPosition)
		{
		  temporaryOccupationNumber[partition[TmpPosition]]++;	      
		}
	      temporaryOccupationNumber[0] = TmpPLevel1 - position;
	      int TmpIndex1 = basis[TmpPLevel1]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
	      for (int k = 0; k <= (this->PLevel + 1); ++k)
		temporaryOccupationNumber[k] = 0x0ul;	  
	      for (TmpPosition = position; TmpPosition < partitionLength; ++TmpPosition)
		{
		  temporaryOccupationNumber[-partition[TmpPosition]]++;	      
		}
	      temporaryOccupationNumber[0] = TmpPLevel2 - partitionLength + position;
	      int TmpIndex2 = basis[TmpPLevel2]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
	      double Tmp;
	      precomputedScalarProduct[TmpPLevel1].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
	      return Tmp;
	    }
	}
    }
  double Tmp = 0l;
  if ((partition[position - 1] + partition[position]) == 0)
    {
      long TmpLength = 0l;
      long Store = partition[position - 1];
      for (int i = position + 1; i < partitionLength; ++i)
	TmpLength += partition[i];
      for (int i = position + 1; i < partitionLength; ++i)
	partition[i - 2] = partition[i];
      Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
	       + (2l * Store) * (weight - TmpLength)) * 
	      this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, weight,
							   precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 2];
      partition[position - 1] = Store;
      partition[position] = -Store;
    }
  else
    {
      long Store1 = partition[position - 1];
      long Store2 = partition[position];
      partition[position - 1] += partition[position];
      for (int i = position + 1; i < partitionLength; ++i)
	partition[i - 1] = partition[i];
      if ((Store1 + Store2) > 0)
	{
	  Tmp += ((Store1 - Store2) 
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position, centralCharge12, weight,
								 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	}
      else
	{
	  Tmp += ((Store1 - Store2) 
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position - 1, centralCharge12, weight,
								 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	}
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 1];
      partition[position] = Store2;
      partition[position - 1] = Store1;
    }

  long Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  Tmp += this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, weight,
						      precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber);
  Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  return Tmp;
}

// compute the matrix elements of any primary field in the Virasoro descendant basis
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// descendantPosition = location of the primary field
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight1 = weight of the primary field that is considered for the left state
// weight2 = weight of the primary field that is considered for the right state
// weight = weight of the primary field whose matrix elements are computed
// return value = matrix element
  
LongRational FQHEMPSClustered2RMatrix::ComputeDescendantMatrixElement (long* partition, int partitionLength, 
								       int descendantPosition, int position, 
								       LongRational& centralCharge12, LongRational& weight1, 
								       LongRational& weight2, LongRational& weight)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  if ((descendantPosition < partitionLength) && (partition[partitionLength - 1] > 0))
    {
      return 0l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if (descendantPosition == partitionLength) 
    {
      LongRational Tmp(1l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = partitionLength - 1; i >= 0; --i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = 0; j < i; ++j)
	    Tmp3 += partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 += Tmp3;	  
	  Tmp *= Tmp2;
	}
      return Tmp;
    }
  if (position == 0)
    {
      LongRational Tmp(1l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = 0; i < partitionLength ; ++i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = i + 1; j < partitionLength; ++j)
	    Tmp3 -= partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 -= Tmp3;
	  Tmp *= Tmp2;
	}
      Tmp *= 1l - (2l * (partitionLength & 1l));
      return Tmp;
    }
  if (descendantPosition < position)
    {
      LongRational Tmp(0l);
      if ((partition[position - 1] + partition[position]) == 0)
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
		   + (2l * Store) * (weight2 - TmpLength)) * 
		  this->ComputeDescendantMatrixElement(partition, partitionLength - 2, descendantPosition, 
						       position - 1, centralCharge12, weight1, weight2, weight));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
      else
	{
 	  long Store1 = partition[position - 1];
	  long Store2 = partition[position];
	  partition[position - 1] += partition[position];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 1] = partition[i];
	  if ((Store1 + Store2) > 0)
	    {
	      Tmp += ((Store1 - Store2) 
		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
							     position, centralCharge12, weight1, weight2, weight));
	    }
	  else
	    {
	      Tmp += ((Store1 - Store2) 
		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
							     position - 1, centralCharge12, weight1, weight2, weight));
	    }
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 1];
	  partition[position] = Store2;
	  partition[position - 1] = Store1;
	}
      
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp += this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition, position + 1, 
						  centralCharge12, weight1, weight2, weight);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;  
    }

  LongRational Tmp1 = this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition - 1, position, 
							   centralCharge12, weight1, weight2, weight);
  LongRational Tmp2 = weight;
  Tmp2 *= partition[position - 1];
  Tmp2 += weight1;
  Tmp2 -= weight2;
  long TmpLength = 0l;
  for (int i = position; i < partitionLength; ++i)
    TmpLength += partition[i];
  for (int i = 0; i < (position - 1); ++i)
    TmpLength += partition[i];
  Tmp2 += TmpLength;
  long Store = partition[position - 1];
  for (int i = position; i < partitionLength; ++i)
    partition[i - 1] = partition[i];
  Tmp2 *= this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition - 1, position - 1, 
					       centralCharge12, weight1, weight2, weight);
  for (int i = partitionLength - 1; i >= position; --i)
    partition[i] = partition[i - 1];
  partition[position - 1] = Store;
  Tmp1 += Tmp2;
  return Tmp1;
}

// compute the matrix elements of any primary field in the Virasoro descendant basis
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// descendantPosition = location of the primary field
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight1 = weight of the primary field that is considered for the left state
// weight2 = weight of the primary field that is considered for the right state
// weight = weight of the primary field whose matrix elements are computed
// precomputedDescendantMatrixElement = matrices where matrix elements computed for previous levels are stored
// precomputedDescendantMatrixElementMaxLeftPLevel = maxixum P level that can be accessed through precomputedDescendantMatrixElement for the left entry
// precomputedDescendantMatrixElementMaxRightPLevel = maxixum P level that can be accessed through precomputedDescendantMatrixElement for the right entry
// basis = basis that related the partitions to their index
// temporaryOccupationNumber = local temporary to store the occupation numbers 
// return value = matrix element
  
LongRational FQHEMPSClustered2RMatrix::ComputeDescendantMatrixElement (long* partition, int partitionLength, 
								       int descendantPosition, int position,
								       LongRational& centralCharge12, LongRational& weight1, 
								       LongRational& weight2, LongRational& weight,
								       LongRationalMatrix** precomputedDescendantMatrixElement, 
								       int precomputedDescendantMatrixElementMaxLeftPLevel, 
								       int precomputedDescendantMatrixElementMaxRightPLevel, 
								       BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  if ((descendantPosition < partitionLength) && (partition[partitionLength - 1] > 0))
    {
      return 0l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if (descendantPosition == partitionLength) 
    {
      LongRational Tmp(1l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = partitionLength - 1; i >= 0; --i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = 0; j < i; ++j)
	    Tmp3 += partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 += Tmp3;	  
	  Tmp *= Tmp2;
	}
      return Tmp;
    }
  if (position == 0)
    {
      LongRational Tmp(1l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = 0; i < partitionLength ; ++i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = i + 1; j < partitionLength; ++j)
	    Tmp3 -= partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 -= Tmp3;
	  Tmp *= Tmp2;
	}
      Tmp *= 1l - (2l * (partitionLength & 1l));
      return Tmp;
    }
  if (descendantPosition == position)
    {
      int TmpPosition = 0;
      while ((TmpPosition < position) && (partition[TmpPosition] >= 0))
	++TmpPosition;
      if (TmpPosition == position)
	{
	  TmpPosition = 1;
	  int TmpPLevel1 = partition[0];
	  bool FlagSorted = true;
	  while (TmpPosition < position)
	    {
	      TmpPLevel1 += partition[TmpPosition];
	      if (partition[TmpPosition - 1] < partition[TmpPosition])
		FlagSorted = false;
	      ++TmpPosition;
	    }
	  if ((TmpPLevel1 <= precomputedDescendantMatrixElementMaxLeftPLevel) && (FlagSorted == true))
	    {
	      int TmpPLevel2 = -partition[TmpPosition];	  
	      FlagSorted = true;
	      ++TmpPosition;
	      while (TmpPosition < partitionLength)
		{
		  TmpPLevel2 -= partition[TmpPosition];
		  if (partition[TmpPosition - 1] < partition[TmpPosition])
		    FlagSorted = false;
		  ++TmpPosition;
		}
	      if ((TmpPLevel2 <= precomputedDescendantMatrixElementMaxRightPLevel) && (FlagSorted == true))
		{
		  for (int k = 0; k <= (this->PLevel + 1); ++k)
		    temporaryOccupationNumber[k] = 0x0ul;	  
		  for (TmpPosition = 0; TmpPosition < position; ++TmpPosition)
		    {
		      temporaryOccupationNumber[partition[TmpPosition]]++;	      
		    }
		  temporaryOccupationNumber[0] = TmpPLevel1 - position;
		  int TmpIndex1 = basis[TmpPLevel1]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
		  for (int k = 0; k <= (this->PLevel + 1); ++k)
		    temporaryOccupationNumber[k] = 0x0ul;	  
		  for (TmpPosition = position; TmpPosition < partitionLength; ++TmpPosition)
		    {
		      temporaryOccupationNumber[-partition[TmpPosition]]++;	      
		    }
		  temporaryOccupationNumber[0] = TmpPLevel2 - partitionLength + position;
		  int TmpIndex2 = basis[TmpPLevel2]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
		  LongRational Tmp;
		  precomputedDescendantMatrixElement[TmpPLevel1][TmpPLevel2].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
		  return Tmp;
		}
	    }
	}
    }
  if (descendantPosition < position)
    {
      LongRational Tmp(0l);
      if ((partition[position - 1] + partition[position]) == 0)
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
		   + (2l * Store) * (weight2 - TmpLength)) * 
		  this->ComputeDescendantMatrixElement(partition, partitionLength - 2, descendantPosition, 
						       position - 1, centralCharge12, weight1, weight2, weight,
						       precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
						       precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
      else
	{
 	  long Store1 = partition[position - 1];
	  long Store2 = partition[position];
	  partition[position - 1] += partition[position];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 1] = partition[i];
	  if ((Store1 + Store2) > 0)
	    {
	      Tmp += ((Store1 - Store2) 
		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
							     position, centralCharge12, weight1, weight2, weight,
							     precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
							     precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
	    }
	  else
	    {
	      Tmp += ((Store1 - Store2) 
		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
							     position - 1, centralCharge12, weight1, weight2, weight,
							     precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
							     precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
	    }
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 1];
	  partition[position] = Store2;
	  partition[position - 1] = Store1;
	}
      
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp += this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition, position + 1, 
						  centralCharge12, weight1, weight2, weight,
						  precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
						  precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;  
    }

  LongRational Tmp1 = this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition - 1, position, 
							   centralCharge12, weight1, weight2, weight,
							   precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
							   precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);


  LongRational Tmp2 = weight;
  Tmp2 *= partition[position - 1];
  Tmp2 += weight1;
  Tmp2 -= weight2;
  long TmpLength = 0l;
  for (int i = position; i < partitionLength; ++i)
    TmpLength += partition[i];
  for (int i = 0; i < (position - 1); ++i)
    TmpLength += partition[i];
  Tmp2 += TmpLength;
  long Store = partition[position - 1];
  for (int i = position; i < partitionLength; ++i)
    partition[i - 1] = partition[i];
  Tmp2 *= this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition - 1, position - 1, 
					       centralCharge12, weight1, weight2, weight,
					       precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
					       precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);
  for (int i = partitionLength - 1; i >= position; --i)
    partition[i] = partition[i - 1];
  partition[position - 1] = Store;
  Tmp1 += Tmp2;
  return Tmp1;
}

// compute the matrix elements of any primary field in the Virasoro descendant basis, using double numbers instead of long rational
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// descendantPosition = location of the primary field
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight1 = weight of the primary field that is considered for the left state
// weight2 = weight of the primary field that is considered for the right state
// weight = weight of the primary field whose matrix elements are computed
// precomputedDescendantMatrixElement = matrices where matrix elements computed for previous levels are stored
// precomputedDescendantMatrixElementMaxLeftPLevel = maxixum P level that can be accessed through precomputedDescendantMatrixElement for the left entry
// precomputedDescendantMatrixElementMaxRightPLevel = maxixum P level that can be accessed through precomputedDescendantMatrixElement for the right entry
// basis = basis that related the partitions to their index
// temporaryOccupationNumber = local temporary to store the occupation numbers 
// return value = matrix element
  
double FQHEMPSClustered2RMatrix::ComputeDescendantMatrixElement (long* partition, int partitionLength, 
								 int descendantPosition, int position,
								 double& centralCharge12, double& weight1, 
								 double& weight2, double& weight,
								 RealMatrix** precomputedDescendantMatrixElement, 
								 int precomputedDescendantMatrixElementMaxLeftPLevel, 
								 int precomputedDescendantMatrixElementMaxRightPLevel, 
								 BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber)
{
  if (partitionLength == 0)
    {
      return 1.0;
    }
  if ((descendantPosition < partitionLength) && (partition[partitionLength - 1] > 0))
    {
      return 0.0;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if (descendantPosition == partitionLength) 
    {
      double Tmp = 1.0;
      double TmpSum = weight1;
      TmpSum -= weight2;
      double Tmp2;
      for (int i = partitionLength - 1; i >= 0; --i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = 0; j < i; ++j)
	    Tmp3 += partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 += Tmp3;	  
	  Tmp *= Tmp2;
	}
      return Tmp;
    }
  if (position == 0)
    {
      double Tmp = 1.0;
      double TmpSum = weight1;
      TmpSum -= weight2;
      double Tmp2;
      for (int i = 0; i < partitionLength ; ++i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = i + 1; j < partitionLength; ++j)
	    Tmp3 -= partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 -= Tmp3;
	  Tmp *= Tmp2;
	}
      Tmp *= 1.0 - (2.0 * (partitionLength & 1l));
      return Tmp;
    }
  if (descendantPosition == position)
    {
      int TmpPosition = 0;
      while ((TmpPosition < position) && (partition[TmpPosition] >= 0))
	++TmpPosition;
      if (TmpPosition == position)
	{
	  TmpPosition = 1;
	  int TmpPLevel1 = partition[0];
	  bool FlagSorted = true;
	  while (TmpPosition < position)
	    {
	      TmpPLevel1 += partition[TmpPosition];
	      if (partition[TmpPosition - 1] < partition[TmpPosition])
		FlagSorted = false;
	      ++TmpPosition;
	    }
	  if ((TmpPLevel1 <= precomputedDescendantMatrixElementMaxLeftPLevel) && (FlagSorted == true))
	    {
	      int TmpPLevel2 = -partition[TmpPosition];	  
	      FlagSorted = true;
	      ++TmpPosition;
	      while (TmpPosition < partitionLength)
		{
		  TmpPLevel2 -= partition[TmpPosition];
		  if (partition[TmpPosition - 1] < partition[TmpPosition])
		    FlagSorted = false;
		  ++TmpPosition;
		}
	      if ((TmpPLevel2 <= precomputedDescendantMatrixElementMaxRightPLevel) && (FlagSorted == true))
		{
		  for (int k = 0; k <= (this->PLevel + 1); ++k)
		    temporaryOccupationNumber[k] = 0x0ul;	  
		  for (TmpPosition = 0; TmpPosition < position; ++TmpPosition)
		    {
		      temporaryOccupationNumber[partition[TmpPosition]]++;	      
		    }
		  temporaryOccupationNumber[0] = TmpPLevel1 - position;
		  int TmpIndex1 = basis[TmpPLevel1]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
		  for (int k = 0; k <= (this->PLevel + 1); ++k)
		    temporaryOccupationNumber[k] = 0x0ul;	  
		  for (TmpPosition = position; TmpPosition < partitionLength; ++TmpPosition)
		    {
		      temporaryOccupationNumber[-partition[TmpPosition]]++;	      
		    }
		  temporaryOccupationNumber[0] = TmpPLevel2 - partitionLength + position;
		  int TmpIndex2 = basis[TmpPLevel2]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
		  double Tmp;
		  precomputedDescendantMatrixElement[TmpPLevel1][TmpPLevel2].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
		  return Tmp;
		}
	    }
	}
    }
  if (descendantPosition < position)
    {
      double Tmp = 0.0;
      if ((partition[position - 1] + partition[position]) == 0)
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
		   + (2l * Store) * (weight2 - TmpLength)) * 
		  this->ComputeDescendantMatrixElement(partition, partitionLength - 2, descendantPosition, 
						       position - 1, centralCharge12, weight1, weight2, weight,
						       precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
						       precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
      else
	{
 	  long Store1 = partition[position - 1];
	  long Store2 = partition[position];
	  partition[position - 1] += partition[position];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 1] = partition[i];
	  if ((Store1 + Store2) > 0)
	    {
	      Tmp += ((Store1 - Store2) 
		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
							     position, centralCharge12, weight1, weight2, weight,
							     precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
							     precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
	    }
	  else
	    {
	      Tmp += ((Store1 - Store2) 
		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
							     position - 1, centralCharge12, weight1, weight2, weight,
							     precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
							     precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
	    }
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 1];
	  partition[position] = Store2;
	  partition[position - 1] = Store1;
	}
      
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp += this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition, position + 1, 
						  centralCharge12, weight1, weight2, weight,
						  precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
						  precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;  
    }

  double Tmp1 = this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition - 1, position, 
							   centralCharge12, weight1, weight2, weight,
							   precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
							   precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);


  double Tmp2 = weight;
  Tmp2 *= partition[position - 1];
  Tmp2 += weight1;
  Tmp2 -= weight2;
  long TmpLength = 0l;
  for (int i = position; i < partitionLength; ++i)
    TmpLength += partition[i];
  for (int i = 0; i < (position - 1); ++i)
    TmpLength += partition[i];
  Tmp2 += TmpLength;
  long Store = partition[position - 1];
  for (int i = position; i < partitionLength; ++i)
    partition[i - 1] = partition[i];
  Tmp2 *= this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition - 1, position - 1, 
					       centralCharge12, weight1, weight2, weight,
					       precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
					       precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);
  for (int i = partitionLength - 1; i >= position; --i)
    partition[i] = partition[i - 1];
  partition[position - 1] = Store;
  Tmp1 += Tmp2;
  return Tmp1;
}

// compute the various arrays required to convert from quantum numbers and local indices to a global linearized index
//
// return value = dimension of the B matrix

int FQHEMPSClustered2RMatrix::ComputeLinearizedIndexArrays()
{
  this->NbrNValuesPerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NInitialValuePerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NLastValuePerPLevelCFTSector = new int* [this->PLevel + 1];     
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->NbrNValuesPerPLevelCFTSector[i] = new int [this->NbrCFTSectors];
      this->NInitialValuePerPLevelCFTSector[i] = new int [this->NbrCFTSectors];
      this->NLastValuePerPLevelCFTSector[i] = new int [this->NbrCFTSectors];
      for (int l = 0; l < this->NbrCFTSectors; ++l)
	{
	  this->ComputeChargeIndexRange(i, l, this->NInitialValuePerPLevelCFTSector[i][l], this->NLastValuePerPLevelCFTSector[i][l]);
	  this->NbrNValuesPerPLevelCFTSector[i][l] =  this->NLastValuePerPLevelCFTSector[i][l] - this->NInitialValuePerPLevelCFTSector[i][l] + 1;
	}
    }

  this->StartingIndexPerPLevelCFTSectorQValue = new int** [this->PLevel + 1];
  this->NbrIndexPerPLevelCFTSectorQValue = new int** [this->PLevel + 1];
  this->StartingIndexPerPLevelCFTSectorQValueU1Sector = new int*** [this->PLevel + 1];
  this->NbrIndexPerPLevelCFTSectorQValueU1Sector = new int*** [this->PLevel + 1];
  int TotalIndex = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->StartingIndexPerPLevelCFTSectorQValue[i] = new int* [this->NbrCFTSectors];
      this->NbrIndexPerPLevelCFTSectorQValue[i] = new int* [this->NbrCFTSectors];
      this->StartingIndexPerPLevelCFTSectorQValueU1Sector[i] = new int** [this->NbrCFTSectors];
      this->NbrIndexPerPLevelCFTSectorQValueU1Sector[i] = new int** [this->NbrCFTSectors];
      for (int l = 0; l < this->NbrCFTSectors; ++l)
	{
	  this->StartingIndexPerPLevelCFTSectorQValue[i][l] = new int [this->NbrNValuesPerPLevelCFTSector[i][l]];
	  this->NbrIndexPerPLevelCFTSectorQValue[i][l] = new int [this->NbrNValuesPerPLevelCFTSector[i][l]];
	  this->StartingIndexPerPLevelCFTSectorQValueU1Sector[i][l] = new int* [this->NbrNValuesPerPLevelCFTSector[i][l]];
	  this->NbrIndexPerPLevelCFTSectorQValueU1Sector[i][l] = new int* [this->NbrNValuesPerPLevelCFTSector[i][l]];
	  for (int j = this->NInitialValuePerPLevelCFTSector[i][l];  j <= this->NLastValuePerPLevelCFTSector[i][l]; ++j)
	    {
	      this->StartingIndexPerPLevelCFTSectorQValue[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]] = TotalIndex;
	      this->StartingIndexPerPLevelCFTSectorQValueU1Sector[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]] = new int [i + 1];
	      this->NbrIndexPerPLevelCFTSectorQValueU1Sector[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]] = new int [i + 1];
	      for (int k = 0; k <= i; ++k)
		{
		  this->StartingIndexPerPLevelCFTSectorQValueU1Sector[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]][k] = TotalIndex;
		  int Tmp = this->U1BasisDimension[k] * this->NeutralSectorDimension[l][i - k];
		  this->NbrIndexPerPLevelCFTSectorQValueU1Sector[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]][k] = Tmp;
		  TotalIndex += Tmp;	      
		}
	      this->NbrIndexPerPLevelCFTSectorQValue[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]] = TotalIndex - this->StartingIndexPerPLevelCFTSectorQValue[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]];
	    }
	}
    }
  return TotalIndex;
}

// get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
//
// cftSector = index of the CFT sector
// return value = Q sector shift

int FQHEMPSClustered2RMatrix::GetQValueCFTSectorShift(int cftSector)
{
  if ((this->RIndex & 1) == 1)
    return 0;
  if (cftSector == 0)
    return 0;
  return ((this->RIndex + 2) / 2);
}

// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSClustered2RMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  if ((pLevel < 0) || (pLevel > this->PLevel))
    return 0;
  if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][0]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][0]))
    {
      int Tmp = this->NbrIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]];
      if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][1]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][1]))
	Tmp += this->NbrIndexPerPLevelCFTSectorQValue[pLevel][1][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][1]];
      return Tmp;
    }
  if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][1]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][1]))
    return this->NbrIndexPerPLevelCFTSectorQValue[pLevel][1][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][1]];
  return 0;
}

// get the range for the bond index when fixing the tuncation level, charge and CFT sector index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = range for the bond index with fixed tuncation level, charge and CFT sector index

int FQHEMPSClustered2RMatrix::GetBondIndexRange(int pLevel, int qValue, int cftSector)
{
  if ((pLevel < 0) || (pLevel > this->PLevel) || (qValue < this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]) || 
      (qValue > this->NLastValuePerPLevelCFTSector[pLevel][cftSector]) || (cftSector > 1) ||  (cftSector < 0))
    return 0;
  return this->NbrIndexPerPLevelCFTSectorQValue[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]];
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int FQHEMPSClustered2RMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][0]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][0]))
    {
      if (localIndex < this->NbrIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]])
	{	  
	  return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]]  + localIndex);
	}
      else
	{
	  return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][1][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][1]]  
		  + (localIndex - this->NbrIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]]));
	}
    }
  else
    {
      return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][1][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][1]] + localIndex);
    }
}

// get the bond index for a fixed truncation level, charge and CFT sector index
//
// localIndex = bond index in the pLevel and qValue and cftSector restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = bond index in the full bond index range

int FQHEMPSClustered2RMatrix::GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector)
{
  return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]]  + localIndex);
}


// get the charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSClustered2RMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][0];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][0];
  if (this->NInitialValuePerPLevelCFTSector[pLevel][1] < minQ)
    minQ = this->NInitialValuePerPLevelCFTSector[pLevel][1];
  if (this->NLastValuePerPLevelCFTSector[pLevel][1] > maxQ)
    maxQ = this->NLastValuePerPLevelCFTSector[pLevel][1];  
  return;
}

// get the charge index range at a given truncation level and in a given CFT sector
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSClustered2RMatrix::GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][cftSector];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][cftSector];
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSClustered2RMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
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
	  columnIndex = this->PLevel - MinQ;
	}
      else
	{
	  rowIndex = 2 * (this->PLevel + this->RIndex) - MinQ;
	  columnIndex = 2 * this->PLevel - MinQ;
	}
    }
}

// load the specific informations from the file header
// 
// file = reference on the input file stream
// return value = true if no error occurred  

bool FQHEMPSClustered2RMatrix::LoadHeader (ifstream& file)
{
    int HeaderSize = 0;
    ReadLittleEndian(file, HeaderSize);
    ReadLittleEndian(file, this->PLevel);
    ReadLittleEndian(file, this->LaughlinIndex);
    ReadLittleEndian(file, this->RIndex);
    ReadLittleEndian(file, this->NbrCFTSectors);
    ReadLittleEndian(file, this->NbrNValue);
    ReadLittleEndian(file, this->CylinderFlag);
    ReadLittleEndian(file, this->Kappa);
    ReadLittleEndian(file, this->UniformChargeIndexRange);

    this->NeutralSectorDimension = new int*[this->NbrCFTSectors];
    for (int x = 0; x < this->NbrCFTSectors; ++x)
        this->NeutralSectorDimension[x] = new int[this->PLevel + 1];

    this->U1BasisDimension = new int[this->PLevel + 1];
    this->NbrNValuesPerPLevelCFTSector = new int*[this->PLevel + 1];
    this->NInitialValuePerPLevelCFTSector = new int*[this->PLevel + 1];
    this->NLastValuePerPLevelCFTSector = new int*[this->PLevel + 1];
    this->StartingIndexPerPLevelCFTSectorQValue = new int**[this->PLevel + 1];
    this->NbrIndexPerPLevelCFTSectorQValue = new int**[this->PLevel + 1];
    this->StartingIndexPerPLevelCFTSectorQValueU1Sector = new int***[this->PLevel + 1];
    this->NbrIndexPerPLevelCFTSectorQValueU1Sector = new int***[this->PLevel + 1];
    for (int p = 0; p <= this->PLevel; ++p)
    {
        ReadLittleEndian(file, this->U1BasisDimension[p]);

        this->NbrNValuesPerPLevelCFTSector[p] = new int[this->NbrCFTSectors];
        this->NInitialValuePerPLevelCFTSector[p] = new int[this->NbrCFTSectors];
        this->NLastValuePerPLevelCFTSector[p] = new int[this->NbrCFTSectors];
        this->StartingIndexPerPLevelCFTSectorQValue[p] = new int*[this->NbrCFTSectors];
        this->NbrIndexPerPLevelCFTSectorQValue[p] = new int*[this->NbrCFTSectors];
        this->StartingIndexPerPLevelCFTSectorQValueU1Sector[p] = new int**[this->NbrCFTSectors];
        this->NbrIndexPerPLevelCFTSectorQValueU1Sector[p] = new int**[this->NbrCFTSectors];
        for (int x = 0; x < this->NbrCFTSectors; ++x)
        {
            ReadLittleEndian(file, this->NeutralSectorDimension[x][p]);

            ReadLittleEndian(file, this->NbrNValuesPerPLevelCFTSector[p][x]);
            ReadLittleEndian(file, this->NInitialValuePerPLevelCFTSector[p][x]);
            ReadLittleEndian(file, this->NLastValuePerPLevelCFTSector[p][x]);

            this->StartingIndexPerPLevelCFTSectorQValue[p][x] = new int[this->NbrNValuesPerPLevelCFTSector[p][x]];
            this->NbrIndexPerPLevelCFTSectorQValue[p][x] = new int[this->NbrNValuesPerPLevelCFTSector[p][x]];
            this->StartingIndexPerPLevelCFTSectorQValueU1Sector[p][x] = new int*[this->NbrNValuesPerPLevelCFTSector[p][x]];
            this->NbrIndexPerPLevelCFTSectorQValueU1Sector[p][x] = new int*[this->NbrNValuesPerPLevelCFTSector[p][x]];
            for (int n = 0; n < this->NbrNValuesPerPLevelCFTSector[p][x]; ++n)
            {
                ReadLittleEndian(file, this->StartingIndexPerPLevelCFTSectorQValue[p][x][n]);
                ReadLittleEndian(file, this->NbrIndexPerPLevelCFTSectorQValue[p][x][n]);
                this->StartingIndexPerPLevelCFTSectorQValueU1Sector[p][x][n] = new int[p + 1];
                this->NbrIndexPerPLevelCFTSectorQValueU1Sector[p][x][n] = new int[p + 1];
                for (int k = 0; k <= p; ++k)
                {
                    ReadLittleEndian(file, this->StartingIndexPerPLevelCFTSectorQValueU1Sector[p][x][n][k]);
                    ReadLittleEndian(file, this->NbrIndexPerPLevelCFTSectorQValueU1Sector[p][x][n][k]);
                }
            }
        }
    }

    return true;
}

// save the specific informations to the file header 
// 
// file = reference on the output file stream
// return value = true if no error occurred  

bool FQHEMPSClustered2RMatrix::SaveHeader (ofstream& file)
{
    int HeaderSize = sizeof(int) * 6 + sizeof(char) * 2 + sizeof(double);
    for (int p = 0; p <= this->PLevel; ++p)
    {
        HeaderSize += sizeof(int);
        for (int x = 0; x < this->NbrCFTSectors; ++x)
            HeaderSize += sizeof(int) * (4 + this->NbrNValuesPerPLevelCFTSector[p][x] * (2 + (p + 1) * 2));
    }

    WriteLittleEndian(file, HeaderSize);
    WriteLittleEndian(file, this->PLevel);
    WriteLittleEndian(file, this->LaughlinIndex);
    WriteLittleEndian(file, this->RIndex);
    WriteLittleEndian(file, this->NbrCFTSectors);
    WriteLittleEndian(file, this->NbrNValue);
    WriteLittleEndian(file, this->CylinderFlag);
    WriteLittleEndian(file, this->Kappa);
    WriteLittleEndian(file, this->UniformChargeIndexRange);

    for (int p = 0; p <= this->PLevel; ++p)
    {
        WriteLittleEndian(file, this->U1BasisDimension[p]);

        for (int x = 0; x < this->NbrCFTSectors; ++x)
        {
            WriteLittleEndian(file, this->NeutralSectorDimension[x][p]);

            WriteLittleEndian(file, this->NbrNValuesPerPLevelCFTSector[p][x]);
            WriteLittleEndian(file, this->NInitialValuePerPLevelCFTSector[p][x]);
            WriteLittleEndian(file, this->NLastValuePerPLevelCFTSector[p][x]);
            for (int n = 0; n < this->NbrNValuesPerPLevelCFTSector[p][x]; ++n)
            {
                WriteLittleEndian(file, this->StartingIndexPerPLevelCFTSectorQValue[p][x][n]);
                WriteLittleEndian(file, this->NbrIndexPerPLevelCFTSectorQValue[p][x][n]);
                for (int k = 0; k <= p; ++k)
                {
                    WriteLittleEndian(file, this->StartingIndexPerPLevelCFTSectorQValueU1Sector[p][x][n][k]);
                    WriteLittleEndian(file, this->NbrIndexPerPLevelCFTSectorQValueU1Sector[p][x][n][k]);
                }
            }
        }
    }

    return true;
}

// compute the charge index range at a given truncation level
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSClustered2RMatrix::ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ)
{
  if (this->UniformChargeIndexRange == true)
    {
      minQ = 0;
      maxQ = this->NbrNValue - 1;
      return;
    }
  minQ = 0;
  maxQ = this->NbrNValue - 1;
  int TmpMinQ = this->NbrNValue - 1;
  int TmpMaxQ = 0;    

  if ((this->RIndex & 1) == 0)
    {
      int NValueShift = this->PLevel;
      int QValue = 1 + (this->RIndex / 2);
      if (cftSector == 0) 
	{
	  for (int Q = 0; Q < this->NbrNValue; ++Q)
	    {
	      int QPrime = Q;
	      int TmpP = 0;
	      int TmpMaxP = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP)
		    TmpMaxP = TmpP;	    
		  QPrime -= (QValue - 1);
		  TmpP += QPrime - (this->RIndex / 2) - NValueShift;
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue - 1);
		      TmpP += QPrime - NValueShift;
		    }
		}
	      QPrime = Q;
	      TmpP = 0;
	      int TmpMaxP2 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP2)
		    TmpMaxP2 = TmpP;	    
		  TmpP -= QPrime - NValueShift;
		  QPrime += (QValue - 1);
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= QPrime  - (this->RIndex / 2) - NValueShift;
		      QPrime += (QValue - 1);
		    }
		}
	      if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		{
		  if (Q < TmpMinQ)
		    TmpMinQ = Q;
		  if (Q > TmpMaxQ)
		    TmpMaxQ = Q;	
		}    
	    }
	}
      else
	{
	  for (int Q = 0; Q < this->NbrNValue; ++Q)
	    {
	      int QPrime = Q;
	      int TmpP = 0;
	      int TmpMaxP3 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP3)
		    TmpMaxP3 = TmpP;	    
		  QPrime -= (QValue - 1);
		  TmpP += QPrime - NValueShift;
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP3)
			TmpMaxP3 = TmpP;	    
		      QPrime -= (QValue - 1);
		      TmpP += QPrime - (this->RIndex / 2) - NValueShift;
		    }
		}
	      
	      QPrime = Q;
	      TmpP = 0;
	      int TmpMaxP4 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP4)
		    TmpMaxP4 = TmpP;	    
		  TmpP -= QPrime - NValueShift - (this->RIndex / 2);
		  QPrime += (QValue - 1);
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP4)
			TmpMaxP4 = TmpP;	    
		      TmpP -= QPrime - NValueShift;
		      QPrime += (QValue - 1);
		    }
		}
	      if (((this->PLevel - TmpMaxP3) >= pLevel) && ((this->PLevel - TmpMaxP4) >= pLevel))
		{
		  if (Q < TmpMinQ)
		    TmpMinQ = Q;
		  if (Q > TmpMaxQ)
		    TmpMaxQ = Q;	    
		}
	    }
	}
    }
  else
    {
      int NValueShift = this->PLevel;
      int QValue = this->RIndex + 2;
      if (cftSector == 0)
	{
	  for (int Q = 0; Q < this->NbrNValue; Q += 2)
	    {
	      int QPrime = Q;
	      int TmpP = 0;
	      int TmpMaxP = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP)
		    TmpMaxP = TmpP;	    
		  QPrime -= (QValue - 2);
		  TmpP += (QPrime - this->RIndex) / 2 - NValueShift;
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue - 2);
		      TmpP += QPrime / 2 - NValueShift;
		    }
		}
	      QPrime = Q;
	      TmpP = 0;
	      int TmpMaxP2 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP2)
		    TmpMaxP2 = TmpP;	    
		  TmpP -= (QPrime / 2) - NValueShift;
		  QPrime += (QValue - 2);
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime - this->RIndex) / 2 - NValueShift;
		      QPrime += (QValue - 2);
		    }
		}
	      if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		{
		  if (Q < TmpMinQ)
		    TmpMinQ = Q;
		  if (Q > TmpMaxQ)
		    TmpMaxQ = Q;	    
		}
	    }
	}
      else
	{
	  for (int Q = 1; Q < this->NbrNValue; Q += 2)
	    {
	      int QPrime = Q;
	      int TmpP = 0;
	      int TmpMaxP = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP)
		    TmpMaxP = TmpP;	    
		  QPrime -= (QValue - 2);
		  TmpP += QPrime / 2 - NValueShift;
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue - 2);
		      TmpP += (QPrime - this->RIndex) / 2 - NValueShift;
		    }
		}
	      QPrime = Q;
	      TmpP = 0;
	      int TmpMaxP2 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP2)
		    TmpMaxP2 = TmpP;	    
		  TmpP -= (QPrime - this->RIndex) / 2 - NValueShift;
		  QPrime += (QValue - 2);
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime / 2) - NValueShift;
		      QPrime += (QValue - 2);
		    }
		}
	      if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		{
		  if (Q < TmpMinQ)
		    TmpMinQ = Q;
		  if (Q > TmpMaxQ)
		    TmpMaxQ = Q;	    
		}
	    }
	}
    }
  minQ = TmpMinQ;
  maxQ = TmpMaxQ;
  cout << "range at p=" << pLevel << ", x=" << cftSector << " : " << minQ << " " << maxQ << " (" << this->NbrNValue << ")" << endl;   
}

// get the number of particles that fit the root configuration once the number of flux quanta is fixed
// 
// nbrFluxQuanta = number of flux quanta
// padding = assume that the state has the extra padding
// return value = number of partciles

int FQHEMPSClustered2RMatrix::GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding)
{
  nbrFluxQuanta += this->RIndex + 1;
  nbrFluxQuanta *= 2;
  return (nbrFluxQuanta / (this->RIndex + 2));
}

