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

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, bool cylinderFlag, double kappa, 
						   AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
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
  
FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool cylinderFlag, double kappa, 
						   AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
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

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int pLevel, int nbrBMatrices, char* fileName, bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
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
      ErrorFlag = StateDefinition.GetAsSingleLongRational("CentralCharge", this->CentralCharge);
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
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag, double kappa)
{
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
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

  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, this->PLevel + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevel + 1];
      RationalMatrixPsi01[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi10[i] = new LongRationalMatrix[this->PLevel + 1];
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
      RationalScalarProductIdentity[i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      RationalScalarProductPsi[i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if (cftDirectory != 0)
	{
	  sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_scalarproducts_identity_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	  sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_scalarproducts_psi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	}
      if ((cftDirectory != 0) && (IsFile(TmpScalarProductIdentityFileName)) && (IsFile(TmpScalarProductPsiFileName)))
	{
	  RationalScalarProductIdentity[i].ReadMatrix(TmpScalarProductIdentityFileName);
	  RationalScalarProductPsi[i].ReadMatrix(TmpScalarProductPsiFileName);
	}
      else
	{
	  if (architecture == 0)
	    {
	      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		for (int m = n; m < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++m)
		  {
		    int PartitionLength = 0;
		    U1BosonBasis[i]->GetOccupationNumber(n, TmpPartition);	    
		    for (int k = 1; k <= i; ++k)
		      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			{
			  ++PartitionLength;
			}
		    int Position = PartitionLength;
		    PartitionLength = 0;
		    for (int k = 1; k <= i; ++k)
		      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			{
			  Partition[Position - PartitionLength - 1] = (long) k;
			  ++PartitionLength;
			}
		    U1BosonBasis[i]->GetOccupationNumber(m, TmpPartition);	    
		    for (int k = 1; k <= i; ++k)
		      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			{
			  Partition[PartitionLength] = -(long) k;
			  ++PartitionLength;		  
			}
		    LongRational Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, this->WeightIdentity,
										     RationalScalarProductIdentity, i - 1, U1BosonBasis, this->TemporaryOccupationNumber);
		    RationalScalarProductIdentity[i].SetMatrixElement(m, n, Tmp);
		    if (n != m)
		      {
			RationalScalarProductIdentity[i].SetMatrixElement(n, m, Tmp);	      
		      }
		    Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, this->WeightPsi,
									RationalScalarProductPsi, i - 1, U1BosonBasis, this->TemporaryOccupationNumber);
		    RationalScalarProductPsi[i].SetMatrixElement(m, n, Tmp);
		    if (n != m)
		      {
			RationalScalarProductPsi[i].SetMatrixElement(n, m, Tmp);	      
		      }
		  }
	      if (cftDirectory != 0)
		{
		  RationalScalarProductIdentity[i].WriteMatrix(TmpScalarProductIdentityFileName);
		  RationalScalarProductPsi[i].WriteMatrix(TmpScalarProductPsiFileName);
		}
	    }    
	  else
	    {
	      if ((cftDirectory != 0) && (IsFile(TmpScalarProductIdentityFileName)))
		{		
		  RationalScalarProductIdentity[i].ReadMatrix(TmpScalarProductIdentityFileName);
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12, 
							 this->WeightIdentity,
							 RationalScalarProductIdentity,  i- 1);
		  Operation1.ApplyOperation(architecture);
		  RationalScalarProductIdentity[i] = Operation1.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalScalarProductIdentity[i].WriteMatrix(TmpScalarProductIdentityFileName);
		    }
		}
	      if ((cftDirectory != 0) && (IsFile(TmpScalarProductPsiFileName)))
		{
		  RationalScalarProductPsi[i].ReadMatrix(TmpScalarProductPsiFileName);
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12, 
							 this->WeightPsi,
							 RationalScalarProductPsi,  i - 1);
		  Operation2.ApplyOperation(architecture);
		  RationalScalarProductPsi[i] = Operation2.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalScalarProductPsi[i].WriteMatrix(TmpScalarProductPsiFileName);
		    }
		}
	    }
	}
      ScalarProductIdentity[i] = RationalScalarProductIdentity[i];      
      ScalarProductPsi[i] = RationalScalarProductPsi[i];
      
      RealSymmetricMatrix TmpMatrix;
      TmpMatrix.Copy(ScalarProductIdentity[i]);
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
	if (fabs(TmpDiag(n, n)) < Error)
	  ++Count;
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

      TmpMatrix.Copy(ScalarProductPsi[i]);
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

  this->StartingIndexPerPLevel = new int* [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->IdentityBasisDimension = new int [this->PLevel + 1];
  this->PsiBasisDimension = new int [this->PLevel + 1];
  this->U1BasisDimension = new int [this->PLevel + 1];	
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->IdentityBasisDimension[i] = OrthogonalBasisIdentityLeft[i].GetNbrColumn();
      this->PsiBasisDimension[i] = OrthogonalBasisPsiLeft[i].GetNbrColumn();
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

  this->NbrNValuesPerPLevel = new int [this->PLevel + 1];
  this->NInitialValuePerPLevel = new int [this->PLevel + 1];
  this->NLastValuePerPLevel = new int [this->PLevel + 1];     
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->ComputeChargeIndexRange(i, this->NInitialValuePerPLevel[i], this->NLastValuePerPLevel[i]);
      this->NbrNValuesPerPLevel[i] =  this->NLastValuePerPLevel[i] - this->NInitialValuePerPLevel[i] + 1;
    }

  this->NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() + OrthogonalBasisPsiLeft[0].GetNbrColumn())) * this->NbrNValuesPerPLevel[0];
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      this->StartingIndexPerPLevel[i] = new int [i + 1];      
      this->StartingIndexPerPLevel[i][0] = this->TotalStartingIndexPerPLevel[i];
      int Tmp = 0;
      int Tmp2;
      for (int j = 0; j < i; ++j)
	{
	  Tmp2 = U1BosonBasis[j]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[i - j].GetNbrColumn() + OrthogonalBasisPsiLeft[i - j].GetNbrColumn()) * this->NbrNValuesPerPLevel[i];
	  this->StartingIndexPerPLevel[i][j + 1] = Tmp2 + this->StartingIndexPerPLevel[i][j];
	  Tmp += Tmp2;
	}
      Tmp += U1BosonBasis[i]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() + OrthogonalBasisPsiLeft[0].GetNbrColumn()) * this->NbrNValuesPerPLevel[i];
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
	  RationalMatrixPsi01[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  RationalMatrixPsi10[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  cout << "Levels = " <<  i << " " << j << endl;
	  if (cftDirectory != 0)
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_matrixelement_identitypsi_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
	      sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_matrixelement_psiidentity_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpScalarProductIdentityFileName)) && (IsFile(TmpScalarProductPsiFileName)))
	    {
	      RationalMatrixPsi01[i][j].ReadMatrix(TmpScalarProductIdentityFileName);
	      RationalMatrixPsi10[i][j].ReadMatrix(TmpScalarProductPsiFileName);
	    }
	  else
	    {
	      if (architecture == 0)
		{
		  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		    for (int m = 0; m < U1BosonBasis[j]->GetHilbertSpaceDimension(); ++m)
		      {
			int PartitionLength = 0;
			U1BosonBasis[i]->GetOccupationNumber(n, TmpPartition);	    
			for (int k = 1; k <= i; ++k)
			  for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			    {
			      ++PartitionLength;
			    }
			int Position = PartitionLength;
			PartitionLength = 0;
			for (int k = 1; k <= i; ++k)
			  for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			    {
			      Partition[Position - PartitionLength - 1] = (long) k;
			      ++PartitionLength;
			    }
			U1BosonBasis[j]->GetOccupationNumber(m, TmpPartition);	    
			for (int k = 1; k <= j; ++k)
			  for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			    {
			      Partition[PartitionLength] = -(long) k;
			      ++PartitionLength;		  
			    }
			LongRational Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, 
										 this->WeightIdentity, this->WeightPsi, this->WeightPrimaryFieldMatrixElement,
										 RationalMatrixPsi01, i - 1, j, U1BosonBasis, this->TemporaryOccupationNumber);
			RationalMatrixPsi01[i][j].SetMatrixElement(n, m, Tmp);
			Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, 
								    this->WeightPsi, this->WeightIdentity, this->WeightPrimaryFieldMatrixElement,
								    RationalMatrixPsi10, i - 1, j, U1BosonBasis, this->TemporaryOccupationNumber);
			RationalMatrixPsi10[i][j].SetMatrixElement(n, m, Tmp);
		      }
		  if (cftDirectory != 0)
		    {
		      RationalMatrixPsi01[i][j].WriteMatrix(TmpScalarProductIdentityFileName);
		      RationalMatrixPsi10[i][j].WriteMatrix(TmpScalarProductPsiFileName);
		    }
		}
	      else
		{
		  if ((cftDirectory != 0) && (IsFile(TmpScalarProductIdentityFileName)))
		    {
		      RationalMatrixPsi01[i][j].ReadMatrix(TmpScalarProductIdentityFileName);
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, j, CentralCharge12, 
							     this->WeightIdentity, this->WeightPsi, this->WeightPrimaryFieldMatrixElement,
							     RationalMatrixPsi01,  i - 1, j);
		      Operation1.ApplyOperation(architecture);
		      RationalMatrixPsi01[i][j] = Operation1.GetMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalMatrixPsi01[i][j].WriteMatrix(TmpScalarProductIdentityFileName);
			}
		    }
		  if ((cftDirectory != 0) && (IsFile(TmpScalarProductPsiFileName)))
		    {
		      RationalMatrixPsi10[i][j].ReadMatrix(TmpScalarProductPsiFileName);
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							     this->WeightPsi, this->WeightIdentity, this->WeightPrimaryFieldMatrixElement,
							     RationalMatrixPsi10,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      RationalMatrixPsi10[i][j] = Operation2.GetMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalMatrixPsi10[i][j].WriteMatrix(TmpScalarProductPsiFileName);
			}
		    }
		}	      
	    }
	  MatrixPsi01[i][j] = RationalMatrixPsi01[i][j];
	  MatrixPsi01[i][j] *= MatrixElementNormalization;
	  MatrixPsi10[i][j] = RationalMatrixPsi10[i][j];
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
		  for (int j = 1; j < this->NbrNValue; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndex(j - 1, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
			    }
			}
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndex(j - 1, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
			    }
			}
		    }
		}
	      else
		{
		  for (int j = 2; j < this->NbrNValue; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndex(j - 2, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
			    }
			}
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndex(j - 2, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
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
		  for (int j = 1; j < this->NbrNValue; ++j)
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
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 1, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
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
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 1, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
		    }
		}
	      else
		{
		  for (int j = 2; j < this->NbrNValue; ++j)
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
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 2, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
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
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 2, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), Tmp);
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
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
				{
				  ++TmpNbrElementPerRow[this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
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
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
				{
				  ++TmpNbrElementPerRow[this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p])];
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
				  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
								this->Get2RMatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), this->StartingIndexPerPLevel[j][q]), Tmp);
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
				  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
								this->Get2RMatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), this->StartingIndexPerPLevel[j][q]), Tmp);
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
    }
  delete[] TmpNbrElementPerRow;
  delete[] U1BosonBasis;
  delete[] MatrixPsi01;
  delete[] MatrixPsi10;
  delete[] OrthogonalBasisIdentityLeft;
  delete[] OrthogonalBasisPsiLeft;
  delete[] OrthogonalBasisIdentityRight;
  delete[] OrthogonalBasisPsiRight;
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

// compute the matrix elements of any primary field in the Virasoro descendant basis
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// descendantPosition = location of the primary field
// position = position in partition up to which  all the indices are positive
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
  
// LongRational FQHEMPSClustered2RMatrix::ComputeDescendantMatrixElementReverse (long* partition, int partitionLength, 
// 									      int descendantPosition, int position,
// 									      LongRational& centralCharge12, LongRational& weight1, 
// 									      LongRational& weight2, LongRational& weight,
// 									      LongRationalMatrix** precomputedDescendantMatrixElement, 
// 									      int precomputedDescendantMatrixElementMaxLeftPLevel, 
// 									      int precomputedDescendantMatrixElementMaxRightPLevel, 
// 									      BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber)
// {
//   if (partitionLength == 0)
//     {
//       return 1l;
//     }
//   if ((descendantPosition < partitionLength) && (partition[partitionLength - 1] > 0))
//     {
//       return 0l;
//     }
//   while ((position < partitionLength) && (partition[position] > 0l))
//     ++position;
//   if (descendantPosition == partitionLength) 
//     {
//       LongRational Tmp(1l);
//       LongRational TmpSum = weight1;
//       TmpSum -= weight2;
//       LongRational Tmp2;
//       for (int i = partitionLength - 1; i >= 0; --i)
// 	{
// 	  Tmp2 = weight;
// 	  Tmp2 *= partition[i];
// 	  long Tmp3 = 0l;
// 	  for (int j = 0; j < i; ++j)
// 	    Tmp3 += partition[j];
// 	  Tmp2 += TmpSum;
// 	  Tmp2 += Tmp3;	  
// 	  Tmp *= Tmp2;
// 	}
//       return Tmp;
//     }
//   if (position == partitionLength)
//     {
//       LongRational Tmp(1l);
//       LongRational TmpSum = weight1;
//       TmpSum -= weight2;
//       LongRational Tmp2;
//       for (int i = 0; i < partitionLength ; ++i)
// 	{
// 	  Tmp2 = weight;
// 	  Tmp2 *= partition[i];
// 	  long Tmp3 = 0l;
// 	  for (int j = i + 1; j < partitionLength; ++j)
// 	    Tmp3 -= partition[j];
// 	  Tmp2 += TmpSum;
// 	  Tmp2 -= Tmp3;
// 	  Tmp *= Tmp2;
// 	}
//       Tmp *= 1l - (2l * (partitionLength & 1l));
//       return Tmp;
//     }
//   if (descendantPosition == position)
//     {
//       int TmpPosition = position + 1;
//       while ((TmpPosition < partitionLength) && (partition[TmpPosition] <= 0))
// 	++TmpPosition;
//       if (TmpPosition == partitionLength)
// 	{
// 	  TmpPosition = 1;
// 	  int TmpPLevel1 = partition[0];
// 	  bool FlagSorted = true;
// 	  while (TmpPosition <= position)
// 	    {
// 	      TmpPLevel1 += partition[TmpPosition];
// 	      if (partition[TmpPosition - 1] < partition[TmpPosition])
// 		FlagSorted = false;
// 	      ++TmpPosition;
// 	    }
// 	  if ((TmpPLevel1 <= precomputedDescendantMatrixElementMaxLeftPLevel) && (FlagSorted == true))
// 	    {
// 	      int TmpPLevel2 = -partition[TmpPosition];	  
// 	      FlagSorted = true;
// 	      ++TmpPosition;
// 	      while (TmpPosition < partitionLength)
// 		{
// 		  TmpPLevel2 -= partition[TmpPosition];
// 		  if (partition[TmpPosition - 1] < partition[TmpPosition])
// 		    FlagSorted = false;
// 		  ++TmpPosition;
// 		}
// 	      if ((TmpPLevel2 <= precomputedDescendantMatrixElementMaxRightPLevel) && (FlagSorted == true))
// 		{
// 		  for (int k = 0; k <= (this->PLevel + 1); ++k)
// 		    temporaryOccupationNumber[k] = 0x0ul;	  
// 		  for (TmpPosition = 0; TmpPosition <= position; ++TmpPosition)
// 		    {
// 		      temporaryOccupationNumber[partition[TmpPosition]]++;	      
// 		    }
// 		  temporaryOccupationNumber[0] = TmpPLevel1 - position;
// 		  int TmpIndex1 = basis[TmpPLevel1]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
// 		  for (int k = 0; k <= (this->PLevel + 1); ++k)
// 		    temporaryOccupationNumber[k] = 0x0ul;	  
// 		  for (TmpPosition = position + 1; TmpPosition < partitionLength; ++TmpPosition)
// 		    {
// 		      temporaryOccupationNumber[-partition[TmpPosition]]++;	      
// 		    }
// 		  temporaryOccupationNumber[0] = TmpPLevel2 - partitionLength + position;
// 		  int TmpIndex2 = basis[TmpPLevel2]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
// 		  LongRational Tmp;
// 		  precomputedDescendantMatrixElement[TmpPLevel1][TmpPLevel2].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
// 		  return Tmp;
// 		}
// 	    }
// 	}
//     }
//   if (descendantPosition > position)
//     {
//       LongRational Tmp(0l);
//       if ((partition[position - 1] + partition[position]) == 0)
// 	{
// 	  long TmpLength = 0l;
// 	  long Store = partition[position - 1];
// 	  for (int i = position + 1; i < partitionLength; ++i)
// 	    TmpLength += partition[i];
// 	  for (int i = position + 1; i < partitionLength; ++i)
// 	    partition[i - 2] = partition[i];
// 	  Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
// 		   + (2l * Store) * (weight2 - TmpLength)) * 
// 		  this->ComputeDescendantMatrixElement(partition, partitionLength - 2, descendantPosition, 
// 						       position - 1, centralCharge12, weight1, weight2, weight,
// 						       precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
// 						       precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
// 	  for (int i = partitionLength - 1; i > position; --i)
// 	    partition[i] = partition[i - 2];
// 	  partition[position - 1] = Store;
// 	  partition[position] = -Store;
// 	}
//       else
// 	{
//  	  long Store1 = partition[position - 1];
// 	  long Store2 = partition[position];
// 	  partition[position - 1] += partition[position];
// 	  for (int i = position + 1; i < partitionLength; ++i)
// 	    partition[i - 1] = partition[i];
// 	  if ((Store1 + Store2) > 0)
// 	    {
// 	      Tmp += ((Store1 - Store2) 
// 		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
// 							     position, centralCharge12, weight1, weight2, weight,
// 							     precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
// 							     precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
// 	    }
// 	  else
// 	    {
// 	      Tmp += ((Store1 - Store2) 
// 		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
// 							     position - 1, centralCharge12, weight1, weight2, weight,
// 							     precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
// 							     precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
// 	    }
// 	  for (int i = partitionLength - 1; i > position; --i)
// 	    partition[i] = partition[i - 1];
// 	  partition[position] = Store2;
// 	  partition[position - 1] = Store1;
// 	}
      
//       long Store1 = partition[position - 1];
//       partition[position - 1] = partition[position];
//       partition[position] = Store1;
//       Tmp += this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition, position + 1, 
// 						  centralCharge12, weight1, weight2, weight,
// 						  precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
// 						  precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);
//       Store1 = partition[position - 1];
//       partition[position - 1] = partition[position];
//       partition[position] = Store1;
//       return Tmp;  
//     }

//   LongRational Tmp1 = this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition - 1, position, 
// 							   centralCharge12, weight1, weight2, weight,
// 							   precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
// 							   precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);
//   LongRational Tmp2 = weight;
//   Tmp2 *= partition[position - 1];
//   Tmp2 += weight1;
//   Tmp2 -= weight2;
//   long TmpLength = 0l;
//   for (int i = position; i < partitionLength; ++i)
//     TmpLength += partition[i];
//   for (int i = 0; i < (position - 1); ++i)
//     TmpLength += partition[i];
//   Tmp2 += TmpLength;
//   long Store = partition[position - 1];
//   for (int i = position; i < partitionLength; ++i)
//     partition[i - 1] = partition[i];
//   Tmp2 *= this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition - 1, position - 1, 
// 					       centralCharge12, weight1, weight2, weight,
// 					       precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
// 					       precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);
//   for (int i = partitionLength - 1; i >= position; --i)
//     partition[i] = partition[i - 1];
//   partition[position - 1] = Store;
//   Tmp1 += Tmp2;
//   return Tmp1;
// }

// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSClustered2RMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  if ((pLevel < 0) || (pLevel > this->PLevel) || (qValue < 0) || (qValue >= this->NbrNValue))
    return 0;
  return this->NbrIndicesPerPLevel[pLevel] / this->NbrNValue;  
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int FQHEMPSClustered2RMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  return (this->TotalStartingIndexPerPLevel[pLevel] + (localIndex * this->NbrNValue + qValue));
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
  ReadLittleEndian(file, this->NbrNValue);
  ReadLittleEndian(file, this->CylinderFlag);
  ReadLittleEndian(file, this->Kappa);
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->TotalStartingIndexPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->NbrIndicesPerPLevel[i]);
    }
  return true;
}

// save the specific informations to the file header 
// 
// file = reference on the output file stream
// return value = true if no error occurred  

bool FQHEMPSClustered2RMatrix::SaveHeader (ofstream& file)
{
  int HeaderSize = (this->PLevel + 1) * (2 * sizeof(int)) + (sizeof(int) * 4) + sizeof(bool) + sizeof(double);
  WriteLittleEndian(file, HeaderSize);
  WriteLittleEndian(file, this->PLevel);
  WriteLittleEndian(file, this->LaughlinIndex);
  WriteLittleEndian(file, this->RIndex);
  WriteLittleEndian(file, this->NbrNValue);
  WriteLittleEndian(file, this->CylinderFlag);
  WriteLittleEndian(file, this->Kappa);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->TotalStartingIndexPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->NbrIndicesPerPLevel[i]);
    }
  return true;
}

