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
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/LongRationalMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"

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

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, bool cylinderFlag, double kappa)
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
  this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
  this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));
  this->MatrixElementNormalization = 1.0 / M_SQRT2;
  this->CreateBMatrices();
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
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
  this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));
  this->MatrixElementNormalization = 1.0 / M_SQRT2;
}

// destructor
//

FQHEMPSClustered2RQuasiholeSectorMatrix::~FQHEMPSClustered2RQuasiholeSectorMatrix()
{
}
  
// create the B matrices for the laughlin state
//

void FQHEMPSClustered2RQuasiholeSectorMatrix::CreateBMatrices ()
{
  bool SelfDualFlag = ((this->RIndex & 1) == 0);
  LongRational CentralCharge12 ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  cout << "central charge = " << CentralCharge12 << endl;
  CentralCharge12 /= 12l;
  double WeightSigmaNumerical = this->WeightSigma.GetNumericalValue();
  double WeightPhiNumerical = this->WeightPhi.GetNumericalValue();
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

  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, this->PLevel + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevel + 1];
      RationalMatrixPsi01[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi10[i] = new LongRationalMatrix[this->PLevel + 1];
    }

  cout << "weight: " <<   this->WeightSigma << " " << this->WeightPhi << endl;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      RationalScalarProductSigma[i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if (SelfDualFlag == false)
	RationalScalarProductPhi[i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
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
 	    LongRational Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, this->WeightSigma,
 									     RationalScalarProductSigma, i - 1, U1BosonBasis);
	    RationalScalarProductSigma[i].SetMatrixElement(m, n, Tmp);
	    if (n != m)
	      {
		RationalScalarProductSigma[i].SetMatrixElement(n, m, Tmp);	      
	      }
	    if (SelfDualFlag == false)
	      {
		Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, this->WeightPhi,
								    RationalScalarProductPhi, i - 1, U1BosonBasis);
		RationalScalarProductPhi[i].SetMatrixElement(m, n, Tmp);
		if (n != m)
		  {
		    RationalScalarProductPhi[i].SetMatrixElement(n, m, Tmp);	      
		  }
	      }
	  }      
      ScalarProductSigma[i] = RationalScalarProductSigma[i];      
      if (SelfDualFlag == false)
	ScalarProductPhi[i] = RationalScalarProductPhi[i];
      
      RealSymmetricMatrix TmpMatrix;
      TmpMatrix.Copy(ScalarProductSigma[i]);
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

      if (SelfDualFlag == false)
	{
	  TmpMatrix.Copy(ScalarProductPhi[i]);
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

  this->StartingIndexPerPLevel = new int* [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->IdentityBasisDimension = new int [this->PLevel + 1];
  this->PsiBasisDimension = new int [this->PLevel + 1];
  this->U1BasisDimension = new int [this->PLevel + 1];	
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->IdentityBasisDimension[i] = OrthogonalBasisSigmaLeft[i].GetNbrColumn();
      if (SelfDualFlag == false)
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

     
  if (SelfDualFlag == false)
    {
      this->NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * OrthogonalBasisSigmaLeft[0].GetNbrColumn()) * this->NbrNValue;
    }
  else
    {
      this->NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * (OrthogonalBasisSigmaLeft[0].GetNbrColumn() + OrthogonalBasisPhiLeft[0].GetNbrColumn())) * this->NbrNValue;
    }
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      this->StartingIndexPerPLevel[i] = new int [i + 1];      
      this->StartingIndexPerPLevel[i][0] = this->TotalStartingIndexPerPLevel[i];
      int Tmp = 0;
      int Tmp2;
      if (SelfDualFlag == false)
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
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int j = 0; j <= this->PLevel; ++j)
	{
	  RationalMatrixPsi01[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  if (SelfDualFlag == false)
	    RationalMatrixPsi10[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  cout << "Levels = " <<  i << " " << j << endl;
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
		LongRational Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, this->WeightSigma, this->WeightPhi, this->WeightPrimaryFieldMatrixElement,
									 RationalMatrixPsi01, i - 1, j - 1, U1BosonBasis);
		RationalMatrixPsi01[i][j].SetMatrixElement(n, m, Tmp);
		if (SelfDualFlag == false)
		  {
		    Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, this->WeightPhi, this->WeightSigma, this->WeightPrimaryFieldMatrixElement,
								RationalMatrixPsi10, i - 1, j - 1, U1BosonBasis);
		    RationalMatrixPsi10[i][j].SetMatrixElement(n, m, Tmp);
		  }
	      }
	  MatrixPsi01[i][j] = RationalMatrixPsi01[i][j];
	  MatrixPsi01[i][j] *= MatrixElementNormalization;
	  if (SelfDualFlag == false)
	    {
	      MatrixPsi10[i][j] = RationalMatrixPsi10[i][j];
	      MatrixPsi10[i][j] *= MatrixElementNormalization;
	    }
	}
    }
  cout << "building B matrices" << endl;


  SparseRealMatrix* BMatrices = new SparseRealMatrix[this->NbrBMatrices];
  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize);
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
		      if (SelfDualFlag == false)
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
		  for (int j = 2; j < this->NbrNValue; ++j)
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
									 + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 2, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
		      if (SelfDualFlag == false)
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
  
  FactorialCoefficient Coef;
  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  BMatrices[1] = SparseRealMatrix(MatrixSize, MatrixSize);
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
			  if (SelfDualFlag == false)
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
			  if (SelfDualFlag == false)
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
			  else
			    {
// 			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
// 				{
// 				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
// 				    {
// 				      double Tmp = 0.0;
// 				      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
// 					{
// 					  double Tmp1 = 0.0;			      
// 					  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
// 					    {
// 					      Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
// 					    }
// 					  Tmp += TmpOrthogonalBasisPhi1(NeutralIndex3, NeutralIndex1) * Tmp1;
// 					}
// 				      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
// 				      if (this->CylinderFlag)
// 					Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightSigmaNumerical + WeightPhiNumerical + ((double) (i + j))
// 										       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))
// 										       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift)) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))));
// 				      BMatrices[1].SetMatrixElement(this->Get2RMatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisSigma1.GetNbrColumn(), this->StartingIndexPerPLevel[i][p]), 
// 								    this->Get2RMatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisSigma2.GetNbrColumn(), this->StartingIndexPerPLevel[j][q]), Tmp);
// 				    }
// 				}
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
    }
  delete[] U1BosonBasis;
  delete[] MatrixPsi01;
  delete[] MatrixPsi10;
  delete[] OrthogonalBasisSigmaLeft;
  delete[] OrthogonalBasisPhiLeft;
  delete[] OrthogonalBasisSigmaRight;
  delete[] OrthogonalBasisPhiRight;
}


