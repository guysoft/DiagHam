////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of MPS matrix for the Read-Rezayi k=3 state             //
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
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3Matrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/LongRationalMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"


// default constructor 
//

FQHEMPSReadRezayi3Matrix::FQHEMPSReadRezayi3Matrix()
{
}

// constructor 
//
// laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital)
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSReadRezayi3Matrix::FQHEMPSReadRezayi3Matrix(int laughlinIndex, int pLevel, int nbrBMatrices, bool cylinderFlag, double kappa)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = 2;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->CreateBMatrices();
}

// constructor from stored B matrices
//
// laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSReadRezayi3Matrix::FQHEMPSReadRezayi3Matrix(int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag, double kappa)
{
  this->RIndex = 2;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->LoadMatrices(fileName);
}


// destructor
//

FQHEMPSReadRezayi3Matrix::~FQHEMPSReadRezayi3Matrix()
{
}
  
// create the B matrices for the laughlin state
//

void FQHEMPSReadRezayi3Matrix::CreateBMatrices ()
{
  LongRational CentralCharge (4l, 5l);
  cout << "central charge = " << CentralCharge << endl;
  LongRational CentralCharge12(CentralCharge);
  CentralCharge12 /= 12l;
  LongRational WeightIdentity (0l, 1l);
  LongRational WeightPsi (2l, 3l);
  LongRational WeightW (3l);
  double WeightIdentityNumerical = WeightIdentity.GetNumericalValue();
  double WeightPsiNumerical = WeightPsi.GetNumericalValue();
  double WeightWNumerical = WeightW.GetNumericalValue();
  long* Partition = new long[2 * (this->PLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevel + 2];
  this->TemporaryOccupationNumber = new unsigned long [this->PLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductW = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductIdentity = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductPsi = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductW = new LongRationalMatrix[this->PLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi11 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi12 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi21 = new RealMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi01 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi10 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi11 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi12 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi21 = new LongRationalMatrix*[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisWLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisWRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[this->PLevel + 1];

  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, this->PLevel + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi11[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi21[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi12[i] = new RealMatrix[this->PLevel + 1];
      RationalMatrixPsi01[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi10[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi11[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi21[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi12[i] = new LongRationalMatrix[this->PLevel + 1];
    }
  
  for (int i = 0; i <= this->PLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      RationalScalarProductIdentity[i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      RationalScalarProductPsi[i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if ((3 + i) <= this->PLevel)
	RationalScalarProductW[3 + i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if (i < 3)
	RationalScalarProductW[i] = LongRationalMatrix();
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
	    LongRational Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightIdentity,
									     RationalScalarProductIdentity, i - 1, U1BosonBasis);
	    RationalScalarProductIdentity[i].SetMatrixElement(m, n, Tmp);
	    if (n != m)
	      {
		RationalScalarProductIdentity[i].SetMatrixElement(n, m, Tmp);	      
	      }
	    if ((3 + i) <= this->PLevel)
	      {
		Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightW);
		RationalScalarProductW[3 + i].SetMatrixElement(m, n, Tmp);
		if (n != m)
		  {
		    RationalScalarProductW[3 + i].SetMatrixElement(n, m, Tmp);	      
		  }
	      }
	    Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightPsi,
								RationalScalarProductPsi, i - 1, U1BosonBasis);
	    RationalScalarProductPsi[i].SetMatrixElement(m, n, Tmp);
	    if (n != m)
	      {
		RationalScalarProductPsi[i].SetMatrixElement(n, m, Tmp);	      
	      }
	  }
      ScalarProductIdentity[i] = RationalScalarProductIdentity[i];      
      ScalarProductPsi[i] = RationalScalarProductPsi[i];
      if ((3 + i) <= this->PLevel)
	ScalarProductW[3 + i] = RationalScalarProductW[3 + i];
      if (i < 3)
	ScalarProductW[i] = RealSymmetricMatrix();
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
      cout << "nbr of null vectors identity sector = " << Count << endl;
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
      if ((3 + i) <= this->PLevel)
	{	  
	  TmpMatrix.Copy(ScalarProductW[3 + i]);
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
	  cout << "nbr of null vectors W_-3 sector = " << Count << endl;
	  if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	    {
	      OrthogonalBasisWLeft[3 + i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	      OrthogonalBasisWRight[3 + i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	      Count = 0;
	      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		if (fabs(TmpDiag(n, n)) > Error)
		  {
		    OrthogonalBasisWLeft[3 + i][Count].Copy(TmpBasis[n]);
		    OrthogonalBasisWRight[3 + i][Count].Copy(TmpBasis[n]);
		    if (TmpDiag(n, n) > 0)
		      {
			OrthogonalBasisWLeft[3 + i][Count] /=  sqrt(TmpDiag(n, n));
			OrthogonalBasisWRight[3 + i][Count] /=  sqrt(TmpDiag(n, n));
		      }
		    else
		      {
			OrthogonalBasisWLeft[3 + i][Count] /=  sqrt(-TmpDiag(n, n));
			OrthogonalBasisWRight[3 + i][Count] /=  -sqrt(-TmpDiag(n, n));
		      }
		    ++Count;
	      }
	    }
	  else
	    {
	      OrthogonalBasisWLeft[3 + i] = RealMatrix();
	      OrthogonalBasisWRight[3 + i] = RealMatrix();
	    }
	}
      if (i < 3)
	{
	  OrthogonalBasisWLeft[i] = RealMatrix();
	  OrthogonalBasisWRight[i] = RealMatrix();
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
      cout << "nbr of null vectors Psi sector = " << Count << endl;
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

  cout << "computing Psi matrix elements" << endl;
  LongRational Weight (WeightPsi);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int j = 0; j <= this->PLevel; ++j)
	{
	  RationalMatrixPsi01[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  RationalMatrixPsi10[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  RationalMatrixPsi11[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  if ((3 + j) <= this->PLevel)
	    {	  
	      RationalMatrixPsi12[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	    }
	  if ((3 + i) <= this->PLevel)
	    {	  
	      RationalMatrixPsi21[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	    }
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
		LongRational Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightIdentity, WeightPsi, Weight,
									 RationalMatrixPsi01, i, j - 1, U1BosonBasis);
		RationalMatrixPsi01[i][j].SetMatrixElement(n, m, Tmp);
		Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightPsi, WeightIdentity, Weight,
							    RationalMatrixPsi10, i, j - 1, U1BosonBasis);
		RationalMatrixPsi10[i][j].SetMatrixElement(n, m, Tmp);
		Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightPsi, WeightPsi, Weight,
							    RationalMatrixPsi11, i, j - 1, U1BosonBasis);
		RationalMatrixPsi11[i][j].SetMatrixElement(n, m, Tmp);
		if ((3 + j) <= this->PLevel)
		  {	  
		    Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightPsi, WeightW, Weight,
								RationalMatrixPsi12, i - 1, j - 1, U1BosonBasis);
		    RationalMatrixPsi12[i][j].SetMatrixElement(n, m, Tmp);
		  }
		if ((3 + i) <= this->PLevel)
		  {	  
		    Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightW, WeightPsi, Weight,
								RationalMatrixPsi21, i - 1, j - 1, U1BosonBasis);
		    RationalMatrixPsi21[i][j].SetMatrixElement(n, m, Tmp);
		  }
	      }
	  MatrixPsi01[i][j] = RationalMatrixPsi01[i][j];
	  MatrixPsi10[i][j] = RationalMatrixPsi10[i][j];
	  MatrixPsi11[i][j] = RationalMatrixPsi11[i][j];	  
	  if ((3 + j) <= this->PLevel)
	    {	
	      MatrixPsi12[i][3 + j] = RationalMatrixPsi12[i][j];
	      MatrixPsi12[i][3 + j] *= -sqrt(26.0) / 9.0;
	    }
	  if ((3 + i) <= this->PLevel)
	    {	  
	      MatrixPsi21[3 + i][j] = RationalMatrixPsi21[i][j];
	      MatrixPsi21[3 + i][j] *= sqrt(26.0) / 9.0;
	    }
	}
    }

  SparseRealMatrix* BMatrices = new SparseRealMatrix[this->NbrBMatrices];

  int** StartingIndexPerPLevel = new int* [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel[0] = 0;
  StartingIndexPerPLevel[0] = new int [1];
  StartingIndexPerPLevel[0][0] = 0;

  
  int NValueShift;
  int QValue;
  int QValueDenominator;
  QValue = 5;
  QValueDenominator = 3;
  this->NbrNValue = QValueDenominator * (2 * this->PLevel) + 4 + 2 + 1;
  NValueShift = QValueDenominator * this->PLevel;

     
  this->NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() 
									   + OrthogonalBasisPsiLeft[0].GetNbrColumn() 
									   + OrthogonalBasisPsiLeft[0].GetNbrColumn() 
									   + OrthogonalBasisWLeft[0].GetNbrColumn())) * this->NbrNValue;
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      StartingIndexPerPLevel[i] = new int [i + 1];      
      StartingIndexPerPLevel[i][0] = this->TotalStartingIndexPerPLevel[i];
      int Tmp = 0;
      int Tmp2;
      for (int j = 0; j < i; ++j)
	{
	  Tmp2 = U1BosonBasis[j]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[i - j].GetNbrColumn() 
								+ OrthogonalBasisPsiLeft[i - j].GetNbrColumn() 
								+ OrthogonalBasisPsiLeft[i - j].GetNbrColumn() 
								+ OrthogonalBasisWLeft[i - j].GetNbrColumn()) * this->NbrNValue;
	  StartingIndexPerPLevel[i][j + 1] = Tmp2 + StartingIndexPerPLevel[i][j];
	  Tmp += Tmp2;
	}
      Tmp += U1BosonBasis[i]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() 
							    + OrthogonalBasisPsiLeft[0].GetNbrColumn()
							    + OrthogonalBasisPsiLeft[0].GetNbrColumn() 
							    + OrthogonalBasisWLeft[0].GetNbrColumn()) * this->NbrNValue;
      this->NbrIndicesPerPLevel[i] =  Tmp;
    }
  
  int MatrixSize = this->NbrIndicesPerPLevel[this->PLevel] + this->TotalStartingIndexPerPLevel[this->PLevel];
  cout << "B matrix size = " << MatrixSize << endl;
  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize);
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
	      for (int j = QValueDenominator; j < this->NbrNValue; ++j)
		{
		  // identity 
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
								     + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								     + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue * QValueDenominator))));
			  BMatrices[0].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							this->GetReadRezayiK3MatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			}
		    }
		  // Psi_{+1} and Psi_{-1}
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
								     + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								     + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			  BMatrices[0].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							this->GetReadRezayiK3MatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			  BMatrices[0].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 3, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							this->GetReadRezayiK3MatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 3, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			}
		    }
		}
	    }
	}
    }

  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= (i - 3); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p - 3];
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisWLeft = OrthogonalBasisWLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisWRight = OrthogonalBasisWRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductW = ScalarProductW[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      for (int j = QValueDenominator; j < this->NbrNValue; ++j)
		{
		  // W
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisWLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisWLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductW(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisWRight(NeutralIndex4, NeutralIndex2);
				}
			      Tmp += TmpOrthogonalBasisWLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  if (this->CylinderFlag)
			    Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
								     + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								     + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			  BMatrices[0].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 5, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							this->GetReadRezayiK3MatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 5, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
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
		  RealMatrix& TmpMatrixPsi11 = MatrixPsi11[i - p][j - q];
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
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
				  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  if (this->CylinderFlag)
				    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
										   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
				  BMatrices[1].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), TmpOrthogonalBasisPsi1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 3, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), TmpOrthogonalBasisPsi2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}

			    }
			  N2 = QValueDenominator * (j - i) + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
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
				  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  if (this->CylinderFlag)
				    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
										   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
				  BMatrices[1].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), TmpOrthogonalBasisPsi1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), TmpOrthogonalBasisPsi2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}
			    }

			  N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi11(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  if (this->CylinderFlag)
				    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightPsiNumerical + WeightPsiNumerical + ((double) (i + j))
										   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
				  BMatrices[1].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 3, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), TmpOrthogonalBasisPsi1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), TmpOrthogonalBasisPsi2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}
			    }
			}
		    }
		}	      
	    }
	}
    }
  
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= (i - 3); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 3];
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		  RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		  RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisW1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi21(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisW1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  if (this->CylinderFlag)
				    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
										   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
				  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 5, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), TmpOrthogonalBasisPsi1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 3, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), TmpOrthogonalBasisPsi2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}

			    }
			}
		    }
		}	      
	    }
	}
    }
  
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= (j - 3); ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 3];
		  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		  RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		  RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisW2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi12(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisW2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  if (this->CylinderFlag)
				    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
										   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
				  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), TmpOrthogonalBasisPsi1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 5, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), TmpOrthogonalBasisPsi2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}
			    }
			}
		    }
		}	      
	    }
	}
    }

  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->RealBMatrices[i] = BMatrices[i];
    }
  delete[] BMatrices;

  delete[] ScalarProductIdentity;
  delete[] ScalarProductPsi;
  delete[] ScalarProductW;
  delete[] RationalScalarProductIdentity;
  delete[] RationalScalarProductPsi;
  delete[] RationalScalarProductW;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      delete[] MatrixPsi01[i];
      delete[] MatrixPsi10[i];
      delete[] MatrixPsi21[i];
      delete[] MatrixPsi12[i];
      delete[] MatrixPsi11[i];
      delete[] RationalMatrixPsi01[i];
      delete[] RationalMatrixPsi10[i];
      delete[] RationalMatrixPsi21[i];
      delete[] RationalMatrixPsi12[i];
      delete[] RationalMatrixPsi11[i];
      delete U1BosonBasis[i];
    }
  delete[] U1BosonBasis;
  delete[] MatrixPsi01;
  delete[] MatrixPsi10;
  delete[] MatrixPsi21;
  delete[] MatrixPsi12;
  delete[] MatrixPsi11;
  delete[] RationalMatrixPsi01;
  delete[] RationalMatrixPsi10;
  delete[] RationalMatrixPsi21;
  delete[] RationalMatrixPsi12;
  delete[] RationalMatrixPsi11;
  delete[] OrthogonalBasisIdentityLeft;
  delete[] OrthogonalBasisPsiLeft;
  delete[] OrthogonalBasisWLeft;
  delete[] OrthogonalBasisIdentityRight;
  delete[] OrthogonalBasisPsiRight;
  delete[] OrthogonalBasisWRight;
}

// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSReadRezayi3Matrix::GetBondIndexRange(int pLevel, int qValue)
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

int FQHEMPSReadRezayi3Matrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  return (this->TotalStartingIndexPerPLevel[pLevel] + (localIndex * this->NbrNValue + qValue));
}

