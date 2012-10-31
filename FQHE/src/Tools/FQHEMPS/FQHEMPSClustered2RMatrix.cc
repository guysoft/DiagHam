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
#include "Matrix/SparseRealMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"


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
// nbrBMatrices = number of B matrices to compute (max occupation per orbital)
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, bool cylinderFlag, double kappa)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->CreateBMatrices();
}

// destructor
//

FQHEMPSClustered2RMatrix::~FQHEMPSClustered2RMatrix()
{
}
  
// create the B matrices for the laughlin state
//

void FQHEMPSClustered2RMatrix::CreateBMatrices ()
{
  LongRational CentralCharge12 ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  cout << "central charge = " << CentralCharge12 << endl;
  CentralCharge12 /= 12l;
  LongRational WeightIdentity (0l, 1l);
  LongRational WeightPsi (this->RIndex, 4l);
  long* Partition = new long[2 * (this->PLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[this->PLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[this->PLevel + 1];

  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, this->PLevel + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevel + 1];
    }
  
  for (int i = 0; i <= this->PLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      ScalarProductIdentity[i] = RealSymmetricMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      ScalarProductPsi[i] = RealSymmetricMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
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
	    LongRational Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightIdentity);
	    ScalarProductIdentity[i].SetMatrixElement(m, n, Tmp.GetNumericalValue());
	    Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightPsi);
	    ScalarProductPsi[i].SetMatrixElement(m, n, Tmp.GetNumericalValue());
	  }
      
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
	  MatrixPsi01[i][j] = RealMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  MatrixPsi10[i][j] = RealMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
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
		LongRational Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightIdentity, WeightPsi, Weight);
		MatrixPsi01[i][j].SetMatrixElement(n, m, Tmp.GetNumericalValue());
		Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightPsi, WeightIdentity, Weight);
		MatrixPsi10[i][j].SetMatrixElement(n, m, Tmp.GetNumericalValue());
	      }
	}
    }


  SparseRealMatrix* BMatrices = new SparseRealMatrix[this->NbrBMatrices];

  int** StartingIndexPerPLevel = new int* [this->PLevel + 1];
  int* TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  int* NbrIndicesPerPLevel = new int [this->PLevel + 1];
  TotalStartingIndexPerPLevel[0] = 0;
  StartingIndexPerPLevel[0] = new int [1];
  StartingIndexPerPLevel[0][0] = 0;

  
  int NbrNValue;
  int NValueShift;
  int QValue;
  int QValueDenominator;
  if ((this->RIndex & 1) == 0)
    {
      QValue = 1 + (this->RIndex / 2);
      NbrNValue = ((2 * this->PLevel) + QValue) + this->RIndex / 2 + 1;
      NValueShift = 2 * this->PLevel - 1;
      QValueDenominator = 1;
    }
  else
    {
      QValue = 2 + this->RIndex;
      NbrNValue = ((4 * this->PLevel) + QValue) + this->RIndex + 1;
      NValueShift = 4 * this->PLevel - 2;
      QValueDenominator = 2;
    }

     
  NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() + OrthogonalBasisPsiLeft[0].GetNbrColumn())) * NbrNValue;
  for (int i = 1; i <= this->PLevel; ++i)
    {
      TotalStartingIndexPerPLevel[i] = TotalStartingIndexPerPLevel[i - 1] + NbrIndicesPerPLevel[i - 1];
      StartingIndexPerPLevel[i] = new int [i + 1];      
      StartingIndexPerPLevel[i][0] = TotalStartingIndexPerPLevel[i];
      int Tmp = 0;
      int Tmp2;
      for (int j = 0; j < i; ++j)
	{
	  Tmp2 = U1BosonBasis[j]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[i - j].GetNbrColumn() + OrthogonalBasisPsiLeft[i - j].GetNbrColumn()) * NbrNValue;
	  StartingIndexPerPLevel[i][j + 1] = Tmp2 + StartingIndexPerPLevel[i][j];
	  Tmp += Tmp2;
	}
      Tmp += U1BosonBasis[i]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() + OrthogonalBasisPsiLeft[0].GetNbrColumn()) * NbrNValue;
      NbrIndicesPerPLevel[i] =  Tmp;
    }
  
  int MatrixSize = NbrIndicesPerPLevel[this->PLevel] + TotalStartingIndexPerPLevel[this->PLevel];
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
	      if ((this->RIndex & 1) == 0)
		{
		  for (int j = 1; j < NbrNValue; ++j)
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
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 1, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
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
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 1, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
		    }
		}
	      else
		{
		  for (int j = 2; j < NbrNValue; ++j)
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
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 2, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
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
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndex(j - 2, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							    this->Get2RMatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
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
				  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndex(N1, ChargedIndex1, NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->Get2RMatrixIndex(N2, ChargedIndex2, NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
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
				  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndex(N1, ChargedIndex1, NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->Get2RMatrixIndex(N2, ChargedIndex2, NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
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

