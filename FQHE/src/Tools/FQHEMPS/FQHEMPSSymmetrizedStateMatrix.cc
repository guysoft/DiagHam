////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix built from symmetrized decoupled            //
//                          copies of a model state                           //
//                                                                            //
//                        last modification : 19/03/2013                      //
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
#include "Tools/FQHEMPS/FQHEMPSSymmetrizedStateMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "MathTools/IntegerAlgebraTools.h"
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

FQHEMPSSymmetrizedStateMatrix::FQHEMPSSymmetrizedStateMatrix()
{
}

// constructor from two MPS matrices
//
// matrix1 = MPS matrices that desribe the first state 
// matrix2 = MPS matrices that desribe the second state 
// antiSymmetrizeFlag - true if anti-symmetrization 

FQHEMPSSymmetrizedStateMatrix::FQHEMPSSymmetrizedStateMatrix(AbstractFQHEMPSMatrix* matrix1, AbstractFQHEMPSMatrix* matrix2, bool antiSymmetrizeFlag)
{
  this->MPSMatrix1 = matrix1;
  this->MPSMatrix2 = matrix2;
  this->PLevel = this->MPSMatrix1->GetTruncationLevel();
  this->NbrCFTSectors = 1;  
  this->TransferMatrixLargestEigenvalueDegeneracy = 1;  
  this->NbrBMatrices = 0;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  if (antiSymmetrizeFlag == false)
    {
      this->NbrBMatrices = this->MPSMatrix1->GetNbrMatrices();
      if (this->NbrBMatrices > this->MPSMatrix2->GetNbrMatrices())
	this->NbrBMatrices = this->MPSMatrix2->GetNbrMatrices();
      this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
      cout << "symmetrizing B matrices" << endl;
      this->RealBMatrices[0] = TensorProduct(this->MPSMatrix1->GetMatrices()[0], this->MPSMatrix2->GetMatrices()[0]);
      this->RealBMatrices[0] /= M_SQRT2;
      BinomialCoefficients Binomial(this->NbrBMatrices);
      for (int j = 1; j <  this->NbrBMatrices; ++j)
	{
	  this->RealBMatrices[j] = TensorProduct(this->MPSMatrix1->GetMatrices()[j], this->MPSMatrix2->GetMatrices()[0]);
	  for (int i = 1; i <= j; ++i)
	    {
	      SparseRealMatrix TmpMatrix = TensorProduct(this->MPSMatrix1->GetMatrices()[j - i], this->MPSMatrix2->GetMatrices()[i]);
	      TmpMatrix *= Binomial.GetNumericalCoefficient(j, i);
	      this->RealBMatrices[j] = this->RealBMatrices[j] + TmpMatrix;
	    }
	  this->RealBMatrices[j] /= M_SQRT2;
	}
    }
  else
    {
      this->NbrBMatrices = 2;
      this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
      cout << "symmetrizing B matrices" << endl;
      int* TmpNbrElementPerRow = new int[this->MPSMatrix1->GetMatrices()[0].GetNbrRow()];
      for (int i = 0; i < this->MPSMatrix1->GetMatrices()[0].GetNbrRow(); ++i)
	{
	  TmpNbrElementPerRow[i] = 1;
	}
      SparseRealMatrix SignMatrix (this->MPSMatrix1->GetMatrices()[0].GetNbrRow(), 
				   this->MPSMatrix1->GetMatrices()[0].GetNbrColumn(), TmpNbrElementPerRow);
      int TmpP;
      int TmpQ;
      for (int i = 0; i < SignMatrix.GetNbrRow(); ++i)
	{
	  this->MPSMatrix1->GetChargeAndPLevelFromMatrixIndex(i, TmpP, TmpQ);
	  cout << i << " " << TmpP << " " << TmpQ << endl;
	  if (((TmpQ / 3) & 1) == 0)
	    SignMatrix.SetMatrixElement(i, i, -1.0);
	  else
	    SignMatrix.SetMatrixElement(i, i, 1.0);
	}
      SparseRealMatrix SignedB0Matrix = Conjugate(SignMatrix, this->MPSMatrix1->GetMatrices()[0], SignMatrix);

      this->RealBMatrices[0] = TensorProduct(this->MPSMatrix1->GetMatrices()[0], this->MPSMatrix2->GetMatrices()[0]);
      this->RealBMatrices[1] = (TensorProduct(this->MPSMatrix1->GetMatrices()[1], this->MPSMatrix2->GetMatrices()[0])
				+ TensorProduct(SignedB0Matrix, this->MPSMatrix2->GetMatrices()[1]));
    }
  this->NbrNValuesPerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NInitialValuePerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NLastValuePerPLevelCFTSector = new int* [this->PLevel + 1];
  for (int p = 0; p <= this->PLevel; ++p)
    {
      this->NbrNValuesPerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      this->NInitialValuePerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      this->NLastValuePerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      for (int currentCFTSector = 0; currentCFTSector < this->NbrCFTSectors; ++currentCFTSector)
	{
	  int MinQ = 1 << 30;
	  int MaxQ = -(1 << 30);
	  for (int j = 0; j <= p; ++j)
	    {
	      int MinQ1 = 0;
	      int MaxQ1 = 0;
	      int MinQ2 = 0;
	      int MaxQ2 = 0;
	      this->MPSMatrix1->GetChargeIndexRange(j, currentCFTSector, MinQ1, MaxQ1);
	      this->MPSMatrix2->GetChargeIndexRange(p - j, currentCFTSector, MinQ2, MaxQ2);
	      if ((MinQ1 + MinQ2) < MinQ)
		{
		  MinQ = MinQ1 + MinQ2;
		}
	      if ((MaxQ1 + MaxQ2) > MaxQ)
		{
		  MaxQ = MaxQ1 + MaxQ2;
		}
	      this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] = MinQ;
	      this->NLastValuePerPLevelCFTSector[p][currentCFTSector] = MaxQ;
	      this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector] = this->NLastValuePerPLevelCFTSector[p][currentCFTSector] - this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] + 1;
	    } 
	}
    }
  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->PhysicalIndices[i] = (unsigned long) i;
    }
}

// create the B matrices for the block state
//

void FQHEMPSSymmetrizedStateMatrix::CreateBMatrices ()
{
}

// get the name describing the B matrices 
// 
// return value = name 

char* FQHEMPSSymmetrizedStateMatrix::GetName()
{
  char* TmpName1 = this->MPSMatrix1->GetName();
  char* TmpName2 = this->MPSMatrix2->GetName();
  char* TmpName3 = new char [strlen(TmpName1) + strlen(TmpName2) + 64];
  sprintf (TmpName3, "symmetrized_%s_%s", TmpName1, TmpName2);
  return TmpName3;
}


// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSSymmetrizedStateMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  cout << "warning : FQHEMPSSymmetrizedStateMatrix::GetBondIndexRange should not be used" << endl;
  return -1;  
}

// get the range for the bond index when fixing the tuncation level, charge and CFT sector index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = range for the bond index with fixed tuncation level, charge and CFT sector index

int FQHEMPSSymmetrizedStateMatrix::GetBondIndexRange(int pLevel, int qValue, int cftSector)
{
  if ((pLevel < 0) || (pLevel > this->PLevel) || (qValue < this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]) || 
      (qValue > this->NLastValuePerPLevelCFTSector[pLevel][cftSector]) || (cftSector > this->NbrCFTSectors) ||  (cftSector < 0))
    return 0;
  return this->NbrIndexPerPLevelCFTSectorQValue[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]];
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int FQHEMPSSymmetrizedStateMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  cout << "FQHEMPSSymmetrizedStateMatrix::GetBondIndexWithFixedChargeAndPLevel should not be used" << endl;
  return -1;
}

// get the bond index for a fixed truncation level, charge and CFT sector index
//
// localIndex = bond index in the pLevel and qValue and cftSector restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = bond index in the full bond index range

int FQHEMPSSymmetrizedStateMatrix::GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector)
{  
  return this->GlobalIndexMapper[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]][localIndex];
}

// get the charge index range
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSSymmetrizedStateMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][0];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][0];
  for (int i = 1; i < this->NbrCFTSectors; ++i)
    {
      if (this->NInitialValuePerPLevelCFTSector[pLevel][i] < minQ)
	minQ = this->NInitialValuePerPLevelCFTSector[pLevel][i];
      if (this->NLastValuePerPLevelCFTSector[pLevel][i] > maxQ)
	maxQ = this->NLastValuePerPLevelCFTSector[pLevel][i];  
    }
  return;
}

// get the charge index range at a given truncation level and in a given CFT sector
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSSymmetrizedStateMatrix::GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][cftSector];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][cftSector];
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSSymmetrizedStateMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int TmpRowIndex1;
  int TmpColumnIndex1;
  this->MPSMatrix1->GetMatrixBoundaryIndices(TmpRowIndex1, TmpColumnIndex1, padding);
  int TmpRowIndex2;
  int TmpColumnIndex2;
  this->MPSMatrix2->GetMatrixBoundaryIndices(TmpRowIndex2, TmpColumnIndex2, padding);
  ++TmpRowIndex1;
  ++TmpRowIndex2;
  ++TmpColumnIndex1;
  ++TmpColumnIndex2;
  rowIndex = TmpRowIndex1 + this->MPSMatrix1->GetMatrices()[0].GetNbrRow() * TmpRowIndex2;
  columnIndex = TmpColumnIndex1 + this->MPSMatrix1->GetMatrices()[0].GetNbrColumn() * TmpColumnIndex2;
}

