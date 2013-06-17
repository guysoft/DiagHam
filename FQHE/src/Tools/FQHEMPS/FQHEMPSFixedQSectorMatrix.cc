////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of MPS matrix built from a fixed charge sector of another MPS   //
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
#include "Tools/FQHEMPS/FQHEMPSFixedQSectorMatrix.h"
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

FQHEMPSFixedQSectorMatrix::FQHEMPSFixedQSectorMatrix()
{
}

// constructor from two MPS matrices (the number of B matrices has to be identical for all of them)
//
// matrix = MPS matrix
// qSector = Q sector that has to be selected (from 0 to qPeriodicity-1)
// qPeriodicity = if set to zero, guess it from the filling factor

FQHEMPSFixedQSectorMatrix::FQHEMPSFixedQSectorMatrix(AbstractFQHEMPSMatrix* matrix, int qSector, int qPeriodicity)
{
  this->MPSMatrix = matrix;
  this->QPeriodicity = qPeriodicity;
  this->QSector = qSector;
  this->PLevel = this->MPSMatrix->GetTruncationLevel();
  this->NbrCFTSectors = this->MPSMatrix->GetNbrCFTSectors();  
  this->TransferMatrixLargestEigenvalueDegeneracy = 1;  
  if (this->QPeriodicity == 0)
    {
      int TmpNumerator;
      this->MPSMatrix->GetFillingFactor(TmpNumerator, this->QPeriodicity);
      this->BMatrixGroupSize = this->QPeriodicity;
//       int GCD = FindGCD(this->QPeriodicity, TmpNumerator);
//       this->QPeriodicity /= GCD;
//      this->TransferMatrixLargestEigenvalueDegeneracy = GCD;
    }
  int NbrGroupBMatrices = 1 << this->BMatrixGroupSize;
  SparseRealMatrix* TmpSparseGroupBMatrices = new SparseRealMatrix[NbrGroupBMatrices];
  cout << "grouping " << this->BMatrixGroupSize << " B matrices (" << NbrGroupBMatrices << " matrices)" << endl;
  for (int i = 0; i < NbrGroupBMatrices; ++i)
    {
      TmpSparseGroupBMatrices[i].Copy(this->MPSMatrix->GetMatrices()[i & 1]);
      for (int j = 1; j < this->BMatrixGroupSize; ++j)
	TmpSparseGroupBMatrices[i].Multiply(this->MPSMatrix->GetMatrices()[(i >> j) & 1]);
    }
  int GroupBMatrixDimension = 0l;
  this->NbrNValuesPerPLevel = new int [this->PLevel + 1];
  this->NInitialValuePerPLevel = new int [this->PLevel + 1];
  this->NLastValuePerPLevel = new int [this->PLevel + 1];
  int MinQ;
  int MaxQ;
  int* QValueCFTSectorShift  = new int [this->NbrCFTSectors];
  QValueCFTSectorShift[0] = 0;
  QValueCFTSectorShift[1] = 2;

  for (int p = 0; p <= this->PLevel; ++p)
    {
      this->MPSMatrix->GetChargeIndexRange(p, MinQ, MaxQ);
      if ((MinQ % this->QPeriodicity) != 0)
	MinQ += this->QPeriodicity - (MinQ % this->QPeriodicity);
      //      MinQ += this->QSector;
      for (int l = 0; l < this->NbrCFTSectors; ++l)
	{
	  int LocalMinQ;
	  int LocalMaxQ;
	  this->MPSMatrix->GetChargeIndexRange(p, l, LocalMinQ, LocalMaxQ);
	  for (int QValue = MinQ; QValue <= MaxQ; QValue += this->QPeriodicity)
	    {
	      if (((QValue + QValueCFTSectorShift[l]) >= LocalMinQ) && ((QValue + QValueCFTSectorShift[l]) <= LocalMaxQ))
		GroupBMatrixDimension += this->MPSMatrix->GetBondIndexRange(p, QValue + QValueCFTSectorShift[l], 0);       
	    }
	}
      if (MaxQ >= MinQ)
	{
	  this->NInitialValuePerPLevel[p] = MinQ / this->QPeriodicity;
	  this->NLastValuePerPLevel[p] = MaxQ / this->QPeriodicity;
	  this->NbrNValuesPerPLevel[p] =  this->NLastValuePerPLevel[p] - this->NInitialValuePerPLevel[p] + 1;
	}
      else
	{
	  this->NInitialValuePerPLevel[p] = MinQ / this->QPeriodicity;
	  this->NLastValuePerPLevel[p] = MaxQ / this->QPeriodicity;
	  this->NbrNValuesPerPLevel[p] =  0;
	}

    }
  int TmpBMatrixDimension = this->MPSMatrix->GetMatrices()[0].GetNbrRow();
  
  int* GlobalIndices = new int [TmpBMatrixDimension];
  for (int i = 0; i < TmpBMatrixDimension; ++i)
    GlobalIndices[i] = -1;
  GroupBMatrixDimension = 0;
  this->NbrIndicesPerPLevelAndQSector = new int*[this->PLevel + 1];
  this->TotalStartingIndexPerPLevelAndQSector = new int*[this->PLevel + 1];
  this->GlobalIndexMapper = new int**[this->PLevel + 1];
  this->MPSMatrix->GetChargeIndexRange(0, MinQ, MaxQ);
  for (int p = 0; p <= this->PLevel; ++p)
    {
      this->NbrIndicesPerPLevelAndQSector[p] = new int[this->NbrNValuesPerPLevel[p]];
      this->TotalStartingIndexPerPLevelAndQSector[p] = new int[this->NbrNValuesPerPLevel[p]];
      this->GlobalIndexMapper[p] = new int*[this->NbrNValuesPerPLevel[p]];      
      this->MPSMatrix->GetChargeIndexRange(p, MinQ, MaxQ);
      if ((MinQ % this->QPeriodicity) != 0)
	MinQ += this->QPeriodicity - (MinQ % this->QPeriodicity);
      //      MinQ += this->QSector;
      int MinQ1, MaxQ1, MinQ2, MaxQ2;
      this->MPSMatrix->GetChargeIndexRange(p, 0, MinQ1, MaxQ1);
      this->MPSMatrix->GetChargeIndexRange(p, 1, MinQ2, MaxQ2);     
      for (int QValue = MinQ; QValue <= MaxQ; QValue += this->QPeriodicity)
	{
	  int MaxLocalIndex = 0;
	  if ((QValue <= MaxQ1) && (QValue >= MinQ1))
	    MaxLocalIndex = this->MPSMatrix->GetBondIndexRange(p, QValue, 0);
	  int MaxLocalIndex2 = 0;
	  if (((QValue + 2) <= MaxQ2) && ((QValue + 2) >= MinQ2))
	    MaxLocalIndex2 = this->MPSMatrix->GetBondIndexRange(p, QValue + 2, 1);
 	  this->GlobalIndexMapper[p][(QValue - MinQ) / this->QPeriodicity] = new int[MaxLocalIndex + MaxLocalIndex2];      
	  this->NbrIndicesPerPLevelAndQSector[p][(QValue - MinQ) / this->QPeriodicity] = MaxLocalIndex + MaxLocalIndex2;
	  this->TotalStartingIndexPerPLevelAndQSector[p][(QValue - MinQ) / this->QPeriodicity] = GroupBMatrixDimension;
	  for (int i = 0; i < MaxLocalIndex; ++i)
	    {
	      GlobalIndices[GroupBMatrixDimension] = this->MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(i, p, QValue, 0);
	      this->GlobalIndexMapper[p][(QValue - MinQ) / this->QPeriodicity][i] = GroupBMatrixDimension;      
	      ++GroupBMatrixDimension;
	    }
	  for (int i = 0; i < MaxLocalIndex2; ++i)
	    {
	      GlobalIndices[GroupBMatrixDimension] = this->MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(i, p, QValue + 2, 1);
	      this->GlobalIndexMapper[p][(QValue - MinQ) / this->QPeriodicity][MaxLocalIndex + i] = GroupBMatrixDimension;      
	      ++GroupBMatrixDimension;
	    }
	}
    }
  cout << "new B matrix dimension = " << GroupBMatrixDimension << endl;
  SparseRealMatrix* TmpSparseGroupBMatrices2 = new SparseRealMatrix[NbrGroupBMatrices];
  this->NbrBMatrices = 0;
  for (int i = 0; i < NbrGroupBMatrices; ++i)
    {
      TmpSparseGroupBMatrices2[i] = TmpSparseGroupBMatrices[i].ExtractMatrix(GroupBMatrixDimension, GroupBMatrixDimension, GlobalIndices, GlobalIndices);
      if (TmpSparseGroupBMatrices2[i].GetNbrRow() > 0)
 	++this->NbrBMatrices;
    }

  cout  << this->NbrBMatrices << " non zero B matrices" << endl;
  delete[] TmpSparseGroupBMatrices;
  delete[] GlobalIndices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->NbrBMatrices = 0;
  for (int i = 0; i < NbrGroupBMatrices; ++i)
    {
      if (TmpSparseGroupBMatrices2[i].GetNbrRow() > 0)
	{
	  this->RealBMatrices[this->NbrBMatrices] = TmpSparseGroupBMatrices2[i];
	  ++this->NbrBMatrices;
	}
      else
	{
	  cout << "throwing away B matrix " << i << endl;
	}
    }
}

// create the B matrices for the block state
//

void FQHEMPSFixedQSectorMatrix::CreateBMatrices ()
{
}

// get the name describing the B matrices 
// 
// return value = name 

char* FQHEMPSFixedQSectorMatrix::GetName()
{
  char* TmpName1 = this->MPSMatrix->GetName();
  char* TmpName2 = new char [strlen(TmpName1) + 64];
  sprintf (TmpName2, "%s_fixedq_%d_%d", TmpName1, this->QSector, this->QPeriodicity);
  return TmpName2;
}


// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSFixedQSectorMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  if ((pLevel < 0) || (pLevel > this->PLevel) || (qValue < this->NInitialValuePerPLevel[pLevel]) || (qValue > this->NLastValuePerPLevel[pLevel]))
    return 0;
  return this->NbrIndicesPerPLevelAndQSector[pLevel][qValue - this->NInitialValuePerPLevel[pLevel]];  
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int FQHEMPSFixedQSectorMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{  
  return this->GlobalIndexMapper[pLevel][qValue - this->NInitialValuePerPLevel[pLevel]][localIndex];
}

// get the charge index range
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSFixedQSectorMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevel[pLevel];
  maxQ = this->NLastValuePerPLevel[pLevel];
  return;
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSFixedQSectorMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int MinQ;
  int MaxQ;
  this->MPSMatrix->GetChargeIndexRange(0, MinQ, MaxQ);
  this->MPSMatrix->GetMatrixBoundaryIndices(rowIndex, columnIndex, padding);
  rowIndex += MinQ;
  columnIndex += MinQ;
  rowIndex /= this->QPeriodicity;
  columnIndex /= this->QPeriodicity;
  this->GetChargeIndexRange(0, MinQ, MaxQ);
  cout << MinQ << " " << MaxQ << endl;
  rowIndex -= MinQ;
  columnIndex -= MinQ;
  rowIndex = 1;
  columnIndex = 1;
}

