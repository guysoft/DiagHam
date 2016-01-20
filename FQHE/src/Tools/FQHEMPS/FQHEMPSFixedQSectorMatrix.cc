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
#include "GeneralTools/ArrayTools.h"

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
    }
  this->TorusFlag = matrix->IsTorus();
  int NbrBMatricesPerOrbital = matrix->GetNbrMatrices();
  int NbrGroupBMatrices = 1;
  for (int i = 0; i < this->BMatrixGroupSize; ++i)
    NbrGroupBMatrices *= NbrBMatricesPerOrbital;
  SparseRealMatrix* TmpSparseGroupBMatrices = new SparseRealMatrix[NbrGroupBMatrices];
  unsigned long** TmpPhysicalIndices = new unsigned long*[NbrGroupBMatrices];
  for (int i = 0; i < NbrGroupBMatrices; ++i)
    TmpPhysicalIndices[i] = new unsigned long [this->BMatrixGroupSize];
  cout << "grouping " << this->BMatrixGroupSize << " B matrices (" << NbrGroupBMatrices << " matrices)" << endl;
  
  int Step = NbrGroupBMatrices / NbrBMatricesPerOrbital;
  int TmpOrbitalIndex = this->BMatrixGroupSize - 1;
  unsigned long* TmpPhysicalIndex = new unsigned long[this->MPSMatrix->GetNbrOrbitals()];
  for (int i = 0; i < NbrGroupBMatrices; i += Step)
    {
      TmpSparseGroupBMatrices[i].Copy(this->MPSMatrix->GetMatrices()[i / Step]);
      this->MPSMatrix->GetPhysicalIndex(i / Step, TmpPhysicalIndex);
      TmpPhysicalIndices[i][TmpOrbitalIndex] = TmpPhysicalIndex[0];
    }
  --TmpOrbitalIndex;
  while (Step > 1)
    {
      int TmpStep = Step / NbrBMatricesPerOrbital;
      for (int i = 0; i < NbrGroupBMatrices; i += Step)
	{
	  for (int j = 1; j < NbrBMatricesPerOrbital; ++j)
	    {
	      TmpSparseGroupBMatrices[i + j * TmpStep].Copy(TmpSparseGroupBMatrices[i]);
	      for (int k = this->BMatrixGroupSize - 1; k > TmpOrbitalIndex; --k)
		TmpPhysicalIndices[i + j * TmpStep][k] = TmpPhysicalIndices[i][k];
	    }
 	  for (int j = 0; j < NbrBMatricesPerOrbital; ++j)
	    {
	      TmpSparseGroupBMatrices[i + j * TmpStep].Multiply(this->MPSMatrix->GetMatrices()[j]);
	      this->MPSMatrix->GetPhysicalIndex(j, TmpPhysicalIndex);
	      TmpPhysicalIndices[i + j * TmpStep][TmpOrbitalIndex] =  TmpPhysicalIndex[0];
	    }	  
	}
      Step = TmpStep;
      --TmpOrbitalIndex;
    }
  delete[] TmpPhysicalIndex;

  int GroupBMatrixDimension = 0l;
  int MinQ;
  int MaxQ;
  int* QValueCFTSectorShift  = new int [this->NbrCFTSectors];
  for (int i = 0; i < this->NbrCFTSectors; ++i)
    {
      QValueCFTSectorShift[i] = this->MPSMatrix->GetQValueCFTSectorShift(i) % this->QPeriodicity;
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
	  int LocalMinQ;
	  int LocalMaxQ;
	  this->MPSMatrix->GetChargeIndexRange(p, currentCFTSector, LocalMinQ, LocalMaxQ);
	  MinQ = (this->QSector + QValueCFTSectorShift[currentCFTSector]) % this->QPeriodicity;
	  while (MinQ < LocalMinQ)
	    MinQ += this->QPeriodicity;
	  int QValue = MinQ;
	  MaxQ = LocalMaxQ;
	  for (; QValue <= MaxQ; QValue += this->QPeriodicity)
	    {
	      GroupBMatrixDimension += this->MPSMatrix->GetBondIndexRange(p, QValue, currentCFTSector);       
	    }
	  QValue -= this->QPeriodicity;
	  MaxQ = QValue;
	  if (MaxQ >= MinQ)
	    {
	      this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] = (MinQ - QValueCFTSectorShift[currentCFTSector]) / this->QPeriodicity;
	      this->NLastValuePerPLevelCFTSector[p][currentCFTSector] = (MaxQ - QValueCFTSectorShift[currentCFTSector]) / this->QPeriodicity;
	      this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector] = this->NLastValuePerPLevelCFTSector[p][currentCFTSector] - this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] + 1;
	    }
	  else
	    {
	      this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] = 0;
	      this->NLastValuePerPLevelCFTSector[p][currentCFTSector] = -1;
	      this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector] = 0;
	    }
	}
    }
  int TmpBMatrixDimension = this->MPSMatrix->GetMatrices()[0].GetNbrRow();
  
  this->GlobalIndices = new int [TmpBMatrixDimension];
  for (int i = 0; i < TmpBMatrixDimension; ++i)
    this->GlobalIndices[i] = -1;
  GroupBMatrixDimension = 0;
  this->NbrIndexPerPLevelCFTSectorQValue = new int**[this->PLevel + 1];
  this->StartingIndexPerPLevelCFTSectorQValue = new int**[this->PLevel + 1];
  this->GlobalIndexMapper = new int***[this->PLevel + 1];
  for (int p = 0; p <= this->PLevel; ++p)
    {
      this->NbrIndexPerPLevelCFTSectorQValue[p] = new int*[this->NbrCFTSectors];
      this->StartingIndexPerPLevelCFTSectorQValue[p] = new int*[this->NbrCFTSectors];
      this->GlobalIndexMapper[p] = new int**[this->NbrCFTSectors];      
      for (int currentCFTSector = 0; currentCFTSector < this->NbrCFTSectors; ++currentCFTSector)
	{
	  this->NbrIndexPerPLevelCFTSectorQValue[p][currentCFTSector] = new int[this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector]];
	  this->StartingIndexPerPLevelCFTSectorQValue[p][currentCFTSector] = new int[this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector]];
	  this->GlobalIndexMapper[p][currentCFTSector] = new int*[this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector]];
	  int LocalMinQ;
	  int LocalMaxQ;
	  this->MPSMatrix->GetChargeIndexRange(p, currentCFTSector, LocalMinQ, LocalMaxQ);
	  MinQ = (this->QSector + QValueCFTSectorShift[currentCFTSector]) % this->QPeriodicity;
	  while (MinQ < LocalMinQ)
	    MinQ += this->QPeriodicity;
	  int QValue = MinQ;
	  MaxQ = LocalMaxQ;
	  for (; QValue <= MaxQ; QValue += this->QPeriodicity)
	    {
	      int LocalQSector = (QValue - MinQ + QValueCFTSectorShift[currentCFTSector]) / this->QPeriodicity;
	      int MaxLocalIndex = this->MPSMatrix->GetBondIndexRange(p, QValue, currentCFTSector);
	      this->GlobalIndexMapper[p][currentCFTSector][LocalQSector] = new int[MaxLocalIndex];      
	      this->NbrIndexPerPLevelCFTSectorQValue[p][currentCFTSector][LocalQSector] = MaxLocalIndex;
	      this->StartingIndexPerPLevelCFTSectorQValue[p][currentCFTSector][LocalQSector] = GroupBMatrixDimension;
	      for (int i = 0; i < MaxLocalIndex; ++i)
		{
		  this->GlobalIndices[GroupBMatrixDimension] = this->MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(i, p, QValue, currentCFTSector);
		  this->GlobalIndexMapper[p][currentCFTSector][LocalQSector][i] = GroupBMatrixDimension;      
		  ++GroupBMatrixDimension;
		}
	    }
	}
    }
  cout << "new B matrix dimension = " << GroupBMatrixDimension << endl;
  SparseRealMatrix* TmpSparseGroupBMatrices2 = new SparseRealMatrix[NbrGroupBMatrices];
  this->NbrBMatrices = 0;
  for (int i = 0; i < NbrGroupBMatrices; ++i)
    {
      TmpSparseGroupBMatrices2[i] = TmpSparseGroupBMatrices[i].ExtractMatrix(GroupBMatrixDimension, GroupBMatrixDimension, this->GlobalIndices, this->GlobalIndices);
      if (TmpSparseGroupBMatrices2[i].GetNbrRow() > 0)
 	++this->NbrBMatrices;
    }

  if (this->TorusFlag == true)
    {
      this->TopologicalSectorIndices = 0;
      this->TopologicalSectorNbrIndices = 0;
      for (int t = 0; t < 2; ++t)
	{
	  int TmpTopologicalSectorNbrIndices = 0;
	  int* TmpTopologicalSectorIndices = this->MPSMatrix->GetTopologicalSectorIndices(t, TmpTopologicalSectorNbrIndices);
	  if (TmpTopologicalSectorNbrIndices > 0)
	    {
	      this->TopologicalSectorNbrIndices = 0;
	      this->TopologicalSectorIndices = new int[TmpTopologicalSectorNbrIndices];
	      for (int i = 0; i < TmpTopologicalSectorNbrIndices; ++i)
		{
		  int TmpPos = SearchInArray(TmpTopologicalSectorIndices[i], this->GlobalIndices, GroupBMatrixDimension);
		  if (TmpPos >= 0)
		    {
		      this->TopologicalSectorIndices[this->TopologicalSectorNbrIndices] = TmpPos;
		      ++this->TopologicalSectorNbrIndices;
		    }
		}
	      if (this->TopologicalSectorNbrIndices == 0)
		{
		  delete[] this->TopologicalSectorIndices;
		  this->TopologicalSectorIndices = 0;
		}
	      else
		{
		  t =2;
		}
	      delete[] TmpTopologicalSectorIndices;
	    }
	}
    }
  else
    {
      this->TopologicalSectorIndices = 0;
      this->TopologicalSectorNbrIndices = 0;
    }
  cout  << this->NbrBMatrices << " non zero B matrices" << endl;
  delete[] TmpSparseGroupBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
  this->NbrBMatrices = 0;
  this->PhysicalIndices[1] = 0x1ul;
  for (int i = 0; i < NbrGroupBMatrices; ++i)
    {
      if (TmpSparseGroupBMatrices2[i].GetNbrRow() > 0)
	{
	  this->RealBMatrices[this->NbrBMatrices] = TmpSparseGroupBMatrices2[i];
	  this->PhysicalIndices[this->NbrBMatrices] = 0x0ul;
	  if (this->GetMaximumOccupation() == 1)
	    {
	      for (int j = 0; j < this->BMatrixGroupSize; ++j)
		this->PhysicalIndices[this->NbrBMatrices] |= TmpPhysicalIndices[i][j] << j;
	    }
	  else
	    {
	      int TmpPos = 0;
	      for (int j = 0; j < this->BMatrixGroupSize; ++j)
		{
		  this->PhysicalIndices[this->NbrBMatrices] |= ((0x1ul << TmpPhysicalIndices[i][j]) - 0x1ul) << TmpPos;
		  TmpPos += TmpPhysicalIndices[i][j];
		  TmpPos++;
		}
	      
	    }
	  ++this->NbrBMatrices;
	}
      else
	{
	  cout << "throwing away B matrix " << i << endl;
	}
    }
  for (int i = 0; i < NbrGroupBMatrices; ++i)
    delete[] TmpPhysicalIndices[i];
  delete[] TmpPhysicalIndices;
}

  // destructor
  //

FQHEMPSFixedQSectorMatrix::~FQHEMPSFixedQSectorMatrix()
{
  if ((this->TorusFlag == true) && (this->TopologicalSectorIndices != 0))
    {
      delete[] this->TopologicalSectorIndices;
      if (this->GlobalIndices != 0)
	delete[] this->GlobalIndices;
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
  cout << "warning : FQHEMPSFixedQSectorMatrix::GetBondIndexRange should not be used" << endl;
  return -1;  
}

// get the range for the bond index when fixing the tuncation level, charge and CFT sector index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = range for the bond index with fixed tuncation level, charge and CFT sector index

int FQHEMPSFixedQSectorMatrix::GetBondIndexRange(int pLevel, int qValue, int cftSector)
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

int FQHEMPSFixedQSectorMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  cout << "FQHEMPSFixedQSectorMatrix::GetBondIndexWithFixedChargeAndPLevel should not be used" << endl;
  return -1;
}

// get the bond index for a fixed truncation level, charge and CFT sector index
//
// localIndex = bond index in the pLevel and qValue and cftSector restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = bond index in the full bond index range

int FQHEMPSFixedQSectorMatrix::GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector)
{  
  return this->GlobalIndexMapper[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]][localIndex];
}

// get the charge index range
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSFixedQSectorMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
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

void FQHEMPSFixedQSectorMatrix::GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][cftSector];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][cftSector];
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

// get a given physical indiex
//
// index = index to retrieve
// configuration = array where the description of the physical index will be stored

void FQHEMPSFixedQSectorMatrix::GetPhysicalIndex(int index, unsigned long* configuration)
{
  if (this->GetMaximumOccupation() == 1)
    {
      for (int i = 0; i < this->GetNbrOrbitals(); ++i)
	{
	  configuration[i] = (this->PhysicalIndices[index] >> i) & 0x1ul;
	}
    }
  else
    {
      unsigned long InitialState = this->PhysicalIndices[index];
      for (int i = 0; i < this->GetNbrOrbitals(); ++i)
	{
	  unsigned long TmpState = (~InitialState - 1ul) ^ (~InitialState);
	  TmpState &= ~(TmpState >> 1);
#ifdef  __64_BITS__
	  unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
	  TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
	  TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
	  TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
	  TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
	  TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
	  unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
	  TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
	  TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
	  TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
	  TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
	  configuration[i] = (unsigned long) TmpPower;
	  ++TmpPower;
	  InitialState >>= TmpPower;
	}
    }
}

