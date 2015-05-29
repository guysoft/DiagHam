////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix built from symmetrized decoupled            //
//                   copies of a model state, using a twist                   //
//                                                                            //
//                        last modification : 20/05/2015                      //
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
#include "Tools/FQHEMPS/FQHEMPSTwistedSymmetrizedStateMatrix.h"
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

FQHEMPSTwistedSymmetrizedStateMatrix::FQHEMPSTwistedSymmetrizedStateMatrix()
{
}

// constructor from the MPS matrix on the doubled system
//
// matrix = MPS matrices that describe the state 
// antiSymmetrizeFlag - true if anti-symmetrization 

FQHEMPSTwistedSymmetrizedStateMatrix::FQHEMPSTwistedSymmetrizedStateMatrix(AbstractFQHEMPSMatrix* matrix, int shift, bool antiSymmetrizeFlag)
{
  this->MPSMatrix = matrix;
  this->PLevel = this->MPSMatrix->GetTruncationLevel();
  this->NbrCFTSectors = 1;  
  this->TransferMatrixLargestEigenvalueDegeneracy = 1;  
  this->NbrBMatrices = 0;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  if (antiSymmetrizeFlag == false)
    {
      this->NbrBMatrices = this->MPSMatrix->GetNbrMatrices();
      this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
      cout << "symmetrizing B matrices with a twist" << endl;
      if (shift == 0)
	{
	  for (int j = 0; j < this->NbrBMatrices; ++j)
	    {
	      this->RealBMatrices[j] = Multiply(this->MPSMatrix->GetMatrices()[j], this->MPSMatrix->GetMatrices()[0]);
	    }
	}
      else
	{
	  for (int j = 0; j < this->NbrBMatrices; ++j)
	    {
	      this->RealBMatrices[j] = Multiply(this->MPSMatrix->GetMatrices()[0], this->MPSMatrix->GetMatrices()[j]);
	    }
	}
    }
  else
    {
      this->NbrBMatrices = 2;
      this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
      cout << "antisymmetrizing B matrices with a twist" << endl;
      if (shift == 0)
	{
	  this->RealBMatrices[0] = Multiply(this->MPSMatrix->GetMatrices()[0], this->MPSMatrix->GetMatrices()[0]);
	  this->RealBMatrices[1] = Multiply(this->MPSMatrix->GetMatrices()[1], this->MPSMatrix->GetMatrices()[0]);
	}
      else
	{
	  this->RealBMatrices[0] = Multiply(this->MPSMatrix->GetMatrices()[0], this->MPSMatrix->GetMatrices()[0]);
	  this->RealBMatrices[1] = Multiply(this->MPSMatrix->GetMatrices()[0], this->MPSMatrix->GetMatrices()[1]);
	}
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
	      this->MPSMatrix->GetChargeIndexRange(j, currentCFTSector, MinQ1, MaxQ1);
	      if (MinQ1 < MinQ)
		{
		  MinQ = MinQ1;
		}
	      if (MaxQ1 > MaxQ)
		{
		  MaxQ = MaxQ1;
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

void FQHEMPSTwistedSymmetrizedStateMatrix::CreateBMatrices ()
{
}

// get the name describing the B matrices 
// 
// return value = name 

char* FQHEMPSTwistedSymmetrizedStateMatrix::GetName()
{
  char* TmpName1 = this->MPSMatrix->GetName();
  char* TmpName2 = new char [strlen(TmpName1) + 64];
  sprintf (TmpName2, "twisted_symmetrized_%s", TmpName1);
  return TmpName2;
}


// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSTwistedSymmetrizedStateMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  cout << "warning : FQHEMPSTwistedSymmetrizedStateMatrix::GetBondIndexRange should not be used" << endl;
  return -1;  
}

// get the range for the bond index when fixing the tuncation level, charge and CFT sector index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = range for the bond index with fixed tuncation level, charge and CFT sector index

int FQHEMPSTwistedSymmetrizedStateMatrix::GetBondIndexRange(int pLevel, int qValue, int cftSector)
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

int FQHEMPSTwistedSymmetrizedStateMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  cout << "FQHEMPSTwistedSymmetrizedStateMatrix::GetBondIndexWithFixedChargeAndPLevel should not be used" << endl;
  return -1;
}

// get the bond index for a fixed truncation level, charge and CFT sector index
//
// localIndex = bond index in the pLevel and qValue and cftSector restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = bond index in the full bond index range

int FQHEMPSTwistedSymmetrizedStateMatrix::GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector)
{  
  return this->GlobalIndexMapper[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]][localIndex];
}

// get the charge index range
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSTwistedSymmetrizedStateMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
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

void FQHEMPSTwistedSymmetrizedStateMatrix::GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][cftSector];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][cftSector];
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSTwistedSymmetrizedStateMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  this->MPSMatrix->GetMatrixBoundaryIndices(rowIndex, columnIndex, padding);
}
