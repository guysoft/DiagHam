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
// qPeriodicity = if set to zero, guess it from the filling factor
// qSector = Q sector that has to be selected (from 0 to qPeriodicity-1)

FQHEMPSFixedQSectorMatrix::FQHEMPSFixedQSectorMatrix(AbstractFQHEMPSMatrix* matrix, int qPeriodicity, int qSector)
{
  this->MPSMatrix = matrix;
  this->QPeriodicity = qPeriodicity;
  this->QSector = qSector;
  this->PLevel = this->MPSMatrix->GetTruncationLevel();
  if (this->QPeriodicity == 0)
    {
      int TmpNumerator;
      this->MPSMatrix->GetFillingFactor(TmpNumerator, this->QPeriodicity);
    }
  int NbrGroupBMatrices = 1 << this->QPeriodicity;
  SparseRealMatrix* TmpSparseGroupBMatrices = new SparseRealMatrix[NbrGroupBMatrices];
  cout << "grouping " << this->QPeriodicity << " B matrices (" << NbrGroupBMatrices << " matrices)" << endl;
  for (int i = 0; i < NbrGroupBMatrices; ++i)
    {
      TmpSparseGroupBMatrices[i].Copy(this->MPSMatrix->GetMatrices()[i & 1]);
      for (int j = 1; j < this->QPeriodicity; ++j)
	TmpSparseGroupBMatrices[i].Multiply(this->MPSMatrix->GetMatrices()[(i >> j) & 1]);
    }
  int MinQ;
  int MaxQ;
  this->MPSMatrix->GetChargeIndexRange(MinQ, MaxQ);
  MinQ = this->QSector;
  int GroupBMatrixDimension = 0l;
  for (int p = 0; p <= this->PLevel; ++p)
    {
      for (int QValue = MinQ; QValue <= MaxQ; QValue += this->QPeriodicity)
	GroupBMatrixDimension += this->MPSMatrix->GetBondIndexRange(p, QValue);
    }
  int TmpBMatrixDimension = this->MPSMatrix->GetMatrices()[0].GetNbrRow();
  bool* GlobalIndices = new bool [TmpBMatrixDimension];
  for (int i = 0; i < TmpBMatrixDimension; ++i)
    GlobalIndices[i] = false;
  GroupBMatrixDimension = 0;
  for (int p = 0; p <= this->PLevel; ++p)
    {
      for (int QValue = MinQ; QValue <= MaxQ; QValue += this->QPeriodicity)
	{
	  int MaxLocalIndex = this->MPSMatrix->GetBondIndexRange(p, QValue);
	  for (int i = 0; i < MaxLocalIndex; ++i)
	    {
	      GlobalIndices[this->MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, p, QValue)] = true;
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

