////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of MPS matrix for the Laughlin state                //
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
#include "Tools/FQHEMPS/FQHEMPSLaughlinMatrix.h"
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

FQHEMPSLaughlinMatrix::FQHEMPSLaughlinMatrix()
{
}

// constructor 
//
// laughlinIndex = power of the Laughlin part (i.e. 1/nu)
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSLaughlinMatrix::FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, bool cylinderFlag, double kappa)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->CreateBMatrices();
}

// constructor from stored B matrices
//
// laughlinIndex = power of the Laughlin part (i.e. 1/nu)
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSLaughlinMatrix::FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag, double kappa)
{
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->LoadMatrices(fileName);
}

// destructor
//

FQHEMPSLaughlinMatrix::~FQHEMPSLaughlinMatrix()
{
  delete[] this->TotalStartingIndexPerPLevel;
  delete[] this->NbrIndicesPerPLevel;
}
  
// create the B matrices for the laughlin state
//

void FQHEMPSLaughlinMatrix::CreateBMatrices ()
{
  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  SparseRealMatrix* BMatrices = new SparseRealMatrix[this->NbrBMatrices];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort(i, i, this->PLevel + 1);
    }
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel[0] = 0;
  this->NbrNValue = ((2 * this->PLevel) + this->LaughlinIndex);
  int NValueShift = this->NbrNValue - 1;
  this->NbrIndicesPerPLevel[0] = U1BosonBasis[0]->GetHilbertSpaceDimension() * this->NbrNValue;
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      this->NbrIndicesPerPLevel[i] = U1BosonBasis[i]->GetHilbertSpaceDimension()  * this->NbrNValue;
    }
  int MatrixSize = this->NbrIndicesPerPLevel[this->PLevel] + this->TotalStartingIndexPerPLevel[this->PLevel];
  cout << "B matrix size = " << MatrixSize << endl;
  int* TmpNbrElementPerRow = new int[MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = 1; j < this->NbrNValue; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    ++TmpNbrElementPerRow[this->GetMatrixIndex(j - 1, k, this->NbrNValue, this->TotalStartingIndexPerPLevel[i])];
	}
    }
  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = 1; j < this->NbrNValue; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    {
	      double Tmp = 1.0;
              if (this->CylinderFlag)
		Tmp *= exp(-this->Kappa * this->Kappa * (((double) i)
							 + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) / (4.0 * this->LaughlinIndex))
							 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) / (4.0 * this->LaughlinIndex))));
	      BMatrices[0].SetMatrixElement(this->GetMatrixIndex(j - 1, k, this->NbrNValue, this->TotalStartingIndexPerPLevel[i]),
					    this->GetMatrixIndex(j, k, this->NbrNValue, this->TotalStartingIndexPerPLevel[i]), Tmp);
	    }
	}
    }

  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  FactorialCoefficient Coef;

  for (int m = 1; m < this->NbrBMatrices; ++m)
    {
      for (int i = 0; i < MatrixSize; ++i)
	TmpNbrElementPerRow[i] = 0;
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
	      int N2 = (2 * (j - i) - this->LaughlinIndex + 1 + NValueShift) / 2;
	      int N1 = N2 + (this->LaughlinIndex - 1);
	      for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
		{
		  for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		    {
		      ++TmpNbrElementPerRow[this->GetMatrixIndex(N1, k1, this->NbrNValue, this->TotalStartingIndexPerPLevel[i])];
		    }
		}
	    }
	}

      BMatrices[m] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);

      for (int i = 0; i <= this->PLevel; ++i)
	{
	  BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
	      int N2 = (2 * (j - i) - this->LaughlinIndex + 1 + NValueShift) / 2;
	      int N1 = N2 + (this->LaughlinIndex - 1);
	      for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
		{
		  TmpSpace1->GetOccupationNumber(k1, Partition1);
		  for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		    {
		      TmpSpace2->GetOccupationNumber(k2, Partition2);
		      double Tmp = this->CreateLaughlinAMatrixElement(this->LaughlinIndex * m * m, 1, Partition1, Partition2, i, j, Coef);
		      if (this->CylinderFlag)
			Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (((double) (i + j))
								       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift)  / (2.0 * this->LaughlinIndex))
								       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift))  / (2.0 * this->LaughlinIndex))));
		      BMatrices[m].SetMatrixElement(this->GetMatrixIndex(N1, k1, this->NbrNValue, this->TotalStartingIndexPerPLevel[i]), 
						    this->GetMatrixIndex(N2, k2, this->NbrNValue, this->TotalStartingIndexPerPLevel[j]), Tmp);
		    }
		}
	    }
	}
    }
  
  delete[] Partition1;
  delete[] Partition2;

  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->RealBMatrices[i] = BMatrices[i];
    }
  delete[] BMatrices;
  for (int i = 0; i <= this->PLevel; ++i)
    delete U1BosonBasis[i];
  delete[] U1BosonBasis;
}

// create the matrix element of the B matrix U(1) part
//
// chargeNumerator = numerator of the charge (in sqrt(q) unit)
// chargeDenominator = denominator of the charge (in sqrt(q) unit)
// partition1 = U(1) partition associated to the left state
// p1Level = length of partition1
// partition2 = U(1) partition associated to the left state
// p1Level = length of partition2
// coef = reference on a temporary factorial coefficient
// return value = matrix element

double FQHEMPSLaughlinMatrix::CreateLaughlinAMatrixElement (int chargeNumerator, int chargeDenominator, 
							    unsigned long* partition1, unsigned long* partition2, 
							    int p1Level, int p2Level, FactorialCoefficient& coef)
{
  double Tmp = 1.0;
  int PMax = p1Level;
  if (p2Level > p1Level)
    PMax = p2Level;
  for (int i = 1; i <= PMax; ++i)
    {
      double Tmp2 = 0.0;
      for (int j = 0; j <= partition1[i]; ++j)
	{
	  int k = partition2[i] + j - partition1[i];
	  if ((k >= 0) && (k <= partition2[i]))
	    {
	      int Sum = k + j;
	      coef.SetToOne();
	      coef.PartialFactorialMultiply(partition1[i] - j + 1, partition1[i]);
	      coef.PartialFactorialMultiply(partition2[i] - k + 1, partition2[i]);
	      coef.FactorialDivide(j);
	      coef.FactorialDivide(k);
	      coef.FactorialDivide(j);
	      coef.FactorialDivide(k);
	      coef.PowerNMultiply(chargeNumerator, Sum);
	      coef.PowerNDivide(chargeDenominator, Sum);
	      coef.PowerNDivide(i, Sum);
	      Tmp2 += ((double) (1 - ((j & 1) << 1))) * sqrt(coef.GetNumericalValue());
	    }
	}
      Tmp *= Tmp2;
    }
  return Tmp;
}

// extract a block with fixed quantum numbers of a given matrix written the MPS basis
//
// matrix = reference on the matrix
// pLevel1 = tuncation level of the block left indices
// q1 = charge index of the block left indices
// pLevel1 = tuncation level of the block right indices
// q2 = charge index of the block left indices
// return value = block corresponding to the quantum numbers

SparseRealMatrix FQHEMPSLaughlinMatrix::ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2)
{
  double Tmp;

  int NbrK1 = this->NbrIndicesPerPLevel[pLevel1] / this->NbrNValue;
  int NbrK2 = this->NbrIndicesPerPLevel[pLevel1] / this->NbrNValue;
  SparseRealMatrix TmpMatrix (NbrK1, NbrK2);
  for (int k1 = 0; k1 < NbrK1; ++k1)
    {
      for (int k2 = 0; k2 < NbrK2; ++k2)
	{
	  matrix.GetMatrixElement(this->GetMatrixIndex(q1, k1, this->NbrNValue, this->TotalStartingIndexPerPLevel[pLevel1]), 
				  this->GetMatrixIndex(q2, k2, this->NbrNValue, this->TotalStartingIndexPerPLevel[pLevel2]), Tmp);
	  if (Tmp != 0.0)
	    TmpMatrix.SetMatrixElement(k1, k2, Tmp);
	}
    }

  return TmpMatrix;
}

// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSLaughlinMatrix::GetBondIndexRange(int pLevel, int qValue)
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

int FQHEMPSLaughlinMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{  
  return (this->TotalStartingIndexPerPLevel[pLevel] + (localIndex * this->NbrNValue + qValue));
}

// get the charge index range
// 
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSLaughlinMatrix::GetChargeIndexRange (int& minQ, int& maxQ)
{
  minQ = 0;
  maxQ = this->NbrNValue - 1;
  return;
}

// load the specific informations from the file header
// 
// file = reference on the input file stream
// return value = true if no error occurred  

bool FQHEMPSLaughlinMatrix::LoadHeader (ifstream& file)
{
  int HeaderSize = 0;
  ReadLittleEndian(file, HeaderSize);
  ReadLittleEndian(file, this->PLevel);
  ReadLittleEndian(file, this->LaughlinIndex);
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

bool FQHEMPSLaughlinMatrix::SaveHeader (ofstream& file)
{
  int HeaderSize = (this->PLevel + 1) * (2 * sizeof(int)) + (sizeof(int) * 3) + sizeof(bool) + sizeof(double);
  WriteLittleEndian(file, HeaderSize);
  WriteLittleEndian(file, this->PLevel);
  WriteLittleEndian(file, this->LaughlinIndex);
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

