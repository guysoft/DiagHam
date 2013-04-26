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
  this->UniformChargeIndexRange = true;
}

// constructor 
//
// laughlinIndex = power of the Laughlin part (i.e. 1/nu)
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// trimChargeIndices = trim the charge indices, assuming an iMPS
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSLaughlinMatrix::FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, bool trimChargeIndices, bool cylinderFlag, double kappa)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->CreateBMatrices();
}

// constructor from stored B matrices
//
// laughlinIndex = power of the Laughlin part (i.e. 1/nu)
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// trimChargeIndices = trim the charge indices, assuming an iMPS
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSLaughlinMatrix::FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, char* fileName, bool trimChargeIndices, bool cylinderFlag, double kappa)
{
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->LoadMatrices(fileName);
}

// destructor
//

FQHEMPSLaughlinMatrix::~FQHEMPSLaughlinMatrix()
{
  delete[] this->TotalStartingIndexPerPLevel;
  delete[] this->NbrIndicesPerPLevel;
  delete[] this->NbrNValuesPerPLevel;
  delete[] this->NInitialValuePerPLevel;
  delete[] this->NLastValuePerPLevel;
}
  
// get the name describing the B matrices 
// 
// return value = name 

char* FQHEMPSLaughlinMatrix::GetName()
{
  char* TmpName = new char[16];
  sprintf (TmpName, "laughlin%d", this->LaughlinIndex);
  return TmpName;
}

// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

void FQHEMPSLaughlinMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  numerator = 1;
  denominator = this->LaughlinIndex;
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
  this->NbrNValuesPerPLevel = new int [this->PLevel + 1];
  this->NInitialValuePerPLevel = new int [this->PLevel + 1];
  this->NLastValuePerPLevel = new int [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel[0] = 0;
  this->NbrNValue = ((2 * this->PLevel) + this->LaughlinIndex);
  int NValueShift = this->NbrNValue - 1;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->ComputeChargeIndexRange(i, this->NInitialValuePerPLevel[i], this->NLastValuePerPLevel[i], this->UniformChargeIndexRange);
      this->NbrNValuesPerPLevel[i] =  this->NLastValuePerPLevel[i] - this->NInitialValuePerPLevel[i] + 1;
    }
  this->NbrIndicesPerPLevel[0] = U1BosonBasis[0]->GetHilbertSpaceDimension() * this->NbrNValuesPerPLevel[0];
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      this->NbrIndicesPerPLevel[i] = U1BosonBasis[i]->GetHilbertSpaceDimension()  * this->NbrNValuesPerPLevel[i];
    }
  int MatrixSize = this->NbrIndicesPerPLevel[this->PLevel] + this->TotalStartingIndexPerPLevel[this->PLevel];
  cout << "B matrix size = " << MatrixSize << endl;
  int* TmpNbrElementPerRow = new int[MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = this->NInitialValuePerPLevel[i] + 1; j <= this->NLastValuePerPLevel[i]; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    ++TmpNbrElementPerRow[this->GetMatrixIndex(j - 1 - this->NInitialValuePerPLevel[i], k, this->NbrNValuesPerPLevel[i], this->TotalStartingIndexPerPLevel[i])];
	}
    }
  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = this->NInitialValuePerPLevel[i] + 1; j <= this->NLastValuePerPLevel[i]; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    {
	      double Tmp = 1.0;
	      if (this->CylinderFlag)
		Tmp *= exp(-this->Kappa * this->Kappa * (((double) i)
							 + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) / (4.0 * this->LaughlinIndex))
							 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) / (4.0 * this->LaughlinIndex))));
	      BMatrices[0].SetMatrixElement(this->GetMatrixIndex(j - 1 - this->NInitialValuePerPLevel[i], k, this->NbrNValuesPerPLevel[i], this->TotalStartingIndexPerPLevel[i]),
					    this->GetMatrixIndex(j - this->NInitialValuePerPLevel[i], k, this->NbrNValuesPerPLevel[i], this->TotalStartingIndexPerPLevel[i]), Tmp);
	    }
	}
    }

  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  FactorialCoefficient Coef;
  
  cout << "NValueShift = " << NValueShift << endl;
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
  	      if (((N1 >= this->NInitialValuePerPLevel[i]) && (N1 <= this->NLastValuePerPLevel[i]))
  		  && ((N2 >= this->NInitialValuePerPLevel[j]) && (N2 <= this->NLastValuePerPLevel[j])))
		{ 
		  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
		    {
		      for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
			{
			  ++TmpNbrElementPerRow[this->GetMatrixIndex(N1 - this->NInitialValuePerPLevel[i], k1, this->NbrNValuesPerPLevel[i], this->TotalStartingIndexPerPLevel[i])];
			}
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
  	      if (((N1 >= this->NInitialValuePerPLevel[i]) && (N1 <= this->NLastValuePerPLevel[i]))
  		  && ((N2 >= this->NInitialValuePerPLevel[j]) && (N2 <= this->NLastValuePerPLevel[j])))
		{ 
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
			  BMatrices[m].SetMatrixElement(this->GetMatrixIndex(N1 - this->NInitialValuePerPLevel[i], k1, this->NbrNValuesPerPLevel[i], this->TotalStartingIndexPerPLevel[i]), 
					 		this->GetMatrixIndex(N2 - this->NInitialValuePerPLevel[j], k2, this->NbrNValuesPerPLevel[j], this->TotalStartingIndexPerPLevel[j]), Tmp);
			}
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

  int NbrK1 = this->NbrIndicesPerPLevel[pLevel1] / this->NbrNValuesPerPLevel[pLevel1];
  int NbrK2 = this->NbrIndicesPerPLevel[pLevel2] / this->NbrNValuesPerPLevel[pLevel2];
  SparseRealMatrix TmpMatrix (NbrK1, NbrK2);
  for (int k1 = 0; k1 < NbrK1; ++k1)
    {
      for (int k2 = 0; k2 < NbrK2; ++k2)
	{
	  matrix.GetMatrixElement(this->GetMatrixIndex(q1 - this->NInitialValuePerPLevel[pLevel1], k1, this->NbrNValuesPerPLevel[pLevel1], this->TotalStartingIndexPerPLevel[pLevel1]),
                                  this->GetMatrixIndex(q2 - this->NInitialValuePerPLevel[pLevel2], k2, this->NbrNValuesPerPLevel[pLevel2], this->TotalStartingIndexPerPLevel[pLevel2]), Tmp);
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
  if ((pLevel < 0) || (pLevel > this->PLevel) || (qValue < this->NInitialValuePerPLevel[pLevel]) || (qValue > this->NLastValuePerPLevel[pLevel]))
    return 0;
  return this->NbrIndicesPerPLevel[pLevel] / this->NbrNValuesPerPLevel[pLevel];  
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int FQHEMPSLaughlinMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{  
  return (this->TotalStartingIndexPerPLevel[pLevel] + (localIndex * this->NbrNValuesPerPLevel[pLevel] + (qValue - this->NInitialValuePerPLevel[pLevel])));
}

// get the charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSLaughlinMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevel[pLevel];
  maxQ = this->NLastValuePerPLevel[pLevel];
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
  ReadLittleEndian(file, this->UniformChargeIndexRange);
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->NInitialValuePerPLevel = new int [this->PLevel + 1];
  this->NLastValuePerPLevel = new int [this->PLevel + 1];
  this->NbrNValuesPerPLevel = new int [this->PLevel + 1];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->TotalStartingIndexPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->NbrIndicesPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->NbrNValuesPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->NInitialValuePerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->NLastValuePerPLevel[i]);
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
  WriteLittleEndian(file, this->UniformChargeIndexRange);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->TotalStartingIndexPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->NbrIndicesPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->NbrNValuesPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->NInitialValuePerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->NLastValuePerPLevel[i]);
    }
  return true;
}

// compute the charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index
// uniformChargeIndexRange = if true, assume full range of indices; else, use trimming and variable range of indices

void FQHEMPSLaughlinMatrix::ComputeChargeIndexRange(int pLevel, int& minQ, int& maxQ, bool uniformChargeIndexRange)
{
  if (uniformChargeIndexRange == true)
    {
      minQ = 0;
      maxQ = this->NbrNValue - 1;
      return;
    }
//   if (this->PLevel == pLevel)
//     {
//       minQ = (this->NbrNValue - 5) / 2;
//       maxQ = (this->NbrNValue + 3) / 2;
//       return;
//     }
//   else
//     {
//       minQ = 0;
//       maxQ = this->NbrNValue - 1;
//       return;
//     }

  minQ = (this->NbrNValue - 1 - (this->LaughlinIndex - 1)) / 2;
  maxQ = (this->NbrNValue - 1 + (this->LaughlinIndex - 1)) / 2;
  if (minQ < 0)
    {
      minQ = 0;
    }
  if (maxQ >= this->NbrNValue)
    maxQ = this->NbrNValue - 1;
  int QShift = minQ;
  for (int LocalP = this->PLevel - 1; LocalP >= pLevel; --LocalP)
    {
      int LocalQ = minQ;
      while ((LocalQ >= 0) && 
	     (((LocalQ >= QShift) && ((((LocalQ - QShift) / (this->LaughlinIndex - 1)) * (2 * (LocalQ - QShift) - (this->LaughlinIndex - 1) * (((LocalQ - QShift) / (this->LaughlinIndex - 1)) + 1))) != (2 * (this->PLevel - LocalP)))) ||
	      ((LocalQ < QShift) && ((((LocalQ - QShift - (this->LaughlinIndex - 2)) / (this->LaughlinIndex - 1)) * (2 * (LocalQ - QShift) - (this->LaughlinIndex - 1) * (((LocalQ - QShift - (this->LaughlinIndex - 2)) / (this->LaughlinIndex - 1)) + 1))) != (2 * (this->PLevel - LocalP))))))

	{
	  --LocalQ;
	}
      if ((LocalQ >= 0) && (LocalQ < minQ))
	minQ = LocalQ;
      LocalQ = maxQ;
      while ((LocalQ < this->NbrNValue) && 
	     (((LocalQ >= QShift) && ((((LocalQ - QShift) / (this->LaughlinIndex - 1)) * (2 * (LocalQ - QShift) - (this->LaughlinIndex - 1) * (((LocalQ - QShift) / (this->LaughlinIndex - 1)) + 1))) != (2 * (this->PLevel - LocalP)))) ||
	      ((LocalQ < QShift) && ((((LocalQ - QShift - (this->LaughlinIndex - 2)) / (this->LaughlinIndex - 1)) * (2 * (LocalQ - QShift) - (this->LaughlinIndex - 1) * (((LocalQ - QShift - (this->LaughlinIndex - 2)) / (this->LaughlinIndex - 1)) + 1))) != (2 * (this->PLevel - LocalP))))))
	{
	  ++LocalQ;
	}
      if ((LocalQ < this->NbrNValue) && (LocalQ > maxQ))
	maxQ = LocalQ;
    }
//   minQ = (this->NbrNValue - 1 - 2 * (this->PLevel - pLevel + 1)) / 2;
//   maxQ = (this->NbrNValue - 1 + 2 * (this->PLevel - pLevel + 1)) / 2;
  cout << "range at " << pLevel << " : " << minQ << " " << maxQ << " (" << this->NbrNValue << ")" << endl;
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSLaughlinMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int MinQ;
  int MaxQ;
  this->GetChargeIndexRange(0, MinQ, MaxQ);
  if (padding == true)
    {
      rowIndex = this->PLevel + (this->LaughlinIndex - 1) / 2 - MinQ;
      columnIndex = rowIndex;
    }
  else
    {
      rowIndex = this->PLevel + this->LaughlinIndex - 1 - MinQ;
      columnIndex = this->PLevel - MinQ;
    }
}
