////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract MPS matrix for the FQHE                 //
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
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"
#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor 
//

AbstractFQHEMPSMatrix::AbstractFQHEMPSMatrix()
{
  this->NbrBMatrices = 0;
  this->RealBMatrices = 0;
  this->ComplexBMatrices = 0;
  this->QuasiholeBMatrices = 0;
}

// destructor
//

AbstractFQHEMPSMatrix::~AbstractFQHEMPSMatrix()
{
  if (this->RealBMatrices != 0)
    delete[] this->RealBMatrices;
  if (this->ComplexBMatrices != 0)
    delete[] this->ComplexBMatrices;
  if (this->QuasiholeBMatrices != 0)
    delete[] this->QuasiholeBMatrices;
}
  
// save the matrices 
// 
// fileName = name of the file where the matrices have to be stored
// return value = true if no error occurred  

bool AbstractFQHEMPSMatrix::SaveMatrices (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->NbrBMatrices);
  this->SaveHeader(File);
  if (this->RealBMatrices != 0)
    {
      int ComplexFlag  = 0;
      WriteLittleEndian(File, ComplexFlag);
      for (int i = 0; i < this->NbrBMatrices; ++i)
	{
	  if (this->RealBMatrices[i].WriteMatrix(File) == false)
	    {
	      File.close();  
	      cout << "error while storing matrix " << i << ", can't create " << fileName << endl;
	      return false;
	    }
	}
    }
  File.close();  
  return true;
}

// load the matrices 
// 
// fileName = name of the file where the matrices are stored
// return value = true if no error occurred  

bool AbstractFQHEMPSMatrix::LoadMatrices (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  ReadLittleEndian(File, this->NbrBMatrices);
  this->LoadHeader(File);
  int ComplexFlag  = 0;
  ReadLittleEndian(File, ComplexFlag);
  if (ComplexFlag == 0)
    {
      this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
      for (int i = 0; i < this->NbrBMatrices; ++i)
	{
	  if (this->RealBMatrices[i].ReadMatrix(File) == false)
	    {
	      File.close();  
	      cout << "error while reading matrix " << i << ", can't read " << fileName << endl;
	      return false;
	    }
	}
    }
  File.close();  
  return true;
}

// get the name describing the B matrices 
// 
// return value = name 

char* AbstractFQHEMPSMatrix::GetName()
{
  char* TmpName = new char[16];
  sprintf(TmpName, "dummy");
  return TmpName;
}

// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

void AbstractFQHEMPSMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  numerator = 0;
  denominator = 0;
}

// extract a block with fixed quantum numbers of a given matrix written the MPS basis
//
// matrix = reference on the matrix
// pLevel1 = tuncation level of the block left indices
// q1 = charge index of the block left indices
// pLevel1 = tuncation level of the block right indices
// q2 = charge index of the block left indices
// return value = block corresponding to the quantum numbers

SparseRealMatrix AbstractFQHEMPSMatrix::ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2)
{
  int BlockNbrRow = this->GetBondIndexRange(pLevel1, q1);
  int BlockNbrColumn = this->GetBondIndexRange(pLevel2, q2);
  SparseRealMatrix TmpMatrix(BlockNbrRow, BlockNbrColumn);
  double Tmp = 0.0;
  for (int i = 0; i < BlockNbrRow; ++i)
    {
      for (int j = 0; j < BlockNbrColumn; ++j)
	{
	  matrix.GetMatrixElement(this->GetBondIndexWithFixedChargeAndPLevel(i, pLevel1, q1),
				  this->GetBondIndexWithFixedChargeAndPLevel(j, pLevel2, q2), Tmp);
	  if (Tmp != 0.0)
	    {
	      TmpMatrix.SetMatrixElement(i, j, Tmp);
	    }
	}
    }
  return TmpMatrix;
}

// get the charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void AbstractFQHEMPSMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
{
  minQ = 1;
  maxQ = 0;
  return;
}

// compute the global charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void AbstractFQHEMPSMatrix::ComputeGlobalChargeIndexRange(int pLevel, int& minQ, int& maxQ)
{
  minQ = 1;
  maxQ = 0;
  return;
}

// load the specific informations from the file header
// 
// file = reference on the input file stream
// return value = true if no error occurred  

bool AbstractFQHEMPSMatrix::LoadHeader (ifstream& file)
{
  int HeaderSize = 0;
  ReadLittleEndian(file, HeaderSize);
  file.seekg (HeaderSize, ios::cur);
  return true;
}

// save the specific informations to the file header 
// 
// file = reference on the output file stream
// return value = true if no error occurred  

bool AbstractFQHEMPSMatrix::SaveHeader (ofstream& file)
{
  int HeaderSize = 0;
  WriteLittleEndian(file, HeaderSize);
  return true;
}

// get the edge matrix for localized quasiholes, with normal ordering
//
// nbrQuasiholes = number of quasiholes
// quasiholePositions = quasihole positions (for cylinder, positions have to be expressed in perimeter units)
// return value = pointer to the edge matrix

SparseComplexMatrix* AbstractFQHEMPSMatrix::GetQuasiholeMatrices(int nbrQuasiholes, Complex* quasiholePositions)
{
  return 0;
}
  
// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int AbstractFQHEMPSMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  return 0;
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int AbstractFQHEMPSMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  return 0;
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void AbstractFQHEMPSMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  rowIndex = -1;
  columnIndex = -1;
}

