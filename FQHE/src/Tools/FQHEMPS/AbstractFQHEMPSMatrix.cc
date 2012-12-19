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
  return SparseRealMatrix();
}

// get the charge index range
// 
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void AbstractFQHEMPSMatrix::GetChargeIndexRange (int& minQ, int& maxQ)
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

// get the B matrices corresponding to localized quasiholes
//
// nbrQuasiholes = number of quasiholes
// quasiholePositions = quasihole positions
// return value = array of nbrQuasiholes matrices corresponding to each quasihole

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
// pLevel = tuncation level of the block
// qValue = charge index of the block
// localIndex = bond index in the pLevel and qValue restricted range
// return value = bond index in the full bond index range

int AbstractFQHEMPSMatrix::GetBondIndexWithFixedChargeAndPLevel(int pLevel, int qValue)
{
  return 0;
}

