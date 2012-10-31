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
}

// destructor
//

AbstractFQHEMPSMatrix::~AbstractFQHEMPSMatrix()
{
  if (this->RealBMatrices != 0)
    delete[] this->RealBMatrices;
  if (this->ComplexBMatrices != 0)
    delete[] this->ComplexBMatrices;
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

