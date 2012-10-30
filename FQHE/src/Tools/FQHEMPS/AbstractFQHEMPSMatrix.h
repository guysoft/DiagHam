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


#ifndef ABSTRACTFQHEMPSMATRIX_H
#define ABSTRACTFQHEMPSMATRIX_H


#include "config.h"


class SparseRealMatrix;
class SparseComplexMatrix;


class AbstractFQHEMPSMatrix
{

 protected:

  // number of B matrices
  int NbrBMatrices;

  // array where the B matrices are stored (for real matrices)
  SparseRealMatrix* RealBMatrices;

  // array where the B matrices are stored (for complex matrices)
  SparseComplexMatrix* ComplexBMatrices;

 public:
  
  // default constructor 
  //
  AbstractFQHEMPSMatrix();

  // destructor
  //
  ~AbstractFQHEMPSMatrix();
  
  // save the matrices 
  // 
  // fileName = name of the file where the matrices have to be stored
  // return value = true if no error occurred  
  virtual bool SaveMatrices (char* fileName);

  // load the matrices 
  // 
  // fileName = name of the file where the matrices are stored
  // return value = true if no error occurred  
  virtual bool LoadMatrices (char* fileName);

  // get the array where the matrices are stored
  //
  // return value = pointer to the array
  virtual SparseRealMatrix* GetMatrices();

 protected:

};

// get the array where the matrices are stored
//
// return value = pointer to the array

inline SparseRealMatrix* AbstractFQHEMPSMatrix::GetMatrices()
{
  return this->RealBMatrices;
}

#endif
