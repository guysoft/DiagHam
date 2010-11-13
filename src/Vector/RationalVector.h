////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class for n dimensional rational vector                   //
//                                                                            //
//                        last modification : 04/11/2010                      //
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


#ifndef RATIONALVECTOR_H
#define RATIONALVECTOR_H


#include "config.h"
#include "Vector/Vector.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>

#ifdef __MPI__
#include <mpi.h>
#endif


using std::ostream;


class RationalVector : public Vector
{

 protected:

  // array for the numerators
  long* Numerators;
  // array for the denominators
  long* Denominators;

  // garbage flag to avoid data duplication
  GarbageFlag Flag;

 public:

  // default constructor  
  // 
  RationalVector();

  // default constructor for an empty vector
  // 
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  RationalVector(int size, bool zeroFlag = false);

  // constructor for an empty rational vector bigger than 2^31
  //
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  RationalVector(long size, bool zeroFlag = false);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  RationalVector(const RationalVector& vector, bool duplicateFlag = false);

  // destructor
  //
  ~RationalVector ();

  // Resize vector
  //
  // dimension = new dimension
  virtual void Resize (int dimension);

  // Resize vector and set to zero all components that have been added
  //
  // dimension = new dimension
  virtual void ResizeAndClean (int dimension);

  // put all vector components to zero
  //
  // return value = reference on current vector
  virtual Vector& ClearVector ();

  // put select vector components to zero
  // start = start index
  // nbrComponent = number of components to set to zero
  // return value = reference on current vector
  virtual Vector& ClearVectorSegment (long start, long nbrComponent);

  // create a new vector with same size and same type but non-initialized components
  //
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  virtual Vector* EmptyClone(bool zeroFlag = false);

  // create an array of new vectors with same size and same type but non-initialized components
  //
  // nbrVectors = number of vectors to sreate
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to the array of new vectors
  virtual Vector* EmptyCloneArray(int nbrVectors, bool zeroFlag = false);

  // get the numerator of a given component
  // 
  // index = component index
  // return value = reference on the denominator
  long& Num(int index);

  // get the numerator of a given component
  // 
  // index = component index
  // return value = reference on the denominator
  long& Num(long index);

  // get the denominator of a given component
  // 
  // index = component index
  // return value = reference on the denominator
  long& Den(int index);

  // get the denominator of a given component
  // 
  // index = component index
  // return value = reference on the denominator
  long& Den(long index);

  // write vector in a file 
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteVector (const char* fileName);

  // write vector in a file in ascii mode
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteAsciiVector (const char* fileName);

  // read vector from a file 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if no error occurs
  bool ReadVector (const char* fileName);

  // read vector from a file, only within a given range of indices
  //
  // fileName = name of the file where the vector has to be read
  // minIndex = index of the first component to read
  // maxIndex = index of the last component to read
  // return value = true if no error occurs
  bool ReadVector (const char* fileName, long minIndex, long maxIndex);

  // read vector dimension from a file, without loading the full vector 
  //
  // fileName = name of the file where the vector has to be read
  // return value = vector dimension
  long ReadVectorDimension (const char* fileName);

  // test if a vector can be read from a file (matching the right type), without loading the full vector 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if the vector can be read
  bool ReadVectorTest (const char* fileName);

  // Output Stream overload
  //
  // str = reference on output stream
  // v = vector to print
  // return value = referenceint GetNumSites(){return this->NSites;} on output stream
  friend ostream& operator << (ostream& str, RationalVector& v);

};
 

// get the numerator of a given component
// 
// index = component index
// return value = reference on the denominator

inline long& RationalVector::Num(int index)
{
  return this->Numerators[index];
}

// get the numerator of a given component
// 
// index = component index
// return value = reference on the denominator

inline long& RationalVector::Num(long index)
{
  return this->Numerators[index];
}

// get the denominator of a given component
// 
// index = component index
// return value = reference on the denominator

inline long& RationalVector::Den(int index)
{
  return this->Denominators[index];
}

// get the denominator of a given component
// 
// index = component index
// return value = reference on the denominator

inline long& RationalVector::Den(long index)
{
  return this->Denominators[index];
}

#endif
