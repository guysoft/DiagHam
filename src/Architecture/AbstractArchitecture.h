////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of Abstract Architecture                      //
//                                                                            //
//                        last modification : 10/04/2002                      //
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


#ifndef ABSTRACTARCHITECTURE_H
#define ABSTRACTARCHITECTURE_H


#include "config.h"


class Vector;
class RealVector;
class ComplexVector;

class AbstractHamiltonian;

class AbstractArchitectureOperation;
class VectorHamiltonianMultiplyOperation;
class AddRealLinearCombinationOperation;
class AddComplexLinearCombinationOperation;
class MultipleRealScalarProductOperation;
class MultipleComplexScalarProductOperation;
class MatrixMatrixMultiplyOperation;
class AbstractPrecalculationOperation;


class AbstractArchitecture
{

public:
  
  // destructor
  //
  virtual ~AbstractArchitecture();
  
  // get typical range of indices on which the local architecture acts
  //
  // minIndex = reference on the minimum index on which the local architecture can act
  // maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
  //            architecture doesn't support this feature)
  virtual void GetTypicalRange (long& minIndex, long& maxIndex);
  
  // get a new real vector with memory alloaction depending on the architecture
  //
  // return value = pointer to the requested vector (zero if an error occurs)
  virtual RealVector* GetNewRealVector ();
  
  // get a new real vector with memory alloaction depending on the architecture
  //
  // dimension = dimension of the requested vector
  // zeroFlag = true if all vector entries has to be set to zero
  // return value = pointer to the requested vector (zero if an error occurs)
  virtual RealVector* GetNewRealVector (long dimension, bool zeroFlag = false);
  
  // get a new complex vector with memory alloaction depending on the architecture
  //
  // return value = pointer to the requested vector (zero if an error occurs)
  virtual ComplexVector* GetNewComplexVector ();
  
  // get a new complex vector with memory alloaction depending on the architecture
  //
  // dimension = dimension of the requested vector
  // zeroFlag = true if all vector entries has to be set to zero
  // return value = pointer to the requested vector (zero if an error occurs)
  virtual ComplexVector* GetNewComplexVector (long dimension, bool zeroFlag = false);
  
  // multiply a vector by an hamiltonian and store the result in another vector
  //
  // hamiltonian = pointer to the hamiltonian to use
  // vSource = vector to multiply 
  // vDestination = vector where result has to be stored 
  virtual void Multiply (AbstractHamiltonian* hamiltonian, Vector& vSource, Vector& vDestination) = 0;

  // execute an architecture-dependent operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (AbstractArchitectureOperation* operation);

  // execute an architecture-dependent vector hamiltonian multiplication operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (VectorHamiltonianMultiplyOperation* operation);
  
  // execute an architecture-dependent add real linear combination operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (AddRealLinearCombinationOperation* operation);
  
  // execute an architecture-dependent add complex linear combination operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (AddComplexLinearCombinationOperation* operation);

  // execute an architecture-dependent multiple real scalar product operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (MultipleRealScalarProductOperation* operation);
  
  // execute an architecture-dependent multiple complex scalar product operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (MultipleComplexScalarProductOperation* operation);
  
  // execute an architecture-dependent matrix matrix multiplication operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (MatrixMatrixMultiplyOperation* operation);
    
  // execute an architecture-dependent abstract hamiltonian precalculation operation
  //
  // operation = pointer to the operation to execute
  // return value = true if operation has been completed successfully
  virtual bool ExecuteOperation (AbstractPrecalculationOperation* operation);

  // request a given amount of memory for an array of Type element
  //
  // architecture = reference on the architecture to which memory will be asked
  // pointer = reference on pointer that will be used for the memory allocation
  // size = number of elements of type Type
  // return value = reference on the pointer
  template <class Type>
  friend Type*& New (AbstractArchitecture& architecture, Type*& pointer,  unsigned long size);
  
  // delete an array that has been requested by the New function
  //
  // architecture = reference on the architecture to which memory has been asked
  // pointer = reference on pointer that of the memory allocation
  template <class Type>
  friend void Delete (AbstractArchitecture& architecture, Type*& pointer);
  
};


// request a given amount of memory for an array of Type element
//
// pointer = reference on pointer that will be used for the memory allocation
// size = number of elements of type Type
// return value = reference on the pointer

template <class Type>
inline Type*& New (AbstractArchitecture& architecture, Type*& pointer, unsigned long size)
{
  pointer = new Type [size];
  return pointer;
}
  
// delete an array that has been requested by the New function
//
// architecture = reference on the architecture to which memory has been asked
// pointer = reference on pointer that of the memory allocation

template <class Type>
inline void Delete (AbstractArchitecture& architecture, Type*& pointer)
{
  delete[] pointer;
}
  

#endif
