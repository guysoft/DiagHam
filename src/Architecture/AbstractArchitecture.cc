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
//                        last modification : 30/04/2002                      //
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
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/MultipleVectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MatrixMatrixMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "Architecture/ArchitectureOperation/AbstractScalarSumOperation.h"
#include "Hamiltonian/AbstractHamiltonian.h"
//#include "Architecture/ArchitectureOperation/GenericOperation.h"

#include "MainTask/AbstractMainTask.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include <sys/time.h>
#include <string.h>


// destructor
//

AbstractArchitecture::~AbstractArchitecture()
{
}

// get typical range of indices on which the local architecture acts
//
// minIndex = reference on the minimum index on which the local architecture can act
// maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
//            architecture doesn't support this feature)

void AbstractArchitecture::GetTypicalRange (long& minIndex, long& maxIndex)
{
  minIndex = 0;
  maxIndex = this->HilbertSpaceDimension - 1;
}
  
// get a new real vector with memory alloaction depending on the architecture
//
// return value = pointer to the requested vector (zero if an error occurs)

RealVector* AbstractArchitecture::GetNewRealVector ()
{
  return new RealVector;
}
  
// get a new real vector with memory alloaction depending on the architecture
//
// dimension = dimension of the requested vector
// zeroFlag = true if all vector entries has to be set to zero
// return value = pointer to the requested vector (zero if an error occurs)

RealVector* AbstractArchitecture::GetNewRealVector (long dimension, bool zeroFlag)
{
  return new RealVector(dimension, zeroFlag);
}
  
// get a new complex vector with memory alloaction depending on the architecture
//
// return value = pointer to the requested vector (zero if an error occurs)

ComplexVector* AbstractArchitecture::GetNewComplexVector ()
{
  return new ComplexVector;
}
  
// get a new complex vector with memory alloaction depending on the architecture
//
// dimension = dimension of the requested vector
// zeroFlag = true if all vector entries has to be set to zero
// return value = pointer to the requested vector (zero if an error occurs)

ComplexVector* AbstractArchitecture::GetNewComplexVector (long dimension, bool zeroFlag)
{
  return new ComplexVector (dimension, zeroFlag);
}
  
// set dimension of the Hilbert space on which the architecture has to work
// 
// dimension = dimension of the Hilbert space

void AbstractArchitecture::SetDimension (long dimension)
{
  this->HilbertSpaceDimension = dimension;
}

// multiply a vector by an hamiltonian and store the result in another vector
//
// hamiltonian = pointer to the hamiltonian to use
// vSource = vector to multiply 
// vDestination = vector where result has to be stored 

void AbstractArchitecture::Multiply (AbstractHamiltonian* hamiltonian, Vector& vSource, Vector& vDestination)
{
  hamiltonian->Multiply(vSource, vDestination);
}

// execute an architecture-dependent operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (AbstractArchitectureOperation* operation)
{
  switch (operation->GetOperationType())
    {
    case AbstractArchitectureOperation::VectorHamiltonianMultiply:
      return this->ExecuteOperation((VectorHamiltonianMultiplyOperation*) operation);
      break;
    case AbstractArchitectureOperation::MultipleVectorHamiltonianMultiply:
      return this->ExecuteOperation((MultipleVectorHamiltonianMultiplyOperation*) operation);
      break;
    case AbstractArchitectureOperation::AddRealLinearCombination:
      return this->ExecuteOperation((AddRealLinearCombinationOperation*) operation);
      break;
    case AbstractArchitectureOperation::AddComplexLinearCombination:
      return this->ExecuteOperation((AddComplexLinearCombinationOperation*) operation);
      break;
    case AbstractArchitectureOperation::MultipleRealScalarProduct:
      return this->ExecuteOperation((MultipleRealScalarProductOperation*) operation);
      break;
    case AbstractArchitectureOperation::MultipleComplexScalarProduct:
      return this->ExecuteOperation((MultipleComplexScalarProductOperation*) operation);
      break;
    case AbstractArchitectureOperation::MatrixMatrixMultiply:
      return this->ExecuteOperation((MatrixMatrixMultiplyOperation*) operation);
      break;
    case AbstractArchitectureOperation::MainTask:
      return this->ExecuteOperation((MainTaskOperation*) operation);
      break;
    case AbstractArchitectureOperation::ScalarSum:
      return this->ExecuteOperation((AbstractScalarSumOperation*) operation);
      break;
//    case AbstractArchitectureOperation::Generic:
//      return this->ExecuteOperation((GenericOperation*) operation);
//      break;
    default:
      return false;
    }
  return false;
}

// execute an architecture-dependent vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (VectorHamiltonianMultiplyOperation* operation)
{
  return false;
}

// execute an architecture-dependent multiple vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (MultipleVectorHamiltonianMultiplyOperation* operation)
{
  return false;
}

// execute an architecture-dependent vector abstact scalar sum operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (AbstractScalarSumOperation* operation)
{
  return false;
}
                                                                                                                                                                                  
// execute an architecture-dependent add real linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (AddRealLinearCombinationOperation* operation)
{
  return false;
}

// execute an architecture-dependent add complex linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (AddComplexLinearCombinationOperation* operation)
{
  return false;
}

// execute an architecture-dependent multiple real scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (MultipleRealScalarProductOperation* operation)
{
  return false;
}

// execute an architecture-dependent multiple complex scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (MultipleComplexScalarProductOperation* operation)
{
  return false;
}

// execute an architecture-dependent matrix matrix multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (MatrixMatrixMultiplyOperation* operation)
{
  return false;
}

// execute an architecture-dependent abstract hamiltonian precalculation operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (AbstractPrecalculationOperation* operation)
{
  return false;
}
    
// execute an architecture-dependent main task operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool AbstractArchitecture::ExecuteOperation (MainTaskOperation* operation)
{
  operation->GetMainTask()->SetArchitecture(this);
  return operation->ApplyOperation();
}
    
// get a temporary file name
//
// return value = string corresponding to a temporary file name

char* AbstractArchitecture::GetTemporaryFileName()
{
  timeval Time;
  gettimeofday (&Time, 0);
  char* TmpString = new char [32];
  sprintf (TmpString, "diagam%d%d.tmp",(int)  Time.tv_sec, (int)  Time.tv_usec);
  return TmpString;
}
  
