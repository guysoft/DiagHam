////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of add linear combination operation                //
//                                                                            //
//                        last modification : 26/05/2003                      //
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
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Vector/ComplexVector.h"
#include "Complex.h"


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// sourceVector = vector to be multiplied by the hamiltonian
// destinationVector = vector where the result has to be stored

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation (ComplexVector* destinationVector, ComplexVector* sourceVector, 
									    int nbrVector, Complex* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = sourceVector;
  this->NbrVector = nbrVector;
  this->Coefficients = coefficients;
  this->RealCoefficients = 0;
  this->DestinationVector = destinationVector;  
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = matrix containing the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation(ComplexVector* destinationVector, ComplexMatrix& sourceVector, 
									   int nbrVector, Complex* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = 0;
  this->SourceVectorMatrix = sourceVector;
  this->NbrVector = nbrVector;
  this->Coefficients = coefficients;
  this->RealCoefficients = 0;
  this->DestinationVector = destinationVector;  
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;  
}

// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// sourceVector = vector to be multiplied by the hamiltonian
// destinationVector = vector where the result has to be stored

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation (ComplexVector* destinationVector, ComplexVector* sourceVector, 
									    int nbrVector, double* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = sourceVector;
  this->NbrVector = nbrVector;
  this->Coefficients = 0;
  this->RealCoefficients = coefficients;
  this->DestinationVector = destinationVector;  
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = matrix containing the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation(ComplexVector* destinationVector, ComplexMatrix& sourceVector, 
									   int nbrVector, double* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = 0;
  this->SourceVectorMatrix = sourceVector;
  this->NbrVector = nbrVector;
  this->Coefficients = 0;
  this->RealCoefficients = coefficients;
  this->DestinationVector = destinationVector;  
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;  
}

// copy constructor 
//
// operation = reference on operation to copy

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation(const AddComplexLinearCombinationOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->NbrVector = operation.NbrVector;
  this->Coefficients = operation.Coefficients;
  this->RealCoefficients = operation.RealCoefficients;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
  this->SourceVector = operation.SourceVector;
  this->DestinationVector = operation.DestinationVector;  
  this->SourceVectorMatrix = operation.SourceVectorMatrix;
}
  
// destructor
//

AddComplexLinearCombinationOperation::~AddComplexLinearCombinationOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void AddComplexLinearCombinationOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* AddComplexLinearCombinationOperation::Clone()
{
  return new AddComplexLinearCombinationOperation (*this);
}
  
// apply operation
//
// return value = true if no error occurs

bool AddComplexLinearCombinationOperation::ApplyOperation()
{
  if (this->RealCoefficients == 0)
    {
      if (this->SourceVector != 0)
	for (int i = 0; i < this->NbrVector; ++i)
	  {
	    this->DestinationVector->AddLinearCombination(Coefficients[i], (this->SourceVector[i]), this->FirstComponent, 
							  this->NbrComponent);
	  }
      else
	for (int i = 0; i < this->NbrVector; ++i)
	  {
	    this->DestinationVector->AddLinearCombination(Coefficients[i], (this->SourceVectorMatrix[i]), this->FirstComponent, 
							  this->NbrComponent);
	  }
    }
  else
    {
      if (this->SourceVector != 0)
	for (int i = 0; i < this->NbrVector; ++i)
	  {
	    this->DestinationVector->AddLinearCombination(RealCoefficients[i], (this->SourceVector[i]), this->FirstComponent, 
							  this->NbrComponent);
	  }
      else
	for (int i = 0; i < this->NbrVector; ++i)
	  {
	    this->DestinationVector->AddLinearCombination(RealCoefficients[i], (this->SourceVectorMatrix[i]), this->FirstComponent, 
							  this->NbrComponent);
	  }
    }
  return true;
}

