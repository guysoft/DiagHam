////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of operator martrix element evaluation operation         //
//                                                                            //
//                        last modification : 12/02/2010                      //
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
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"
#include "Operator/AbstractOperator.h"


// constructor 
//
// oper = operator which has to be evaluated  
// leftVector = pointer to the vector to use for the left hand side of the matrix element
// rightVector = pointer to the vector to use for the right hand side of the matrix element

OperatorMatrixElementOperation::OperatorMatrixElementOperation(AbstractOperator* oper, RealVector& leftVector, RealVector& rightVector)
{
  this->RealLeftVector = leftVector;
  this->RealRightVector = rightVector;
  this->Operator = oper->Clone();
  this->OperationType = AbstractArchitectureOperation::OperatorMatrixElement;  
  this->FirstComponent = 0;
  this->NbrComponent = rightVector.GetVectorDimension();
  this->LargeFirstComponent = 0l;
  this->LargeNbrComponent = rightVector.GetLargeVectorDimension();
  this->Scalars = 0;
  this->NbrScalars = 1;
}

// constructor 
//
// oper = operator which has to be evaluated  
// leftVector = pointer to the vector to use for the left hand side of the matrix element
// rightVector = pointer to the vector to use for the right hand side of the matrix element

OperatorMatrixElementOperation::OperatorMatrixElementOperation(AbstractOperator* oper, ComplexVector& leftVector, ComplexVector& rightVector)
{
  this->ComplexLeftVector = leftVector;
  this->ComplexRightVector = rightVector;
  this->Operator = oper->Clone();
  this->OperationType = AbstractArchitectureOperation::OperatorMatrixElement;
  this->FirstComponent = 0;
  this->NbrComponent = leftVector.GetVectorDimension();
  this->LargeFirstComponent = 0l;
  this->LargeNbrComponent = leftVector.GetLargeVectorDimension();
  this->Scalars = 0;
  this->NbrScalars = 1;
}
 
// copy constructor 
//
// operation = reference on operation to copy

OperatorMatrixElementOperation::OperatorMatrixElementOperation(const OperatorMatrixElementOperation& operation)
{
  this->RealRightVector = operation.RealRightVector;
  this->RealLeftVector = operation.RealLeftVector;
  this->ComplexRightVector = operation.ComplexRightVector;
  this->ComplexLeftVector = operation.ComplexLeftVector;
  this->Operator = operation.Operator->Clone();
  this->OperationType = operation.OperationType;
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LargeFirstComponent = operation.LargeFirstComponent;
  this->LargeNbrComponent = operation.LargeNbrComponent;
  this->Scalars = 0;
  this->NbrScalars = 1;
}
  
// destructor
//

OperatorMatrixElementOperation::~OperatorMatrixElementOperation()
{
  delete this->Operator;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* OperatorMatrixElementOperation::Clone()
{
  return new OperatorMatrixElementOperation (*this);
}
  
// apply operation(architecture independent)
//
// return value = true if no error occurs

bool OperatorMatrixElementOperation::RawApplyOperation()
{
  if (this->RealRightVector.GetVectorDimension() > 0)
    {
      this->Scalar = this->Operator->PartialMatrixElement(this->RealLeftVector, this->RealRightVector, this->LargeFirstComponent, this->LargeNbrComponent);
    }
  else
    {
      this->Scalar = this->Operator->PartialMatrixElement(this->ComplexLeftVector, this->ComplexRightVector, this->LargeFirstComponent, this->LargeNbrComponent);
    }
  return true;
}
