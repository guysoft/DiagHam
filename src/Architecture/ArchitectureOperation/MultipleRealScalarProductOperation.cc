////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of multiple real scalar product operation             //
//                                                                            //
//                        last modification : 24/10/2002                      //
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
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"


// constructor 
//
// leftVector = pointer to the vector to use for the left hand side of the scalar product
// rightVectors = array of vectors to use for the right hand side of the scalar product
// nbrScalarProduct = number of scalar products that have to be evaluated
// scalarProducts = array where scalar products have to be stored

MultipleRealScalarProductOperation::MultipleRealScalarProductOperation(RealVector* leftVector, RealVector* rightVectors, int nbrScalarProduct, double* scalarProducts)
{
  this->FirstComponent = 0;
  this->NbrComponents = leftVector->GetVectorDimension();
  this->NbrScalarProduct = nbrScalarProduct;
  this->ScalarProducts = scalarProducts;
  this->RightVectors = rightVectors; 
  this->RightVectorsByPointers = 0;
  this->LeftVector = leftVector;
  this->OperationType = AbstractArchitectureOperation::MultipleRealScalarProduct;
}

// constructor 
//
// leftVector = pointer to the vector to use for the left hand side of the scalar product
// rightVectors = array of pointers to the vectors to use for the right hand side of the scalar product
// nbrScalarProduct = number of scalar products that have to be evaluated
// scalarProducts = array where scalar products have to be stored

MultipleRealScalarProductOperation::MultipleRealScalarProductOperation(RealVector* leftVector, RealVector** rightVectors, int nbrScalarProduct, double* scalarProducts)
{
  this->FirstComponent = 0;
  this->NbrComponents = leftVector->GetVectorDimension();
  this->NbrScalarProduct = nbrScalarProduct;
  this->ScalarProducts = scalarProducts;
  this->RightVectors = 0; 
  this->RightVectorsByPointers = rightVectors;
  this->LeftVector = leftVector;
  this->OperationType = AbstractArchitectureOperation::MultipleRealScalarProduct;
}

// constructor 
//
// leftVector = pointer to the vector to use for the left hand side of the scalar product
// rightVectors = real matrix where vectors to use for the right hand side of the scalar product are stored
// nbrScalarProduct = number of scalar products that have to be evaluated
// scalarProducts = array where scalar products have to be stored

MultipleRealScalarProductOperation::MultipleRealScalarProductOperation(RealVector* leftVector, RealMatrix& rightVectors, int nbrScalarProduct, double* scalarProducts)
{
  this->FirstComponent = 0;
  this->NbrComponents = leftVector->GetVectorDimension();
  this->NbrScalarProduct = nbrScalarProduct;
  this->ScalarProducts = scalarProducts;
  this->RightVectors = 0; 
  this->RightVectorsByPointers = 0;
  this->RightVectorMatrix = rightVectors;
  this->LeftVector = leftVector;
  this->OperationType = AbstractArchitectureOperation::MultipleRealScalarProduct;
}

// copy constructor 
//
// operation = reference on operation to copy

MultipleRealScalarProductOperation::MultipleRealScalarProductOperation(const MultipleRealScalarProductOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponents = operation.NbrComponents;
  this->NbrScalarProduct = operation.NbrScalarProduct;
  this->ScalarProducts = operation.ScalarProducts;
  this->RightVectors = operation.RightVectors; 
  this->RightVectorsByPointers = operation.RightVectorsByPointers; 
  this->RightVectorMatrix = operation.RightVectorMatrix;
  this->LeftVector = operation.LeftVector;
  this->OperationType = AbstractArchitectureOperation::MultipleRealScalarProduct;
}
  
// destructor
//

MultipleRealScalarProductOperation::~MultipleRealScalarProductOperation()
{
}
  
// set the array where scalar products have to be stored
//
// scalarProducts = array where scalar products have to be stored

void MultipleRealScalarProductOperation::SetScalarProducts (double* scalarProducts)
{
  this->ScalarProducts = scalarProducts;
}

// set index range of scalar product that have to be calculated
// 
// firstComponent = index of the first component of each partial scalar product
// nbrComponent = number of component to take into account for each partial scalar product

void MultipleRealScalarProductOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponents = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* MultipleRealScalarProductOperation::Clone()
{
  return new MultipleRealScalarProductOperation (*this);
}
  
// apply operation
//
// return value = true if no error occurs

bool MultipleRealScalarProductOperation::ApplyOperation()
{
  if (this->RightVectors != 0)
    {
      for (int i = 0; i < this->NbrScalarProduct; ++i)
	{
	  this->ScalarProducts[i] = this->LeftVector->PartialScalarProduct(this->RightVectors[i], this->FirstComponent, this->NbrComponents);
	}
    }
  else
    if (this->RightVectorsByPointers != 0)
      {
	for (int i = 0; i < this->NbrScalarProduct; ++i)
	  {
	    this->ScalarProducts[i] = this->LeftVector->PartialScalarProduct((*(this->RightVectorsByPointers[i])), this->FirstComponent, this->NbrComponents);
	  }
      }
    else
      {
	for (int i = 0; i < this->NbrScalarProduct; ++i)
	  {
	    this->ScalarProducts[i] = this->LeftVector->PartialScalarProduct(this->RightVectorMatrix[i], this->FirstComponent, this->NbrComponents);
	  }
      }
  return true;
}
