////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of multiple complex scalar product operation            //
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


#ifndef MULTIPLECOMPLEXSCALARPRODUCTOPERATION_H
#define MULTIPLECOMPLEXSCALARPRODUCTOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"


class MultipleComplexScalarProductOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first scalar product to evaluate 
  int FirstScalarProduct;
  // number of component 
  int NbrScalarProduct;

  // array containing the scalar product
  Complex* ScalarProducts;

  // array of vectors to use for the right hand side of the scalar product
  ComplexVector* RightVectors;
  // array of pointers to the vectors to use for the right hand side of the scalar product
  ComplexVector** RightVectorsByPointers;
  // real matrix where vectors to use for the right hand side of the scalar product are stored (can be used instead of RightVectors)
  ComplexMatrix RightVectorMatrix;

  // pointer to the vector to use for the left hand side of the scalar product
  ComplexVector* LeftVector;  

 public:
  
  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = array of vectors to use for the right hand side of the scalar product
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  MultipleComplexScalarProductOperation(ComplexVector* leftVector, ComplexVector* rightVectors, int nbrScalarProduct, Complex* scalarProducts);

  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = array of  pointers to the vectors to use for the right hand side of the scalar product
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  MultipleComplexScalarProductOperation(ComplexVector* leftVector, ComplexVector** rightVectors, int nbrScalarProduct, Complex* scalarProducts);

  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = real matrix where vectors to use for the right hand side of the scalar product are stored
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  MultipleComplexScalarProductOperation(ComplexVector* leftVector, ComplexMatrix& rightVectors, int nbrScalarProduct, Complex* scalarProducts);

  // copy constructor 
  //
  // operation = reference on operation to copy
  MultipleComplexScalarProductOperation(const MultipleComplexScalarProductOperation& operation);
  
  // destructor
  //
  ~MultipleComplexScalarProductOperation();
  
  // set index range of scalar product that have to be calculated
  // 
  // firstScalarProduct = index of the first scalar product to evaluate
  // nbrScalarProduct = number of scalar products that have to be evaluated
  void SetIndicesRange (const int& firstScalarProduct, const int& nbrScalarProduct);

  // get the number of scalar products that have to be evaluated
  // 
  // return value = number of scalar products
  int GetNbrScalarProduct ();

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation
  //
  // return value = true if no error occurs
  bool ApplyOperation();
  
};

// get the number of scalar products that have to be evaluated
// 
// return value = number of scalar products

inline int MultipleComplexScalarProductOperation::GetNbrScalarProduct ()
{
  return this->NbrScalarProduct;
}

#endif
