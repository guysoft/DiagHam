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


#ifndef MULTIPLEREALSCALARPRODUCTOPERATION_H
#define MULTIPLEREALSCALARPRODUCTOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"


class MultipleRealScalarProductOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component of each partial scalar product
  int FirstComponent;
  // number of component to take into account for each partial scalar product
  int NbrComponents;

  // number of scalar products
  int NbrScalarProduct;

  // array containing the scalar products
  double* ScalarProducts;

  // array of vectors to use for the right hand side of the scalar product
  RealVector* RightVectors;
  // array of pointers to the vectors to use for the right hand side of the scalar product
  RealVector** RightVectorsByPointers;
  // real matrix where vectors to use for the right hand side of the scalar product are stored (can be used instead of RightVectors)
  RealMatrix RightVectorMatrix;

  // pointer to the vector to use for the left hand side of the scalar product
  RealVector* LeftVector;  

 public:
  
  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = array of vectors to use for the right hand side of the scalar product
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  MultipleRealScalarProductOperation(RealVector* leftVector, RealVector* rightVectors, int nbrScalarProduct, double* scalarProducts);

  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = array of pointers to the vectors to use for the right hand side of the scalar product
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  MultipleRealScalarProductOperation(RealVector* leftVector, RealVector** rightVectors, int nbrScalarProduct, double* scalarProducts);

  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = real matrix where vectors to use for the right hand side of the scalar product are stored
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  MultipleRealScalarProductOperation(RealVector* leftVector, RealMatrix& rightVectors, int nbrScalarProduct, double* scalarProducts);

  // copy constructor 
  //
  // operation = reference on operation to copy
  MultipleRealScalarProductOperation(const MultipleRealScalarProductOperation& operation);
  
  // destructor
  //
  ~MultipleRealScalarProductOperation();
  
  // set the array where scalar products have to be stored
  //
  // scalarProducts = array where scalar products have to be stored
  void SetScalarProducts (double* scalarProducts);

  // get the array where scalar products have to be stored
  //
  // return value = array where scalar products have to be stored
  double* GetScalarProducts ();

  // get the vector used as the left hand side of the scalar product
  //
  // return value = pointer to the vector
  RealVector* GetLeftVector();

  // set index range of scalar product that have to be calculated
  // 
  // firstComponent = index of the first component of each partial scalar product
  // nbrComponent = number of component to take into account for each partial scalar product
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

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

inline int MultipleRealScalarProductOperation::GetNbrScalarProduct ()
{
  return this->NbrScalarProduct;
}

// get the array where scalar products have to be stored
//
// return value = array where scalar products have to be stored

inline double* MultipleRealScalarProductOperation::GetScalarProducts ()
{
  return this->ScalarProducts;
}

// get the vector used as the left hand side of the scalar product
//
// return value = pointer to the vector

inline RealVector* MultipleRealScalarProductOperation::GetLeftVector()
{
  return this->LeftVector;
}

#endif
