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


#ifndef OPERATORMATRIXELEMENTOPERATION_H
#define OPERATORMATRIXELEMENTOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractScalarSumOperation.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"


class OperatorMatrixElementOperation: public AbstractScalarSumOperation
{

 protected:

  // index of the first component of each partial scalar product
  int FirstComponent;
  // number of component to take into account for each partial scalar product
  int NbrComponents;

  // number of scalar products
  int NbrScalarProduct;

  // strategy used to do the scalar products
  int Strategy;

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

  enum StrategyType
    {
      // each process achieves part of each scalar product
      VectorSubdivision = 0x01,
      // each process achieves full scalar product for a given group of vectors
      GroupSubdivision = 0x02
    };
  
  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = array of vectors to use for the right hand side of the scalar product
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  OperatorMatrixElementOperation(RealVector* leftVector, RealVector* rightVectors, int nbrScalarProduct, double* scalarProducts);

  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = array of pointers to the vectors to use for the right hand side of the scalar product
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  OperatorMatrixElementOperation(RealVector* leftVector, RealVector** rightVectors, int nbrScalarProduct, double* scalarProducts);

  // constructor 
  //
  // leftVector = pointer to the vector to use for the left hand side of the scalar product
  // rightVectors = real matrix where vectors to use for the right hand side of the scalar product are stored
  // nbrScalarProduct = number of scalar products that have to be evaluated
  // scalarProducts = array where scalar products have to be stored
  OperatorMatrixElementOperation(RealVector* leftVector, RealMatrix& rightVectors, int nbrScalarProduct, double* scalarProducts);

  // copy constructor 
  //
  // operation = reference on operation to copy
  OperatorMatrixElementOperation(const OperatorMatrixElementOperation& operation);
  
  // destructor
  //
  ~OperatorMatrixElementOperation();
  
  // apply operation(architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
};

#endif
