////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class for n dimensional real vector                     //
//            whose memory allocation is done only on one process             //
//                                                                            //
//                        last modification : 15/06/2004                      //
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


#ifndef DELOCALIZEDREALVECTOR_H
#define DELOCALIZEDREALVECTOR_H


#include "config.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"
#include "Architecture/ClusterArchitecture/AbstractClusterArchitecture.h"

#include <iostream>
#include <fstream>


using std::ostream;
using std::ifstream;
using std::ofstream;


class Complex;
class ComplexVector;
class MathematicaOutput;
class BlockDiagonalMatrix;
class Matrix;
class RealMatrix;
class RealAntisymmetricMatrix;
class RealDiagonalMatrix;
class RealSymmetricMatrix;
class RealTriDiagonalSymmetricMatrix;
class ComplexVector;


class DelocalizedRealVector : public RealVector
{

 protected:
  
  // id of the node where the vector is localized
  int LocalizationId;
  // id of the local node
  int LocalId;
  // pointer to the cluster architecture in use
  AbstractClusterArchitecture* Architecture;

 public:

  // default constructor
  //
  DelocalizedRealVector();

  // constructor for an empty real vector 
  //
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  DelocalizedRealVector(int size, bool zeroFlag = false);

  // constructor from an array of doubles
  //
  // array = array of doubles with real in even position and imaginary part in odd position
  // size = Vector Dimension  
  DelocalizedRealVector(double* array, int size);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  DelocalizedRealVector(const DelocalizedRealVector& vector, bool duplicateFlag = false);

  // copy constructor from a real vector
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  DelocalizedRealVector(const RealVector& vector, bool duplicateFlag = false);

  // copy constructor from a complex vector (keep only real part and datas are duplicated)
  //
  // vector = vector to copy
  DelocalizedRealVector(const ComplexVector& vector);

  // copy constructor from a vector (duplicate datas if necessary)
  //
  // vector = vector to copy
  DelocalizedRealVector(const Vector& vector);

  // destructor
  //
  ~DelocalizedRealVector ();

  // assignement
  //
  // vector = vector to assign
  // return value = reference on current vector
  DelocalizedRealVector& operator = (const DelocalizedRealVector& vector);

  // assignement from a real vector
  //
  // vector = vector to assign
  // return value = reference on current vector
  DelocalizedRealVector& operator = (const RealVector& vector);

  // assignement from a complex vector (keep only real part and datas are duplicated)
  //
  // vector = vector to assign
  // return value = reference on current vector
  DelocalizedRealVector& operator = (const ComplexVector& vector);

  // assignement from a vector (duplicate datas if necessary)
  //
  // vector = vector to assign
  // return value = reference on current vector
  DelocalizedRealVector& operator = (const Vector& vector);

  // Resize vector
  //
  // dimension = new dimension
  void Resize (int dimension);

  // Resize vector and set to zero all components that have been added
  //
  // dimension = new dimension
  void ResizeAndClean (int dimension);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  DelocalizedRealVector& Copy (DelocalizedRealVector& vector, double coefficient = 1.0);

  // create a new vector with same size and same type but without duplicating datas
  //
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  Vector* EmptyClone(bool zeroFlag = false);

  // put all vector components to zero
  //
  // return value = reference on current vector
  Vector& ClearVector ();

  // change sign of a vector
  //
  // return value = reference on current vector
  DelocalizedRealVector& operator - ();

  // return a new vector with opposite sign form a given source vector
  //
  // V1 = source vector
  // return value = new vector
  friend DelocalizedRealVector operator - (const DelocalizedRealVector& V1);

  // scalar product between two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = result of scalar product
  friend double operator * (const DelocalizedRealVector& V1, const DelocalizedRealVector& V2);

  // do part of the scalar product between two vectors in a given range of indices
  //
  // vRight = right vector of the scalar product
  // firstComponent = index of the first component to consider
  // nbrComponent = number of components to consider
  // step = increment between to consecutive indices to consider
  // return value = result of the partial scalar product
  double PartialScalarProduct (const DelocalizedRealVector& vRight, int firstComponent, int nbrComponent, int step = 1);

  // sum two vectors
  //
  // V1 = vector to add
  // return value = reference on current vector
  DelocalizedRealVector& operator += (const DelocalizedRealVector& V1);

  // sum two vectors
  //
  // vector = vector to add
  // return value = reference on current vector
  Vector& operator += (const Vector& vector);

  // substract two vectors
  //
  // V1 = first vector
  // return value = reference on current vector
  DelocalizedRealVector& operator -= (const DelocalizedRealVector& V1);

  // sum two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = resulting vector
  friend DelocalizedRealVector operator + (const DelocalizedRealVector& V1, const DelocalizedRealVector& V2);

  // substract two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = resulting vector
  friend DelocalizedRealVector operator - (const DelocalizedRealVector& V1, const DelocalizedRealVector& V2);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  DelocalizedRealVector& AddLinearCombination (const double& x, const DelocalizedRealVector& V);

  // add a linear combination to a given vector, for a given range of indices
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  DelocalizedRealVector& AddLinearCombination (double x, const DelocalizedRealVector& V, int firstComponent, int nbrComponent);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  DelocalizedRealVector& AddLinearCombination (double x1, const DelocalizedRealVector& v1, double x2, const DelocalizedRealVector& v2);

  // add a linear combination of two vectors to a given vector, for a given range of indices
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  DelocalizedRealVector& AddLinearCombination (double x1, const DelocalizedRealVector& v1, double x2, 
				    const DelocalizedRealVector& v2, int firstComponent, int nbrComponent);

  // multiply a vector with a real number on the right hand side
  //
  // V1 = vector to multiply
  // d = real to use
  // return value = resulting vector
  friend DelocalizedRealVector operator * (const DelocalizedRealVector& V1, double d);

  // multiply a vector with a real number on the left hand side
  //
  // V1 = vector to multiply
  // d = real to use
  // return value = resulting vector
  friend DelocalizedRealVector operator * (double d, const DelocalizedRealVector& V1);

  // multiply a vector with a real number on the right hand side
  //
  // d = real to use
  // return value = reference on current vector
  DelocalizedRealVector& operator *= (double d);

  // divide a vector by a real number on the right hand side
  //
  // d = real to use
  // return value = reference on current vector
  DelocalizedRealVector& operator /= (double d);

  // left multiply a vector with a real symmetric matrix (without using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  DelocalizedRealVector& operator *= (const RealSymmetricMatrix&  M);

  // left multiply a vector with a real tridiagonal symmetric matrix (without using temporary vector)
  //
  // M = matrix to use
  // return value = reference on current vector
  DelocalizedRealVector& operator *= (const RealTriDiagonalSymmetricMatrix&  M);

  // left multiply a vector with a symmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply  
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealSymmetricMatrix&  M, DelocalizedRealVector& V);

  // do a partial left multication of a vector with a real symmetric matrix and store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealSymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a symmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealSymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, int destStart, int destStep);

  // left multiply a vector with a symmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealSymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealSymmetricMatrix&  M, DelocalizedRealVector& V);

  // do a partial left multication of a vector with a real symmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealSymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealSymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, 
			   int sourceStep, int destStart, int destStep);

  // left multiply a vector with a antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealSymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, 
			   int sourceNbrComponent, int sourceStep, int destStart, int destStep);

  // left multiply a vector with an antisymmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply  
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealAntisymmetricMatrix&  M, DelocalizedRealVector& V);

  // do a partial left multication of a vector with a real antisymmetric matrix and store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealAntisymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with an antisymmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealAntisymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			int destStart, int destStep);

  // left multiply a vector with an antisymmetric matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealAntisymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with an antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealAntisymmetricMatrix&  M, DelocalizedRealVector& V);

  // do a partial left multication of a vector with a real antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealAntisymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with an antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealAntisymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, int destStart, int destStep);

  // left multiply a vector with an antisymmetric matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealAntisymmetricMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			   int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a real diagonal matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply  
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealDiagonalMatrix&  M, DelocalizedRealVector& V);

  // do a partial left multication of a vector with a real matrix and store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealDiagonalMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a real diagonal matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealDiagonalMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			int destStart, int destStep);

  // left multiply a vector with a real diagonal matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealDiagonalMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a real diagonal matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealDiagonalMatrix&  M, DelocalizedRealVector& V);
  
  // do a partial left multication of a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealDiagonalMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceNbrComponent);
    
  // left multiply a vector with a real diagonal matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealDiagonalMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, int destStart, int destStep);

  // left multiply a vector with a real diagonal matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply  
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealDiagonalMatrix&  M, DelocalizedRealVector& V, int sourceStart, 
			   int sourceNbrComponent, int sourceStep, int destStart, int destStep);

  // left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealMatrix&  M, DelocalizedRealVector& V);

  // do a partial left multication of a vector with a real matrix and store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const RealMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector  
  DelocalizedRealVector& Multiply (const RealMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, int destStart, int destStep);

  // left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector  
  DelocalizedRealVector& Multiply (const RealMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			int sourceNbrComponent, int destStart, int destStep);

  // do a partial left multication of a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceNbrComponent = number of component to take into account in the source vector
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, int destStart, int destStep);

  // left multiply a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const RealMatrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			   int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a block-diagonal matrix and use to store result 
  // in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const BlockDiagonalMatrix&  M, DelocalizedRealVector& V);

  // left multiply a vector with a block-diagonal matrix and add result to the current vector
  // in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const BlockDiagonalMatrix&  M, DelocalizedRealVector& V);

  // left multiply a vector with a block-diagonal matrix and use to store result 
  // in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector  
  DelocalizedRealVector& Multiply (const BlockDiagonalMatrix&  M, DelocalizedRealVector& V, int sourceStart, 
			int sourceStep, int destStart, int destStep);

  // left multiply a vector with a block-diagonal matrix and use to store result 
  // in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector  
  DelocalizedRealVector& Multiply (const BlockDiagonalMatrix&  M, DelocalizedRealVector& V, int sourceStart, 
			int sourceStep, int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a block-diagonal matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const BlockDiagonalMatrix&  M, DelocalizedRealVector& V, int sourceStart, 
			   int sourceStep, int destStart, int destStep);

  // left multiply a vector with a block-diagonal matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const BlockDiagonalMatrix&  M, DelocalizedRealVector& V, int sourceStart, 
			   int sourceStep, int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a matrix and use to store result in current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const Matrix&  M, DelocalizedRealVector& V);

  // left multiply a vector with a matrix and add result to current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const Matrix&  M, DelocalizedRealVector& V);

  // do a partial left multication of a vector with a real matrix and store result in current vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const Matrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceNbrComponent);

  // do a partial left multication of a vector with a real matrix and add result to the current vector
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const Matrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceNbrComponent);

  // left multiply a vector with a matrix and use to store result in current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const Matrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			int destStart, int destStep);

  // left multiply a vector with a matrix and add current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const Matrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			   int destStart, int destStep);

  // left multiply a vector with a matrix and use to store result in current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& Multiply (const Matrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			int sourceNbrComponent, int destStart, int destStep);

  // left multiply a vector with a matrix and add current 
  // vector (without creating temporary vector)
  //
  // M = matrix to use
  // V = vector to multiply
  // sourceStart = source vector first coordinate to modify
  // sourceStep = step to add to go to the following source vector coordinate
  // sourceNbrComponent = number of component to take into account in the source vector
  // destStart = destination vector first coordinate to modify
  // destStep = step to add to go to the following destination vector coordinate
  // return value = reference on current vector
  DelocalizedRealVector& AddMultiply (const Matrix&  M, DelocalizedRealVector& V, int sourceStart, int sourceStep, 
			   int sourceNbrComponent, int destStart, int destStep);

  // return vector i-th coordinate (without testing if position is valid)
  //
  // i = coordinate position
  double& operator [] (int i);

  // get vector norm
  //
  // return value = vector norm
  double Norm();
  
  // get square of vector norm
  //
  // return value = square of vector norm
  double SqrNorm ();
  
  // normalize vector
  //
  // return value = reference on current vector
  DelocalizedRealVector& Normalize();

  // orthonormalized a vector with respect to a set of orthonormalized vectors
  //
  // vectors = vector array corresponding to the set
  // nbrVectors = number of vectors in the set
  // return value = resulting vector norm (can be used to see if vector is can be decomposed on vector set)
  double Orthonormalized (DelocalizedRealVector* vectors, int nbrVectors);

  // Extract a subvector from a given vector
  //
  // firstCoordinate = Coordinate where extraction has to begin
  // lastCoordinate = Coordinate where extraction has to stop (extract also include this last coordinate)
  // step = distance to the next coordinate (1 means to take the folowing)
  // return value = return corresponding subvector
  DelocalizedRealVector Extract(int firstCoordinate, int lastCoordinate, int step = 1);
  
  // Merge a subvector into a given vector
  //
  // V = vector to merge
  // firstCoordinate = Coordinate where merge has to begin
  // step = distance to the next coordinate in the destination vector (1 means to take the following)
  // return value = reference to the current Vector
  DelocalizedRealVector& Merge(const DelocalizedRealVector& V, int firstCoordinate, int step = 1);
  
  // write vector in a file 
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteVector (char* fileName);

  // write vector in a file in ascii mode
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteAsciiVector (char* fileName);

  // read vector from a file 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if no error occurs
  bool ReadVector (char* fileName);

  // Output Stream overload
  //
  // str = reference on output stream
  // v = vector to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const DelocalizedRealVector& v);

  // output file stream overload
  //
  // file = reference on output file stream
  // vector = reference on vector to save
  // return value = reference on output file stream
//  friend ofstream& operator << (ofstream& file, const DelocalizedRealVector& vector);

  // input file stream overload
  //
  // file = reference on output file stream
  // vector = reference on vector to load
  // return value = reference on output file stream
  friend ifstream& operator >> (ifstream& file, DelocalizedRealVector& vector);

#ifdef __MPI__

  // send a vector to a given MPI process
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // return value = reference on the current vector
  Vector& SendVector(MPI::Intracomm& communicator, int id);

  // broadcast a vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // return value = reference on the current vector
  Vector& BroadcastVector(MPI::Intracomm& communicator,  int id);

  // broadcast part of vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // firstComponent = index of the first component (useless if the method is not called by the MPI process which broadcasts the vector)
  // nbrComponent = number of component (useless if the method is not called by the MPI process which broadcasts the vector)
  // return value = reference on the current vector
  Vector& BroadcastPartialVector(MPI::Intracomm& communicator, int id, int firstComponent = 0, int nbrComponent = 0);

  // receive a vector from a MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the source MPI process
  // return value = reference on the current vector
  Vector& ReceiveVector(MPI::Intracomm& communicator, int id);

  // add current vector to the current vector of a given MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current vector
  Vector& SumVector(MPI::Intracomm& communicator, int id);

#endif

};

// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline double& DelocalizedRealVector::operator [] (int i)
{
  if (this->LocalizationId == this->LocalId)
    {
      return this->Components[i];
    }
  else
    {
      return this->Architecture->RequestRealVectorElement(this->VectorId, i);
    }
}
 

#endif

