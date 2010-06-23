////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of abstract hamiltonian                      //
//                                                                            //
//                        last modification : 28/02/2001                      //
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


#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "MathTools/Complex.h"
#include "BitmapTools/BitmapPicture/AbstractBitmapPicture.h"
#include "BitmapTools/BitmapPicture/TgaFormat.h" 
#include "BitmapTools/Color/PicRGB.h"


#include <stdlib.h>


using std::cout;
using std::endl;


// destructor
//

AbstractHamiltonian::~AbstractHamiltonian() 
{
}

// save precalculations in a file
//
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool AbstractHamiltonian::SavePrecalculation (char* fileName)
{
  return false;
}

// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

HermitianMatrix& AbstractHamiltonian::GetHamiltonian (HermitianMatrix& M)
{
  ComplexVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      if (this->IsHermitian())
	this->HermitianLowLevelMultiply(TmpV1, TmpV2);
      else
	this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = i; j < this->GetHilbertSpaceDimension(); j++)
	{
	  M.SetMatrixElement(i, j, TmpV2[j]);
	}
      TmpV1[i] = 0.0;
    }
  return M;  
}
  
// store real part of Hamiltonian into a real symmetric matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding real symmetric matrix 

RealSymmetricMatrix& AbstractHamiltonian::GetHamiltonian (RealSymmetricMatrix& M)
{
  RealVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      if (this->IsHermitian())
	this->HermitianLowLevelMultiply(TmpV1, TmpV2);
      else
	this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = i; j < this->GetHilbertSpaceDimension(); j++)
	M.SetMatrixElement(i, j, TmpV2[j]);
      TmpV1[i] = 0.0;	
    }
  return M;
}
  
// store Hamiltonian into a matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding matrix 

Matrix& AbstractHamiltonian::GetHamiltonian (Matrix& M)
{
  switch (M.GetMatrixType())
    {
      case (Matrix::RealElements | Matrix::Symmetric):
	return this->GetHamiltonian((RealSymmetricMatrix&) M);
      break;
      case (Matrix::ComplexElements | Matrix::Hermitian):
	return this->GetHamiltonian((HermitianMatrix&) M);
      break;
    default:
      return M;
      break;
    }
  return M;
}
  
// return matrix representation of current Hamiltonian
//
// return value = reference to representation

Matrix* AbstractHamiltonian::GetHamiltonian ()
{
  HermitianMatrix* TmpH = new HermitianMatrix(this->GetHilbertSpaceDimension());
  this->GetHamiltonian(*TmpH);
  return TmpH;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
{
  RealVector TmpVector (V2.GetVectorDimension(), true);
  this->LowLevelMultiply(V2, TmpVector);
  Complex Tmp = V1 * TmpVector;
  return Tmp;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  ComplexVector TmpVector (V2.GetVectorDimension(), true);
  this->LowLevelMultiply(V2, TmpVector);
  Complex Tmp = V1 * TmpVector;
  return Tmp;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AbstractHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination)
{
  if (this->IsHermitian())
    return this->HermitianLowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
  else
    return this->LowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						  int firstComponent, int nbrComponent) 
{
  vDestination.ClearVector();
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::LowLevelMultiply"<<endl;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& AbstractHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::LowLevelAddMultiply"<<endl;
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  if (this->IsHermitian())
    return this->HermitianLowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
  else
    return this->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
                                                                     int firstComponent, int nbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    vDestinations[i].ClearVector();
  return LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// multiply a set of vector by the current hamiltonian for a given range of indices 
// and store result in another set of vector, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->LowLevelMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			   destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vector sat which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
							     int firstComponent, int nbrComponent)
{
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}
 
// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result in another set of vectors, low level function (no architecture optimization)
//
// vSource = array of vectors to be multiplied
// vDestination = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
							     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->LowLevelAddMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			      destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  if (this->IsHermitian())
    return this->HermitianLowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
  else
    return this->LowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  vDestination.ClearVector();
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::LowLevelMultiply"<<endl;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int firstComponent, int nbrComponent)
{
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::LowLevelAddMultiply"<<endl;
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  if (this->IsHermitian())
    return this->HermitianLowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
  else
    return this->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							     int firstComponent, int nbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    vDestinations[i].ClearVector();
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->LowLevelMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			   destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								int firstComponent, int nbrComponent)
{
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}
 

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
								int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->LowLevelAddMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			      destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}




// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AbstractHamiltonian::ConjugateLowLevelMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->ConjugateLowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::ConjugateLowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
							   int firstComponent, int nbrComponent)
{
  vDestination.ClearVectorSegment((long)firstComponent, (long)nbrComponent);
  return this->ConjugateLowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::ConjugateLowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::ConjugateLowLevelMultiply"<<endl;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& AbstractHamiltonian::ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->ConjugateLowLevelAddMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  return this->ConjugateLowLevelAddMultiply(vSource, vDestination, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::ConjugateLowLevelAddMultiply"<<endl;
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::ConjugateLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  return this->ConjugateLowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::ConjugateLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							  int firstComponent, int nbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    vDestinations[i].ClearVector();
  return ConjugateLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// multiply a set of vector by the current hamiltonian for a given range of indices 
// and store result in another set of vector, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::ConjugateLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->ConjugateLowLevelMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			   destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vector sat which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::ConjugateLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  return this->ConjugateLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::ConjugateLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
							     int firstComponent, int nbrComponent)
{
  return this->ConjugateLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}
 
// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result in another set of vectors, low level function (no architecture optimization)
//
// vSource = array of vectors to be multiplied
// vDestination = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::ConjugateLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
							     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->ConjugateLowLevelAddMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			      destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractHamiltonian::ConjugateLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->ConjugateLowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::ConjugateLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  vDestination.ClearVectorSegment((long)firstComponent, (long)nbrComponent);
  return this->ConjugateLowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::ConjugateLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::ConjugateLowLevelMultiply"<<endl;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractHamiltonian::ConjugateLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->ConjugateLowLevelAddMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::ConjugateLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int firstComponent, int nbrComponent)
{
  return this->ConjugateLowLevelAddMultiply(vSource, vDestination, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::ConjugateLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::ConjugateLowLevelAddMultiply"<<endl;
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::ConjugateLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  return this->ConjugateLowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::ConjugateLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							     int firstComponent, int nbrComponent)
{
  for (int i=0; i<nbrVectors; ++i)
    vDestinations[i].ClearVectorSegment((long)firstComponent, (long)nbrComponent);
  return this->ConjugateLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::ConjugateLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->ConjugateLowLevelMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
				    destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::ConjugateLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  return this->ConjugateLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::ConjugateLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								int firstComponent, int nbrComponent)
{
  return this->ConjugateLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}
 

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::ConjugateLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
								int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->ConjugateLowLevelAddMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			      destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}



// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AbstractHamiltonian::HermitianLowLevelMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->HermitianLowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::HermitianLowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						  int firstComponent, int nbrComponent)
{
  vDestination.ClearVector();
  return this->HermitianLowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::HermitianLowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::HermitianowLevelMultiply"<<endl;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& AbstractHamiltonian::HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->HermitianLowLevelAddMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  return this->HermitianLowLevelAddMultiply(vSource, vDestination, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

RealVector& AbstractHamiltonian::HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::HermitianLowLevelAddMultiply"<<endl;
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::HermitianLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  return this->HermitianLowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::HermitianLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							  int firstComponent, int nbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    vDestinations[i].ClearVector();
  return HermitianLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// multiply a set of vector by the current hamiltonian for a given range of indices 
// and store result in another set of vector, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::HermitianLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							  int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							  int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->HermitianLowLevelMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			   destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vector sat which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  return this->HermitianLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
							     int firstComponent, int nbrComponent)
{
  return this->HermitianLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}
 
// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result in another set of vectors, low level function (no architecture optimization)
//
// vSource = array of vectors to be multiplied
// vDestination = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractHamiltonian::HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors,
							     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->HermitianLowLevelAddMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			      destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractHamiltonian::HermitianLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->HermitianLowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::HermitianLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  vDestination.ClearVector();
  return this->HermitianLowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::HermitianLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
						     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::HermitianLowLevelMultiply"<<endl;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->HermitianLowLevelAddMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int firstComponent, int nbrComponent)
{
  return this->HermitianLowLevelAddMultiply(vSource, vDestination, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  cout << "Attention, using dummy method AbstractHamiltonian::HermitianLowLevelAddMultiply"<<endl;
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::HermitianLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  return this->HermitianLowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::HermitianLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							     int firstComponent, int nbrComponent)
{
  for (int i=0; i<nbrVectors; ++i)
    vDestinations[i].ClearVector();
  return this->HermitianLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::HermitianLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							     int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
							     int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->HermitianLowLevelMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			   destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  return this->HermitianLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, 0, this->GetHilbertSpaceDimension());
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								int firstComponent, int nbrComponent)
{
  return this->HermitianLowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, 1, 0, nbrComponent, 0, 1, 0, this->GetHilbertSpaceDimension());
}
 

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
								int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->HermitianLowLevelAddMultiply(vSources[i], vDestinations[i], sourceStart, sourceStep, sourceShift, sourceNbrComponent, sourceNbrComponent, 
			      destinationStep, destinationShift, destinationNbrComponent);
  return vDestinations;
}


// store Hamiltonian into a picture (drawing non zero element with a color scale)
//
// error = absolute minimum value to be considered as non zero element
// return value = pointer to the picture associated to the matrix

AbstractBitmapPicture* AbstractHamiltonian::GetHamiltonianPicture (double error)
{
  RealVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  Color BlackColor (0.0, 0.0, 0.0);
  Color WhiteColor (1.0, 1.0, 1.0);
  TgaFormat* TmpPicture = new TgaFormat (this->GetHilbertSpaceDimension(), this->GetHilbertSpaceDimension());
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < this->GetHilbertSpaceDimension(); j++)
	{
	  if (fabs(TmpV2[j]) < error)
	    {
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, WhiteColor);
	    }
	  else
	    {
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, BlackColor);
	    }
	}
      TmpV1[i] = 0.0;	
    }
  return TmpPicture;
}
  
// store Hamiltonian into a picture (drawing non zero element in black)
//
// error = absolute minimum value to be considered as non zero element
// return value = pointer to the picture associated to the matrix

AbstractBitmapPicture* AbstractHamiltonian::GetHamiltonianColorPicture (double error)
{
  RealVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  Color RedColor (1.0, 0.0, 0.0);
  Color GreenColor (0.0, 1.0, 0.0);
  Color BlueColor (0.0, 0.0, 1.0);
  Color BlackColor (0.0, 0.0, 0.0);
  TgaFormat* TmpPicture = new TgaFormat (this->GetHilbertSpaceDimension(), this->GetHilbertSpaceDimension());
  double Max = 0.0;
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < this->GetHilbertSpaceDimension(); j++)
	{
	  if (fabs(TmpV2[j]) > error)
	    {
	      if (fabs(TmpV2[j]) > Max)
		{
		  Max = fabs(TmpV2[j]);
		}
	    }
	}
      TmpV1[i] = 0.0;	
    }
  Max = 4.0 / Max;
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < this->GetHilbertSpaceDimension(); j++)
	{
	  if (fabs(TmpV2[j]) < error)
	    {
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, BlackColor);
	    }
	  else
	    {
	      double Fac = (fabs(TmpV2[j]) * Max);
	      Color TmpColor;
	      if (Fac >= 4.0)
		{
		  TmpColor = RedColor;
		}
	      else
		if (Fac >= 3.0)
		  {
		    TmpColor = GreenColor * (4.0 - Fac) + RedColor;
		  }
		else
		  if (Fac >= 2.0)
		    {
		      TmpColor = GreenColor + RedColor * (Fac - 2.0);
		    }
		  else
		    if (Fac >= 1.0)
		      {
			TmpColor = BlueColor * (2.0 - Fac) + GreenColor;
		      }
		    else
		      {
			TmpColor = GreenColor * (Fac - 1.0) + BlueColor;
		      }
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, TmpColor);
	    }
	}
      TmpV1[i] = 0.0;	
    }
  return TmpPicture;
}
  
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractHamiltonian::LeftInteractionOperators()
{
  return List<Matrix*>();
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractHamiltonian::RightInteractionOperators()
{
  return List<Matrix*>();
}

// get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
// nbrThreads = number of threads requested
// segmentIndices = array returning the reference to an array of the first index of each of the segments
//
bool AbstractHamiltonian::GetLoadBalancing(int nbrTasks, long* &segmentIndices)
{
  return false;
}

// set the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
// nbrThreads = number of threads requested
// segmentIndices = array returning the first index of each of the segments
//
bool AbstractHamiltonian::SetLoadBalancing(int nbrTasks, long* segmentIndices)
{
  return false;
}

// ask if Hamiltonian implements methods using hermitian symmetry 
//
bool AbstractHamiltonian::IsHermitian()
{
  return false;
}

// ask if Hamiltonian implements methods applying the conjugate of the Hamiltonian
//
bool AbstractHamiltonian::IsConjugate()
{
  return false;
}


// multiply a vector by the current hamiltonian and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

Vector& AbstractHamiltonian::Multiply(Vector& vSource, Vector& vDestination)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    {
      return vDestination;
    }
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelMultiply((RealVector&) vSource, (RealVector&) vDestination);
    }
  else
    {
      return this->LowLevelMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination);
    }
  return vDestination;
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractHamiltonian::Multiply(Vector& vSource, Vector& vDestination, 
				      int firstComponent, int nbrComponent)
{
   if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    return vDestination;
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelMultiply((RealVector&) vSource, (RealVector&) vDestination, firstComponent, nbrComponent);
    }
  else
    {
      return this->LowLevelMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination, firstComponent, nbrComponent);
    }
  return vDestination;
}


// multiply a vector by the current hamiltonian and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

Vector& AbstractHamiltonian::ConjugateMultiply(Vector& vSource, Vector& vDestination)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    {
      return vDestination;
    }
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->ConjugateLowLevelMultiply((RealVector&) vSource, (RealVector&) vDestination);
    }
  else
    {
      return this->ConjugateLowLevelMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination);
    }
  return vDestination;
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractHamiltonian::ConjugateMultiply(Vector& vSource, Vector& vDestination, 
				      int firstComponent, int nbrComponent)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    return vDestination;
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->ConjugateLowLevelMultiply((RealVector&) vSource, (RealVector&) vDestination, firstComponent, nbrComponent);
    }
  else
    {
      return this->ConjugateLowLevelMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination, firstComponent, nbrComponent);
    }
  return vDestination;
}


// multiply a vector by the current hamiltonian and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

Vector& AbstractHamiltonian::HermitianMultiply(Vector& vSource, Vector& vDestination)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    {
      return vDestination;
    }
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->HermitianLowLevelMultiply((RealVector&) vSource, (RealVector&) vDestination);
    }
  else
    {
      return this->HermitianLowLevelMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination);
    }
  return vDestination;
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractHamiltonian::HermitianMultiply(Vector& vSource, Vector& vDestination, 
				      int firstComponent, int nbrComponent)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    return vDestination;
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->HermitianLowLevelMultiply((RealVector&) vSource, (RealVector&) vDestination, firstComponent, nbrComponent);
    }
  else
    {
      return this->HermitianLowLevelMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination, firstComponent, nbrComponent);
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

Vector& AbstractHamiltonian::AddMultiply(Vector& vSource, Vector& vDestination)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    {
      return vDestination;
    }
  if (vSource.GetVectorType() == Vector::RealDatas)
    {
      return this->LowLevelAddMultiply((RealVector&) vSource, (RealVector&) vDestination);
    }
  else
    {
      return this->LowLevelAddMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination);
    }
  return vDestination;
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractHamiltonian::AddMultiply(Vector& vSource, Vector& vDestination, 
					 int firstComponent, int nbrComponent)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    return vDestination;
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelAddMultiply((RealVector&) vSource, (RealVector&) vDestination, firstComponent, nbrComponent);
    }
  else
    {
      return this->LowLevelAddMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination, firstComponent, nbrComponent);
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::MultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelMultipleMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors);
    }
  else
    {
      return this->LowLevelMultipleMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::MultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors, 
					      int firstComponent, int nbrComponent)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  cout << "check 0" << endl;
  for (int i = 1; i < nbrVectors; ++i)
    {
      if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
	return vDestinations;
    }
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      cout << "check real " << endl;
      return this->LowLevelMultipleMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
  else
    {
      return this->LowLevelMultipleMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::MultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelMultipleAddMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors);
    }
  else
    {
      return this->LowLevelMultipleAddMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::MultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors,
						 int firstComponent, int nbrComponent)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->LowLevelMultipleAddMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
  else
    {
      return this->LowLevelMultipleAddMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

Vector& AbstractHamiltonian::ConjugateAddMultiply(Vector& vSource, Vector& vDestination)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    {
      return vDestination;
    }
  if (vSource.GetVectorType() == Vector::RealDatas)
    {
      return this->ConjugateLowLevelAddMultiply((RealVector&) vSource, (RealVector&) vDestination);
    }
  else
    {
      return this->ConjugateLowLevelAddMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination);
    }
  return vDestination;
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractHamiltonian::ConjugateAddMultiply(Vector& vSource, Vector& vDestination, 
					 int firstComponent, int nbrComponent)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    return vDestination;
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->ConjugateLowLevelAddMultiply((RealVector&) vSource, (RealVector&) vDestination, firstComponent, nbrComponent);
    }
  else
    {
      return this->ConjugateLowLevelAddMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination, firstComponent, nbrComponent);
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::ConjugateMultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->ConjugateLowLevelMultipleMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors);
    }
  else
    {
      return this->ConjugateLowLevelMultipleMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::ConjugateMultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors, 
					      int firstComponent, int nbrComponent)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  cout << "check 0" << endl;
  for (int i = 1; i < nbrVectors; ++i)
    {
      if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
	return vDestinations;
    }
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      cout << "check real " << endl;
      return this->ConjugateLowLevelMultipleMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
  else
    {
      return this->ConjugateLowLevelMultipleMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::ConjugateMultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->ConjugateLowLevelMultipleAddMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors);
    }
  else
    {
      return this->ConjugateLowLevelMultipleAddMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::ConjugateMultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors,
						 int firstComponent, int nbrComponent)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->ConjugateLowLevelMultipleAddMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
  else
    {
      return this->ConjugateLowLevelMultipleAddMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

Vector& AbstractHamiltonian::HermitianAddMultiply(Vector& vSource, Vector& vDestination)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    {
      return vDestination;
    }
  if (vSource.GetVectorType() == Vector::RealDatas)
    {
      return this->HermitianLowLevelAddMultiply((RealVector&) vSource, (RealVector&) vDestination);
    }
  else
    {
      return this->HermitianLowLevelAddMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination);
    }
  return vDestination;
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractHamiltonian::HermitianAddMultiply(Vector& vSource, Vector& vDestination, 
					 int firstComponent, int nbrComponent)
{
  if ((vSource.GetVectorType() & Vector::DataTypeMask) != (vDestination.GetVectorType() & Vector::DataTypeMask))
    return vDestination;
  if ((vSource.GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->HermitianLowLevelAddMultiply((RealVector&) vSource, (RealVector&) vDestination, firstComponent, nbrComponent);
    }
  else
    {
      return this->HermitianLowLevelAddMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination, firstComponent, nbrComponent);
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::HermitianMultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->HermitianLowLevelMultipleMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors);
    }
  else
    {
      return this->HermitianLowLevelMultipleMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::HermitianMultipleMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors, 
					      int firstComponent, int nbrComponent)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  cout << "check 0" << endl;
  for (int i = 1; i < nbrVectors; ++i)
    {
      if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
	return vDestinations;
    }
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      cout << "check real " << endl;
      return this->HermitianLowLevelMultipleMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
  else
    {
      return this->HermitianLowLevelMultipleMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::HermitianMultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->HermitianLowLevelMultipleAddMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors);
    }
  else
    {
      return this->HermitianLowLevelMultipleAddMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors);
    }
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector* AbstractHamiltonian::HermitianMultipleAddMultiply(Vector* vSources, Vector* vDestinations, int nbrVectors,
						 int firstComponent, int nbrComponent)
{
  if (vSources[0].GetVectorType() != vDestinations[0].GetVectorType())
    return vDestinations;
  for (int i = 1; i < nbrVectors; ++i)
    if ((vSources[0].GetVectorType() != vSources[i].GetVectorType()) || (vDestinations[0].GetVectorType() != vDestinations[i].GetVectorType()))
      return vDestinations;
  if ((vSources[0].GetVectorType() & Vector::DataTypeMask) == Vector::RealDatas)
    {
      return this->HermitianLowLevelMultipleAddMultiply((RealVector*) vSources, (RealVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
  else
    {
      return this->HermitianLowLevelMultipleAddMultiply((ComplexVector*) vSources, (ComplexVector*) vDestinations, nbrVectors, firstComponent, nbrComponent);
    }
}
