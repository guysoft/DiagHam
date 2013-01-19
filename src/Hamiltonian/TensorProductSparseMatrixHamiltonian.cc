////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of hamiltonian defined as a linear combination of tensor products   //
//                                                                            //
//                        last modification : 08/11/2012                      //
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


#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"
#include "MathTools/Complex.h" 
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"


#include <iostream>


using std::cout;
using std::endl;


// default contructor 
//

TensorProductSparseMatrixHamiltonian::TensorProductSparseMatrixHamiltonian()
{
}

// contructor 
//
// nbrTensorProducts = number of tensor products whose linear combination defined the Hamiltonian 
// leftMatrices = left matrices of each tensor product
// rightMatrices = right matrices of each tensor product
// coefficients = coefficients of the ensor product linear combination

TensorProductSparseMatrixHamiltonian::TensorProductSparseMatrixHamiltonian(int nbrTensorProducts, SparseRealMatrix* leftMatrices,  SparseRealMatrix* rightMatrices, double* coefficients)
{
  this->NbrTensorProducts = nbrTensorProducts;
  this->LeftMatrices = new SparseRealMatrix[this->NbrTensorProducts];
  this->RightMatrices = new SparseRealMatrix[this->NbrTensorProducts];
  this->Coefficients = new double[this->NbrTensorProducts];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      this->LeftMatrices[i] = leftMatrices[i];
      this->RightMatrices[i] = rightMatrices[i];
      this->Coefficients[i] = coefficients[i];
    }
  this->HamiltonianShift = 0.0;
  long HamiltonianDimension = this->LeftMatrices[0].GetNbrRow() * this->RightMatrices[0].GetNbrRow();
  this->HilbertSpace = new UndescribedHilbertSpace(HamiltonianDimension);
  this->LeftHamiltonianVectorMultiplicationFlag = true;
}

// destructor
//

TensorProductSparseMatrixHamiltonian::~TensorProductSparseMatrixHamiltonian() 
{
  if (this->LeftMatrices != 0)
    {
      delete[] this->LeftMatrices;
      delete[] this->RightMatrices;
      delete[] this->Coefficients;
    }
  delete this->HilbertSpace;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TensorProductSparseMatrixHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace) 
{
  this->HilbertSpace = hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TensorProductSparseMatrixHamiltonian::GetHilbertSpace ()
{
  return this->HilbertSpace;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TensorProductSparseMatrixHamiltonian::GetHilbertSpaceDimension () 
{
  return this->HilbertSpace->GetHilbertSpaceDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void TensorProductSparseMatrixHamiltonian::ShiftHamiltonian (double shift) 
{
  this->HamiltonianShift = shift;
}


// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex TensorProductSparseMatrixHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  return Complex();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex TensorProductSparseMatrixHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& TensorProductSparseMatrixHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent - 1;
  int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
      int AMatrixStartingIndex = firstComponent / this->LeftMatrices[0].GetNbrRow();
      int BMatrixStartingIndex = firstComponent % this->LeftMatrices[0].GetNbrRow();
      int TotalIndex = firstComponent;
      for (; AMatrixStartingIndex <=  AMatrixLastIndex; ++AMatrixStartingIndex)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[AMatrixStartingIndex];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[AMatrixStartingIndex];
	      int TmpBMatrixLastIndex = TmpRightMatrix.GetNbrRow() - 1;
	      if (AMatrixStartingIndex == AMatrixLastIndex)
		TmpBMatrixLastIndex = BMatrixLastIndex;
	      for (; BMatrixStartingIndex <= TmpBMatrixLastIndex; ++BMatrixStartingIndex)
		{
		  TmpBRowPointer = TmpRightMatrix.RowPointers[BMatrixStartingIndex];
		  if (TmpBRowPointer >= 0l)
		    {
		      TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[BMatrixStartingIndex];
		      double Tmp= 0.0;
		      for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
			{
			  double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  int TmpIndex = TmpLeftMatrix.ColumnIndices[k] * IndexStep;
			  for (long j = TmpBRowPointer; j <= TmpBRowLastPointer; ++j)
			    {
			      Tmp += Tmp2 * TmpRightMatrix.MatrixElements[j] * vSource[TmpIndex + TmpRightMatrix.ColumnIndices[j]];
			    }
			}
		      vDestination[AMatrixStartingIndex * IndexStep + BMatrixStartingIndex] += Tmp;
		    }
		}
	    }
	  BMatrixStartingIndex = 0;
	}
    }
  if (this->HamiltonianShift != 0.0)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	vDestination[i] += this->HamiltonianShift * vSource[i];
    }
  return vDestination;
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

RealVector* TensorProductSparseMatrixHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent;
  int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  double* Tmp = new double[nbrVectors];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
      int AMatrixStartingIndex = firstComponent / this->LeftMatrices[0].GetNbrRow();
      int BMatrixStartingIndex = firstComponent % this->LeftMatrices[0].GetNbrRow();
      int TotalIndex = firstComponent;
      for (; AMatrixStartingIndex <  AMatrixLastIndex; ++AMatrixStartingIndex)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[AMatrixStartingIndex];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[AMatrixStartingIndex];
	      int TmpBMatrixLastIndex = TmpRightMatrix.GetNbrRow();
	      if (AMatrixStartingIndex == (AMatrixLastIndex - 1))
		TmpBMatrixLastIndex = BMatrixLastIndex;
	      for (; BMatrixStartingIndex <  TmpBMatrixLastIndex; ++BMatrixStartingIndex)
		{
		  TmpBRowPointer = TmpRightMatrix.RowPointers[BMatrixStartingIndex];
		  if (TmpBRowPointer >= 0l)
		    {
		      TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[BMatrixStartingIndex];
		      for (int l = 0; l < nbrVectors; ++l)
			Tmp[l] = 0.0;
		      for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
			{
			  double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  int TmpIndex = TmpLeftMatrix.ColumnIndices[k] * IndexStep;
			  for (long j = TmpBRowPointer; j <= TmpBRowLastPointer; ++j)
			    {
			      int InputIndex = TmpIndex + TmpRightMatrix.ColumnIndices[j];
			      for (int l = 0; l < nbrVectors; ++l)			      
				Tmp[l] += Tmp2 * TmpRightMatrix.MatrixElements[j] * vSources[l][InputIndex];
			    }
			}
		      int OutputIndex = AMatrixStartingIndex * IndexStep + BMatrixStartingIndex;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][OutputIndex] += Tmp[l];
		    }
		}
	    }
	  BMatrixStartingIndex = 0;
	}
    }
  delete[] Tmp;
  if (this->HamiltonianShift != 0.0)
    {
      for (int k= 0; k < nbrVectors; ++k)
	{
	  RealVector& TmpDestination = vDestinations[k];
	  RealVector& TmpSource = vSources[k];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
	}
    }
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& TensorProductSparseMatrixHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent;
  int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
      int AMatrixStartingIndex = firstComponent / this->LeftMatrices[0].GetNbrRow();
      int BMatrixStartingIndex = firstComponent % this->LeftMatrices[0].GetNbrRow();
      int TotalIndex = firstComponent;
      for (; AMatrixStartingIndex <  AMatrixLastIndex; ++AMatrixStartingIndex)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[AMatrixStartingIndex];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[AMatrixStartingIndex];
	      int TmpBMatrixLastIndex = TmpRightMatrix.GetNbrRow();
	      if (AMatrixStartingIndex == (AMatrixLastIndex - 1))
		TmpBMatrixLastIndex = BMatrixLastIndex;
	      for (; BMatrixStartingIndex <  TmpBMatrixLastIndex; ++BMatrixStartingIndex)
		{
		  TmpBRowPointer = TmpRightMatrix.RowPointers[BMatrixStartingIndex];
		  if (TmpBRowPointer >= 0l)
		    {
		      TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[BMatrixStartingIndex];
		      Complex Tmp= 0.0;
		      for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
			{
			  Complex Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  int TmpIndex = TmpLeftMatrix.ColumnIndices[k] * IndexStep;
			  for (long j = TmpBRowPointer; j <= TmpBRowLastPointer; ++j)
			    {
			      Tmp += Tmp2 * TmpRightMatrix.MatrixElements[j] * vSource[TmpIndex + TmpRightMatrix.ColumnIndices[j]];
			    }
			}
		      vDestination[AMatrixStartingIndex * IndexStep + BMatrixStartingIndex] += Tmp;
		    }
		}
	    }
	  BMatrixStartingIndex = 0;
	}
    }
  if (this->HamiltonianShift != 0.0)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	vDestination[i] += this->HamiltonianShift * vSource[i];
    }
  return vDestination;
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

ComplexVector* TensorProductSparseMatrixHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent;
  int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  Complex* Tmp = new Complex[nbrVectors];
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
      SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
      int AMatrixStartingIndex = firstComponent / this->LeftMatrices[0].GetNbrRow();
      int BMatrixStartingIndex = firstComponent % this->LeftMatrices[0].GetNbrRow();
      int TotalIndex = firstComponent;
      for (; AMatrixStartingIndex <  AMatrixLastIndex; ++AMatrixStartingIndex)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[AMatrixStartingIndex];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[AMatrixStartingIndex];
	      int TmpBMatrixLastIndex = TmpRightMatrix.GetNbrRow();
	      if (AMatrixStartingIndex == (AMatrixLastIndex - 1))
		TmpBMatrixLastIndex = BMatrixLastIndex;
	      for (; BMatrixStartingIndex <  TmpBMatrixLastIndex; ++BMatrixStartingIndex)
		{
		  TmpBRowPointer = TmpRightMatrix.RowPointers[BMatrixStartingIndex];
		  if (TmpBRowPointer >= 0l)
		    {
		      TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[BMatrixStartingIndex];
		      for (int l = 0; l < nbrVectors; ++l)
			Tmp[l] = 0.0;
		      for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
			{
			  Complex Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  int TmpIndex = TmpLeftMatrix.ColumnIndices[k] * IndexStep;
			  for (long j = TmpBRowPointer; j <= TmpBRowLastPointer; ++j)
			    {
			      int InputIndex = TmpIndex + TmpRightMatrix.ColumnIndices[j];
			      for (int l = 0; l < nbrVectors; ++l)			      
				Tmp[l] += Tmp2 * TmpRightMatrix.MatrixElements[j] * vSources[l][InputIndex];
			    }
			}
		      int OutputIndex = AMatrixStartingIndex * IndexStep + BMatrixStartingIndex;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][OutputIndex] += Tmp[l];
		    }
		}
	    }
	  BMatrixStartingIndex = 0;
	}
    }
  delete[] Tmp;
  if (this->HamiltonianShift != 0.0)
    {
      for (int k= 0; k < nbrVectors; ++k)
	{
	  ComplexVector& TmpDestination = vDestinations[k];
	  ComplexVector& TmpSource = vSources[k];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
	}
    }
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& TensorProductSparseMatrixHamiltonian::ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									       int firstComponent, int nbrComponent)
{
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
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

RealVector* TensorProductSparseMatrixHamiltonian::ConjugateLowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										    int firstComponent, int nbrComponent)
{
  return this->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& TensorProductSparseMatrixHamiltonian::ConjugateLowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									       int firstComponent, int nbrComponent)
{
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
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

ComplexVector* TensorProductSparseMatrixHamiltonian::ConjugateLowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										       int firstComponent, int nbrComponent)
{
  return this->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

