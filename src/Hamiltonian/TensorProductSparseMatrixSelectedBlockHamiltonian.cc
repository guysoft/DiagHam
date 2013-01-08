////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of hamiltonian defined as a linear combination of tensor products   //
//               focusing on a single block of  tensor product                //
//                                                                            //
//                        last modification : 07/01/2013                      //
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


#include "Hamiltonian/TensorProductSparseMatrixSelectedBlockHamiltonian.h"
#include "MathTools/Complex.h" 
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "GeneralTools/ArrayTools.h"


#include <iostream>


using std::cout;
using std::endl;


// contructor 
//
// nbrTensorProducts = number of tensor products whose linear combination defined the Hamiltonian 
// leftMatrices = left matrices of each tensor product
// rightMatrices = right matrices of each tensor product
// coefficients = coefficients of the ensor product linear combination
// blockSize = number of indices in the selected block
// blockIndices = pairs of indices (for resp. the left and right matrix) that define the selected block 

TensorProductSparseMatrixSelectedBlockHamiltonian::TensorProductSparseMatrixSelectedBlockHamiltonian(int nbrTensorProducts, SparseRealMatrix* leftMatrices,  
												     SparseRealMatrix* rightMatrices, double* coefficients,
												     long blockSize, int* blockIndices)
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
  this->HilbertSpace = new UndescribedHilbertSpace(blockSize);
  this->LeftHamiltonianVectorMultiplicationFlag = true;
  this->BlockIndices = blockIndices;
  this->BlockSize = blockSize;
}

// destructor
//

TensorProductSparseMatrixSelectedBlockHamiltonian::~TensorProductSparseMatrixSelectedBlockHamiltonian() 
{
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
      for (int j = firstComponent; j < LastComponent; ++j)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[this->BlockIndices[j] / IndexStep];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[this->BlockIndices[j] / IndexStep];
	      TmpBRowPointer = TmpRightMatrix.RowPointers[this->BlockIndices[j] % IndexStep];
	      if (TmpBRowPointer >= 0l)
		{
		  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[this->BlockIndices[j] % IndexStep];
		  double Tmp= 0.0;
		  for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
		    {
		      double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
		      int TmpIndex = TmpLeftMatrix.ColumnIndices[k] * IndexStep;
		      for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
			{
			  int TmpIndex2 = TmpIndex + TmpRightMatrix.ColumnIndices[l];
			  int TmpIndex3 = SearchInArray<int>((TmpIndex + TmpRightMatrix.ColumnIndices[l]), this->BlockIndices, this->BlockSize);
			  if (TmpIndex3 >= 0)
			    Tmp += Tmp2 * TmpRightMatrix.MatrixElements[l] * vSource[TmpIndex3];
			}
		    }
		  vDestination[j] += Tmp;
		}
	    }
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

RealVector* TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
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

ComplexVector& TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
      for (int j = firstComponent; j < LastComponent; ++j)
	{
	  TmpARowPointer = TmpLeftMatrix.RowPointers[this->BlockIndices[j] / IndexStep];
	  if (TmpARowPointer >= 0l)
	    {
	      TmpARowLastPointer = TmpLeftMatrix.RowLastPointers[this->BlockIndices[j] / IndexStep];
	      TmpBRowPointer = TmpRightMatrix.RowPointers[this->BlockIndices[j] % IndexStep];
	      if (TmpBRowPointer >= 0l)
		{
		  TmpBRowLastPointer = TmpRightMatrix.RowLastPointers[this->BlockIndices[j] % IndexStep];
		  Complex Tmp= 0.0;
		  for (long k = TmpARowPointer; k <= TmpARowLastPointer; ++k)
		    {
		      double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
		      int TmpIndex = TmpLeftMatrix.ColumnIndices[k] * IndexStep;
		      for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
			{
			  int TmpIndex2 = TmpIndex + TmpRightMatrix.ColumnIndices[l];
			  int TmpIndex3 = SearchInArray<int>((TmpIndex + TmpRightMatrix.ColumnIndices[l]), this->BlockIndices, this->BlockSize);
			  if (TmpIndex3 >= 0)
			    Tmp += Tmp2 * TmpRightMatrix.MatrixElements[l] * vSource[TmpIndex3];
			}
		    }
		  vDestination[j] += Tmp;
		}
	    }
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

ComplexVector* TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
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
