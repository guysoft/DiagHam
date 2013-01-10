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
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/GenericHamiltonianPrecalculationOperation.h"


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
// architecture = architecture to use for precalculation
// memory = amount of memory that can be used for precalculations (in bytes)

TensorProductSparseMatrixSelectedBlockHamiltonian::TensorProductSparseMatrixSelectedBlockHamiltonian(int nbrTensorProducts, SparseRealMatrix* leftMatrices,  
												     SparseRealMatrix* rightMatrices, double* coefficients,
												     int blockSize, long* blockIndices, AbstractArchitecture* architecture, long memory)
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
  this->HilbertSpace = new UndescribedHilbertSpace(blockSize);
  this->LeftHamiltonianVectorMultiplicationFlag = true;
  this->BlockIndices = blockIndices;
  long TmpMatrixSize = this->LeftMatrices[0].GetNbrRow();
  this->BlockIndexProductTable = new long* [TmpMatrixSize];
  this->BlockIndexProductTableNbrElements = new int [TmpMatrixSize];
  this->BlockIndexProductTableShift = new int [TmpMatrixSize];
  long* TmpBlockIndices = new long [TmpMatrixSize];
  for (long i = 0; i < TmpMatrixSize; ++i)
    {
      this->BlockIndexProductTableNbrElements[i] = 0;
      this->BlockIndexProductTableShift[i] = -1;
      for (int j = 0; j < blockSize; ++j)
	{
	  if ((this->BlockIndices[j] / TmpMatrixSize) == i)
	    {
	      if (this->BlockIndexProductTableShift[i] < 0)
		this->BlockIndexProductTableShift[i] = j;
	      TmpBlockIndices[this->BlockIndexProductTableNbrElements[i]] = this->BlockIndices[j];
	      ++this->BlockIndexProductTableNbrElements[i];
	    }
	}
      this->BlockIndexProductTable[i] = new long[this->BlockIndexProductTableNbrElements[i]];
      for (int j = 0; j < this->BlockIndexProductTableNbrElements[i]; ++j)
	{
	  this->BlockIndexProductTable[i][j] = TmpBlockIndices[j];
	}
    }
  delete[] TmpBlockIndices;
  this->BlockSize = blockSize;
  this->Architecture = architecture;
  if (memory > 0l)
    {
      this->TemporaryRowPointers = new long[this->HilbertSpace->GetHilbertSpaceDimension()];
      long NbrNonZeroMatrixElements = this->FastMultiplicationMemory() >> 3;
      cout << "nbr non-zero matrix elements = " << NbrNonZeroMatrixElements << endl;
      if (memory > (NbrNonZeroMatrixElements << 3))
	{
	  this->TemporaryRowLastPointers = new long[this->HilbertSpace->GetHilbertSpaceDimension()];
	  this->FastMultiplicationFlag = true;
	  long TmpPointer = 0;
	  for (int i = 0; i < this->HilbertSpace->GetHilbertSpaceDimension(); ++i)
	    {
	      long Tmp = this->TemporaryRowPointers[i]; 
	      if (this->TemporaryRowPointers[i] > 0l)
		{     
		  this->TemporaryRowPointers[i] = TmpPointer;
		  this->TemporaryRowLastPointers[i] = TmpPointer + Tmp - 1l;
		}
	      else
		{
		  this->TemporaryRowPointers[i] = -1l;
		  this->TemporaryRowLastPointers[i] = -1l;
		}
	      TmpPointer += Tmp;
	    }
	  this->TemporaryMatrixElements = new double[NbrNonZeroMatrixElements];
	  this->TemporaryMatrixColumnIndices = new int[NbrNonZeroMatrixElements];
	  this->EnableFastMultiplication(); 
	  cout << "using ";
	  PrintMemorySize(cout, (NbrNonZeroMatrixElements << 3)) << " for hamiltonian precalculations" << endl;
	}
      else
	{
	  cout << "cannot use hamiltonian precalculations, ";
	  PrintMemorySize(cout, (NbrNonZeroMatrixElements << 3)) << " are required" << endl;
	  this->FastMultiplicationFlag = false;
	  delete[] this->TemporaryRowPointers;
	}
    }
  else
    {
      this->FastMultiplicationFlag = false;
    }
}

// contructor providing an efficient block index scheme
//
// nbrTensorProducts = number of tensor products whose linear combination defined the Hamiltonian 
// leftMatrices = left matrices of each tensor product
// rightMatrices = right matrices of each tensor product
// coefficients = coefficients of the ensor product linear combination
// blockSize = number of indices in the selected block
// blockIndices = pairs of indices (for resp. the left and right matrix) that define the selected block 
// blockIndexProductTable = table that contains all the linearized indices attached to one left matrix index
// blockIndexProductTableNbrElements = number of  linearized indices attached to one left matrix index
// blockIndexProductTableShift = first index in the hilbert space where a given left matrix index occurs
// architecture = architecture to use for precalculation
// memory = amount of memory that can be used for precalculations (in bytes)

TensorProductSparseMatrixSelectedBlockHamiltonian::TensorProductSparseMatrixSelectedBlockHamiltonian(int nbrTensorProducts, SparseRealMatrix* leftMatrices,  
												     SparseRealMatrix* rightMatrices, double* coefficients,
												     int blockSize, long* blockIndices, 
												     long** blockIndexProductTable, int* blockIndexProductTableNbrElements,
												     int* blockIndexProductTableShift,
												     AbstractArchitecture* architecture, long memory)
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
  this->HilbertSpace = new UndescribedHilbertSpace(blockSize);
  this->LeftHamiltonianVectorMultiplicationFlag = true;
  this->BlockIndices = blockIndices;
  long TmpMatrixSize = this->LeftMatrices[0].GetNbrRow();
  this->BlockIndexProductTable = blockIndexProductTable;
  this->BlockIndexProductTableNbrElements = blockIndexProductTableNbrElements;
  this->BlockIndexProductTableShift = blockIndexProductTableShift;
  this->BlockSize = blockSize;
  this->Architecture = architecture;
  if (memory > 0l)
    {
      this->TemporaryRowPointers = new long[this->HilbertSpace->GetHilbertSpaceDimension()];
      long NbrNonZeroMatrixElements = this->FastMultiplicationMemory() >> 3;
      cout << "nbr non-zero matrix elements = " << NbrNonZeroMatrixElements << endl;
      if (memory > (NbrNonZeroMatrixElements << 3))
	{
	  this->TemporaryRowLastPointers = new long[this->HilbertSpace->GetHilbertSpaceDimension()];
	  this->FastMultiplicationFlag = true;
	  long TmpPointer = 0;
	  for (int i = 0; i < this->HilbertSpace->GetHilbertSpaceDimension(); ++i)
	    {
	      long Tmp = this->TemporaryRowPointers[i]; 
	      if (this->TemporaryRowPointers[i] > 0l)
		{     
		  this->TemporaryRowPointers[i] = TmpPointer;
		  this->TemporaryRowLastPointers[i] = TmpPointer + Tmp - 1l;
		}
	      else
		{
		  this->TemporaryRowPointers[i] = -1l;
		  this->TemporaryRowLastPointers[i] = -1l;
		}
	      TmpPointer += Tmp;
	    }
	  this->TemporaryMatrixElements = new double[NbrNonZeroMatrixElements];
	  this->TemporaryMatrixColumnIndices = new int[NbrNonZeroMatrixElements];
	  this->EnableFastMultiplication(); 
	  cout << "using ";
	  PrintMemorySize(cout, (NbrNonZeroMatrixElements << 3)) << " for hamiltonian precalculations" << endl;
	}
      else
	{
	  cout << "cannot use hamiltonian precalculations, ";
	  PrintMemorySize(cout, (NbrNonZeroMatrixElements << 3)) << " are required" << endl;
	  this->FastMultiplicationFlag = false;
	  delete[] this->TemporaryRowPointers;
	}
    }
  else
    {
      this->FastMultiplicationFlag = false;
    }
}

// destructor
//

TensorProductSparseMatrixSelectedBlockHamiltonian::~TensorProductSparseMatrixSelectedBlockHamiltonian() 
{
  long TmpMatrixSize = this->LeftMatrices[0].GetNbrRow();
  for (long i = 0; i < TmpMatrixSize; ++i)
    delete[] this->BlockIndexProductTable[i];
  delete[] this->BlockIndexProductTable;
  delete[] this->BlockIndexProductTableNbrElements;
  delete[] this->BlockIndexProductTableShift;
  if (this->FastMultiplicationFlag == true)
    {
      delete[] this->TemporaryMatrixElements;
      delete[] this->TemporaryMatrixColumnIndices;
      delete[] this->TemporaryRowPointers;
      delete[] this->TemporaryRowLastPointers;
    }
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
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == true)
    {
      for (int j = firstComponent; j < LastComponent; ++j)
	{
	  long TmpRowPointer = this->TemporaryRowPointers[j];
	  if (TmpRowPointer >= 0l)
	    {
	      long TmpRowLastPointer = this->TemporaryRowLastPointers[j];
	      double Tmp = 0.0;
	      for (; TmpRowPointer <= TmpRowLastPointer; ++TmpRowPointer)
		{
		  Tmp += (this->TemporaryMatrixElements[TmpRowPointer] 
			  * vSource[this->TemporaryMatrixColumnIndices[TmpRowPointer]]);
		}		
	      vDestination[j] += Tmp;
	    }
	}
      return vDestination;
    }
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
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
		      int TmpLeftMatrixColumnIndex = TmpLeftMatrix.ColumnIndices[k];
		      int LocalBlockSize = this->BlockIndexProductTableNbrElements[TmpLeftMatrixColumnIndex];
		      if (LocalBlockSize > 0)
			{
			  double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  long TmpIndex = ((long) TmpLeftMatrixColumnIndex) * IndexStep;
			  long* LocalBlockIndices = this->BlockIndexProductTable[TmpLeftMatrixColumnIndex];
			  int LocalShift = this->BlockIndexProductTableShift[TmpLeftMatrixColumnIndex];
			  for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
			    {
			      int TmpIndex2 = TmpIndex + TmpRightMatrix.ColumnIndices[l];
			      int TmpIndex3 = SearchInArray<long>(TmpIndex + TmpRightMatrix.ColumnIndices[l], LocalBlockIndices, LocalBlockSize);
			      if (TmpIndex3 >= 0)
				Tmp += Tmp2 * TmpRightMatrix.MatrixElements[l] * vSource[LocalShift + TmpIndex3];
			    }
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
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == true)
    {
      for (int j = firstComponent; j < LastComponent; ++j)
	{
	  long TmpRowPointer = this->TemporaryRowPointers[j];
	  if (TmpRowPointer >= 0l)
	    {
	      long TmpRowLastPointer = this->TemporaryRowLastPointers[j];
	      Complex Tmp = 0.0;
	      for (; TmpRowPointer <= TmpRowLastPointer; ++TmpRowPointer)
		{
		  Tmp += (this->TemporaryMatrixElements[TmpRowPointer] 
			  * vSource[this->TemporaryMatrixColumnIndices[TmpRowPointer]]);
		}		
	      vDestination[j] += Tmp;
	    }
	}
      return vDestination;
    }
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
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
		      int TmpLeftMatrixColumnIndex = TmpLeftMatrix.ColumnIndices[k];
		      int LocalBlockSize = this->BlockIndexProductTableNbrElements[TmpLeftMatrixColumnIndex];
		      if (LocalBlockSize > 0)
			{
			  double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  long TmpIndex = ((long) TmpLeftMatrixColumnIndex) * IndexStep;
			  long* LocalBlockIndices = this->BlockIndexProductTable[TmpLeftMatrixColumnIndex];
			  int LocalShift = this->BlockIndexProductTableShift[TmpLeftMatrixColumnIndex];
			  for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
			    {
			      int TmpIndex2 = TmpIndex + TmpRightMatrix.ColumnIndices[l];
			      int TmpIndex3 = SearchInArray<long>(TmpIndex + TmpRightMatrix.ColumnIndices[l], LocalBlockIndices, LocalBlockSize);
			      if (TmpIndex3 >= 0)
				Tmp += Tmp2 * TmpRightMatrix.MatrixElements[l] * vSource[LocalShift + TmpIndex3];
			    }
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

ComplexVector* TensorProductSparseMatrixSelectedBlockHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, 
											      int nbrVectors, int firstComponent, int nbrComponent)
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

// test the amount of memory needed for fast multiplication algorithm
//
// return value = amount of memory needed

long TensorProductSparseMatrixSelectedBlockHamiltonian::FastMultiplicationMemory()
{
  GenericHamiltonianPrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);
  long NbrNonZeroMatrixElements = 0l;
  for (int i = 0; i < this->HilbertSpace->GetHilbertSpaceDimension(); ++i)
    {
      NbrNonZeroMatrixElements += this->TemporaryRowPointers[i];
    }
  return (NbrNonZeroMatrixElements << 3);
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that has to be precalcualted
// return value = number of non-zero matrix elements that have to be stored

long TensorProductSparseMatrixSelectedBlockHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent;
  int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  long NbrNonZeroMatrixElements = 0l;
  int TmpDimension = this->HilbertSpace->GetHilbertSpaceDimension();
  double* TmpNonZeroMatrixElements = new double [TmpDimension];
  for (int j = firstComponent; j < LastComponent; ++j)
    {
      for (int i = 0; i < TmpDimension; ++i)
	{
	  TmpNonZeroMatrixElements[i] = 0.0;
	}
      for (int i = 0; i < this->NbrTensorProducts; ++i)
	{
	  SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
	  SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
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
		      int TmpLeftMatrixColumnIndex = TmpLeftMatrix.ColumnIndices[k];
		      int LocalBlockSize = this->BlockIndexProductTableNbrElements[TmpLeftMatrixColumnIndex];
		      if (LocalBlockSize > 0)
			{
			  double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  long TmpIndex = ((long) TmpLeftMatrixColumnIndex) * IndexStep;
			  long* LocalBlockIndices = this->BlockIndexProductTable[TmpLeftMatrixColumnIndex];
			  int LocalShift = this->BlockIndexProductTableShift[TmpLeftMatrixColumnIndex];
			  for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
			    {
			      int TmpIndex2 = TmpIndex + TmpRightMatrix.ColumnIndices[l];
			      int TmpIndex3 = SearchInArray<long>(TmpIndex + TmpRightMatrix.ColumnIndices[l], LocalBlockIndices, LocalBlockSize);
			      if (TmpIndex3 >= 0)
				{
				  TmpNonZeroMatrixElements[LocalShift + TmpIndex3] += Tmp2 * TmpRightMatrix.MatrixElements[l];				  
				}
			    }
			}
		    }
		}
	    }
	}
      if (this->HamiltonianShift != 0.0)
	{
	  TmpNonZeroMatrixElements[j] += this->HamiltonianShift;
	}
      long Tmp = 0;
      for (int i = 0; i < TmpDimension; ++i)
	{
	  if (TmpNonZeroMatrixElements[i] != 0.0)
	    ++Tmp;
	}
      this->TemporaryRowPointers[j] = Tmp;
      NbrNonZeroMatrixElements += Tmp;
    }
  delete[] TmpNonZeroMatrixElements;
  return NbrNonZeroMatrixElements;
}

// enable fast multiplication algorithm
//
  
void TensorProductSparseMatrixSelectedBlockHamiltonian::EnableFastMultiplication()
{
  GenericHamiltonianPrecalculationOperation Operation(this, false);
  Operation.ApplyOperation(this->Architecture);
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that has to be precalcualted

void TensorProductSparseMatrixSelectedBlockHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{
  int IndexStep = this->LeftMatrices[0].GetNbrColumn();
  int LastComponent = firstComponent + nbrComponent;
  int AMatrixLastIndex = LastComponent / this->LeftMatrices[0].GetNbrRow();
  int BMatrixLastIndex = LastComponent % this->LeftMatrices[0].GetNbrRow();
  long TmpARowPointer;
  long TmpARowLastPointer;
  long TmpBRowPointer;
  long TmpBRowLastPointer;
  long NbrNonZeroMatrixElements = 0l;
  int TmpDimension = this->HilbertSpace->GetHilbertSpaceDimension();
  double* TmpNonZeroMatrixElements = new double [TmpDimension];
  for (int j = firstComponent; j < LastComponent; ++j)
    {
      for (int i = 0; i < TmpDimension; ++i)
	{
	  TmpNonZeroMatrixElements[i] = 0.0;
	}
      for (int i = 0; i < this->NbrTensorProducts; ++i)
	{
	  SparseRealMatrix& TmpLeftMatrix = this->LeftMatrices[i];
	  SparseRealMatrix& TmpRightMatrix = this->RightMatrices[i];
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
		      int TmpLeftMatrixColumnIndex = TmpLeftMatrix.ColumnIndices[k];
		      int LocalBlockSize = this->BlockIndexProductTableNbrElements[TmpLeftMatrixColumnIndex];
		      if (LocalBlockSize > 0)
			{
			  double Tmp2 = TmpLeftMatrix.MatrixElements[k] * this->Coefficients[i];
			  long TmpIndex = ((long) TmpLeftMatrixColumnIndex) * IndexStep;
			  long* LocalBlockIndices = this->BlockIndexProductTable[TmpLeftMatrixColumnIndex];
			  int LocalShift = this->BlockIndexProductTableShift[TmpLeftMatrixColumnIndex];
			  for (long l = TmpBRowPointer; l <= TmpBRowLastPointer; ++l)
			    {
			      int TmpIndex2 = TmpIndex + TmpRightMatrix.ColumnIndices[l];
			      int TmpIndex3 = SearchInArray<long>(TmpIndex + TmpRightMatrix.ColumnIndices[l], 
								  LocalBlockIndices, LocalBlockSize);
			      if (TmpIndex3 >= 0)
				{
				  TmpNonZeroMatrixElements[LocalShift + TmpIndex3] += Tmp2 * TmpRightMatrix.MatrixElements[l];				  
				}
			    }
			}
		    }
		}
	    }
	}
      if (this->HamiltonianShift != 0.0)
	{
	  TmpNonZeroMatrixElements[j] += this->HamiltonianShift;
	}
       long Tmp = 0;
       long Shift = this->TemporaryRowPointers[j];
       for (int i = 0; i < TmpDimension; ++i)
	 {
	   if (TmpNonZeroMatrixElements[i] != 0.0)
	     {
	       this->TemporaryMatrixElements[Shift] = TmpNonZeroMatrixElements[i];
	       this->TemporaryMatrixColumnIndices[Shift] = i;
	       ++Shift;
	     }
	 }     
       NbrNonZeroMatrixElements += Shift - this->TemporaryRowPointers[j];
    }
  delete[] TmpNonZeroMatrixElements;
}
  
