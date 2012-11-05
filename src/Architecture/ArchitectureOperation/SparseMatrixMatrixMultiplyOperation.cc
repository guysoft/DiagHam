////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of matrix matrix multiplication operation             //
//                                                                            //
//                        last modification : 15/01/2003                      //
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
#include "Architecture/ArchitectureOperation/SparseMatrixMatrixMultiplyOperation.h"
#include "Matrix/SparseRealMatrix.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"


// constructor 
//
// leftMatrix = pointer to the matrix used as left matrix for the multiplication and where the result will be stored
// rightMatrix = pointer to the matrix used as right matrix for the multiplication
// tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
// nbrTmpMatrixElements = maximum number of elements available in tmpMatrixElements

SparseMatrixMatrixMultiplyOperation::SparseMatrixMatrixMultiplyOperation (SparseRealMatrix* leftMatrix, SparseRealMatrix* rightMatrix,
									  double* tmpMatrixElements, int* tmpColumnIndices, 
									  long nbrTmpMatrixElements)
{
  this->FirstComponent = 0;
  this->NbrComponent = leftMatrix->GetNbrRow();
  this->LeftMatrix = leftMatrix;
  this->RightMatrix = rightMatrix;
  this->LocalTmpMatrixElements = tmpMatrixElements;
  this->LocalTmpColumnIndices = tmpColumnIndices;
  this->MaxTmpMatrixElements = nbrTmpMatrixElements;
  this->LocalNbrMatrixElements = 0l;
  this->DestinationMatrix = SparseRealMatrix (this->LeftMatrix->NbrRow, this->RightMatrix->NbrColumn, 0);
  this->OperationType = AbstractArchitectureOperation::SparseMatrixMatrixMultiply;
}

// copy constructor 
//
// operation = reference on operation to copy

SparseMatrixMatrixMultiplyOperation::SparseMatrixMatrixMultiplyOperation(const SparseMatrixMatrixMultiplyOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LeftMatrix = operation.LeftMatrix;
  this->RightMatrix = operation.RightMatrix;
  this->LocalTmpMatrixElements = operation.LocalTmpMatrixElements;
  this->LocalTmpColumnIndices = operation.LocalTmpColumnIndices;
  this->MaxTmpMatrixElements = operation.MaxTmpMatrixElements;
  this->LocalNbrMatrixElements = operation.LocalNbrMatrixElements;
  this->DestinationMatrix = operation.DestinationMatrix;
  this->OperationType = AbstractArchitectureOperation::SparseMatrixMatrixMultiply;
}
  
// destructor
//

SparseMatrixMatrixMultiplyOperation::~SparseMatrixMatrixMultiplyOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void SparseMatrixMatrixMultiplyOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set the temporary arrays required during the matrix product
//
// tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements

void SparseMatrixMatrixMultiplyOperation::SetTemporaryArrays (double* tmpMatrixElements, int* tmpColumnIndices)
{
  this->LocalTmpMatrixElements = tmpMatrixElements;
  this->LocalTmpColumnIndices = tmpColumnIndices;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* SparseMatrixMatrixMultiplyOperation::Clone()
{
  return new SparseMatrixMatrixMultiplyOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool SparseMatrixMatrixMultiplyOperation::RawApplyOperation()
{
  long TmpNbrMatrixElements = 0l;
  long PreviousTmpNbrMatrixElements = 0l;
  double* TmpElements = new double [this->RightMatrix->NbrColumn];
  for (int i = 0; i < this->RightMatrix->NbrColumn; ++i)
    {
      TmpElements[i] = 0.0;
    }
  int LastComponent = this->FirstComponent + this->NbrComponent;
  for (int i = this->FirstComponent; i < this->LeftMatrix->NbrRow; ++i)
    {
      long MinPos =  this->LeftMatrix->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->LeftMatrix->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = this->LeftMatrix->ColumnIndices[MinPos];
	      long MinPos2 = this->RightMatrix->RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = this->LeftMatrix->MatrixElements[MinPos];
		  long MaxPos2 = this->RightMatrix->RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      TmpElements[this->RightMatrix->ColumnIndices[MinPos2]] += Tmp * this->RightMatrix->MatrixElements[MinPos2];
		    }      
		}
	    }	 
   
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < this->RightMatrix->NbrColumn; ++j)
	    if (TmpElements[j] != 0.0)
	      {
		this->LocalTmpMatrixElements[TmpNbrMatrixElements] = TmpElements[j];
		this->LocalTmpColumnIndices[TmpNbrMatrixElements] = j;
		TmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	  if (TmpNbrMatrixElements == PreviousTmpNbrMatrixElements)
	    {
	      this->DestinationMatrix.RowPointers[i] = -1l;
	      this->DestinationMatrix.RowLastPointers[i] = 1l;
	    }
	  else
	    {
	      this->DestinationMatrix.RowPointers[i] = PreviousTmpNbrMatrixElements;
	      this->DestinationMatrix.RowLastPointers[i] = TmpNbrMatrixElements - 1;
	    }
	}
      else
	{
	  this->DestinationMatrix.RowPointers[i] = -1l;
	  this->DestinationMatrix.RowLastPointers[i] = -1;
	}
    }
  delete[] TmpElements;
  this->LocalNbrMatrixElements = TmpNbrMatrixElements;
  return true;
}

// apply operation for mono processor architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool SparseMatrixMatrixMultiplyOperation::ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture)
{
  this->FirstComponent = 0;
  this->NbrComponent = this->LeftMatrix->GetNbrRow();
  this->RawApplyOperation();

  this->DestinationMatrix.NbrMatrixElements = this->LocalNbrMatrixElements; 
  this->DestinationMatrix.MatrixElements = new double[this->DestinationMatrix.NbrMatrixElements];
  this->DestinationMatrix.ColumnIndices = new int[this->DestinationMatrix.NbrMatrixElements];
  long TmpNbrMatrixElements = this->DestinationMatrix.NbrMatrixElements;
  for (long j = 0l; j < TmpNbrMatrixElements; ++j)
    {
      this->DestinationMatrix.MatrixElements[j] = this->LocalTmpMatrixElements[j];
      this->DestinationMatrix.ColumnIndices[j] = this->LocalTmpColumnIndices[j];
    }
  return true;
}
  
// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool SparseMatrixMatrixMultiplyOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->LeftMatrix->GetNbrRow() / architecture->GetNbrThreads();
  int TmpFirstComponent = 0;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  SparseMatrixMatrixMultiplyOperation** TmpOperations = new SparseMatrixMatrixMultiplyOperation* [architecture->GetNbrThreads()];
  long StepTemporaryArrays = this->MaxTmpMatrixElements / architecture->GetNbrThreads();
  long TmpNbrMatrixElements = 0l;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (SparseMatrixMatrixMultiplyOperation*) this->Clone();
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      TmpOperations[i]->SetTemporaryArrays(this->LocalTmpMatrixElements + TmpNbrMatrixElements, this->LocalTmpColumnIndices + TmpNbrMatrixElements);
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpFirstComponent += Step;
      TmpNbrMatrixElements += StepTemporaryArrays;
    }
  TmpOperations[ReducedNbrThreads] = (SparseMatrixMatrixMultiplyOperation*) this->Clone();
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->LeftMatrix->GetNbrRow() - TmpFirstComponent);  
  TmpOperations[ReducedNbrThreads]->SetTemporaryArrays(this->LocalTmpMatrixElements + TmpNbrMatrixElements, this->LocalTmpColumnIndices + TmpNbrMatrixElements);
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);

  architecture->SendJobs();

  this->DestinationMatrix.NbrMatrixElements = 0l; 
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      this->DestinationMatrix.NbrMatrixElements += TmpOperations[i]->LocalNbrMatrixElements;
    }
  this->DestinationMatrix.MatrixElements = new double[this->DestinationMatrix.NbrMatrixElements];
  this->DestinationMatrix.ColumnIndices = new int[this->DestinationMatrix.NbrMatrixElements];
  long TmpIndex = 0l;
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      double* TmpMatrixElements = TmpOperations[i]->LocalTmpMatrixElements;
      int* TmpColumnIndices = TmpOperations[i]->LocalTmpColumnIndices;
      long TmpNbrMatrixElements = TmpOperations[i]->LocalNbrMatrixElements;
      for (long j = 0l; j < TmpNbrMatrixElements; ++j)
	{
	  this->DestinationMatrix.MatrixElements[TmpIndex] = TmpMatrixElements[j];
	  this->DestinationMatrix.ColumnIndices[TmpIndex] = TmpColumnIndices[j];
	  ++TmpIndex;
	}
    }
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}

