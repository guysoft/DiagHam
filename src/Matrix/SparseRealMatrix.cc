////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of real matrix with sparse storage                  //
//                                                                            //
//                        last modification : 17/10/2012                      //
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


#include "Matrix/SparseRealMatrix.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <cstdlib>


using std::endl;
using std::cout;




// default constructor
//

SparseRealMatrix::SparseRealMatrix() 
{
  this->MatrixElements = 0;
  this->ColumnIndices = 0;
  this->RowPointers = 0;
  this->RowLastPointers = 0;
  this->NbrMatrixElements = 0l;
  this->MaximumNbrMatrixElements = 0l;
  this->NbrMatrixElementPacketSize = 0l;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;  
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
}

// constructor for a sparse matrix without any specific struture
//
// nbrRow = number of rows
// nbrColumn = number of columns

SparseRealMatrix::SparseRealMatrix(int nbrRow, int nbrColumn)
{
  this->Flag.Initialize();
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->NbrMatrixElements = 0;
  this->RowPointers = new long[this->NbrRow];
  this->RowLastPointers = new long[this->NbrRow];
  this->NbrMatrixElementPacketSize = 1024l;
  this->MaximumNbrMatrixElements = this->NbrMatrixElementPacketSize;
  this->MatrixElements = new double[this->MaximumNbrMatrixElements];
  this->ColumnIndices = new int[this->MaximumNbrMatrixElements];
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->RowPointers[i] = -1l;
      this->RowLastPointers[i] = -1l;
    }
}

// constructor for a sparse matrix without any specific struture but a given number of non-zero matrix elements
//
// nbrRow = number of rows
// nbrColumn = number of columns
// nbrMatrixElements = number of non-zero matrix elements
// zero = true if matrix elements have to be set to zero

SparseRealMatrix::SparseRealMatrix(int nbrRow, int nbrColumn, long nbrMatrixElements, bool zero)
{
  this->Flag.Initialize();
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->NbrMatrixElements = nbrMatrixElements;
  this->MaximumNbrMatrixElements = nbrMatrixElements;
  this->NbrMatrixElementPacketSize = 0l;
  this->RowPointers = new long[this->NbrRow];
  this->RowLastPointers = new long[this->NbrRow];
  this->MatrixElements = new double[this->NbrMatrixElements];
  this->ColumnIndices = new int[this->NbrMatrixElements];
  if (zero == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	{
	  this->RowPointers[i] = -1l;
	  this->RowLastPointers[i] = -1l;
	}
      for (long i = 0l; i < this->NbrMatrixElements; ++i)
	{
	  this->MatrixElements[i] = 0.0;
	  this->ColumnIndices[i] = -1;
	}
   }
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

SparseRealMatrix::SparseRealMatrix(const SparseRealMatrix& M) 
{
  this->Flag = M.Flag;
  this->MatrixElements = M.MatrixElements;
  this->ColumnIndices = M.ColumnIndices;
  this->RowPointers = M.RowPointers;
  this->RowLastPointers = M.RowLastPointers;
  this->NbrMatrixElements = M.NbrMatrixElements;
  this->MaximumNbrMatrixElements = M.MaximumNbrMatrixElements;
  this->NbrMatrixElementPacketSize = M.NbrMatrixElementPacketSize;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;  
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
}

// copy constructor (duplicating all datas)
//
// M = matrix to copy

SparseRealMatrix::SparseRealMatrix(Matrix& M)
{
  if ((M.GetNbrRow() == 0) || (M.GetNbrColumn() == 0))
    {
      this->MatrixElements = 0;
      this->ColumnIndices = 0;
      this->RowPointers = 0;
      this->RowLastPointers = 0;
      this->NbrMatrixElements = 0l;
      this->MaximumNbrMatrixElements = 0l;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::RealElements | Matrix::Sparse;
    }
  else
    {
      this->Flag.Initialize();
      this->NbrColumn = M.GetNbrColumn();
      this->NbrRow = M.GetNbrRow();
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      this->NbrMatrixElements = 0l;
      double Tmp;
      this->RowPointers = new long[this->NbrRow];
      this->RowLastPointers = new long[this->NbrRow];
      for (int i = 0; i < this->NbrRow; i++)
	{
	  long PreviousNbrMatrixElements = this->NbrMatrixElements;
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      M.GetMatrixElement(i, j, Tmp);
	      if (Tmp != 0.0)
		++this->NbrMatrixElements;
	    }
	  if (PreviousNbrMatrixElements == this->NbrMatrixElements)
	    {
	      this->RowPointers[i] = -1l;
	      this->RowLastPointers[i] = -1l;
	    }
	  else
	    {
	      this->RowPointers[i] = PreviousNbrMatrixElements;	      
	      this->RowLastPointers[i] = this->NbrMatrixElements - 1l;
	    }
	}
      this->MatrixElements = new double[this->NbrMatrixElements];
      this->ColumnIndices = new int[this->NbrMatrixElements];
      this->NbrMatrixElements = 0l;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      M.GetMatrixElement(i, j, Tmp);
	      if (Tmp != 0.0)
		{
		  this->MatrixElements[this->NbrMatrixElements] = Tmp;
		  this->ColumnIndices[this->NbrMatrixElements] = j;
		  ++this->NbrMatrixElements;
		}
	    }
	}
      this->MaximumNbrMatrixElements = this->NbrMatrixElements;
      this->NbrMatrixElementPacketSize = 0l;
      this->MatrixType = Matrix::RealElements | Matrix::Sparse;
    }
}

// destructor
//

SparseRealMatrix::~SparseRealMatrix() 
{
  if ((this->MatrixElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->MatrixElements;
	delete[] this->ColumnIndices;
	delete[] this->RowPointers;
	delete[] this->RowLastPointers;
      }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

SparseRealMatrix& SparseRealMatrix::operator = (const SparseRealMatrix& M) 
{
  if ((this->MatrixElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->MatrixElements;
	delete[] this->ColumnIndices;
	delete[] this->RowPointers;
	delete[] this->RowLastPointers;
      }
  this->Flag = M.Flag;
  this->MatrixElements = M.MatrixElements;
  this->ColumnIndices = M.ColumnIndices;
  this->RowPointers = M.RowPointers;
  this->RowLastPointers = M.RowLastPointers;
  this->NbrMatrixElements = M.NbrMatrixElements;
  this->MaximumNbrMatrixElements = M.MaximumNbrMatrixElements;
  this->NbrMatrixElementPacketSize = M.NbrMatrixElementPacketSize;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;  
  this->MatrixType = Matrix::RealElements | Matrix::Sparse;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* SparseRealMatrix::Clone ()
{
  return ((Matrix*) new SparseRealMatrix (*this));
}

// copy a matrix into another (duplicating data)
//
// matrix = matrix to copy
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::Copy (SparseRealMatrix& matrix)
{
  if ((this->MatrixElements != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->MatrixElements;
	delete[] this->ColumnIndices;
	delete[] this->RowPointers;
	delete[] this->RowLastPointers;
      }
  this->Flag.Initialize();
  this->NbrMatrixElements = matrix.NbrMatrixElements;
  this->NbrRow = matrix.NbrRow;
  this->NbrColumn = matrix.NbrColumn;
  this->TrueNbrRow = matrix.TrueNbrRow;
  this->TrueNbrColumn = matrix.TrueNbrColumn;  
  this->RowPointers = new long[this->NbrRow];
  this->RowLastPointers = new long[this->NbrRow];
  this->MatrixElements = new double[this->NbrMatrixElements];
  this->ColumnIndices = new int[this->NbrMatrixElements];
  for (int j = 0; j < this->NbrRow; ++j)
    {
      this->RowPointers[j] = matrix.RowPointers[j];
      this->RowLastPointers[j] = matrix.RowLastPointers[j];
    }
  for (long j = 0l; j < this->NbrMatrixElements; ++j)
    {
      this->MatrixElements[j] = matrix.MatrixElements[j];
      this->ColumnIndices[j] = matrix.ColumnIndices[j];
    }
  return *this;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void SparseRealMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  if (this->NbrMatrixElementPacketSize == 0l)
    {
      if (this->RowPointers[i] >= 0l)
	{
	  long TmpIndex = this->FindColumnIndexPosition(j, this->RowPointers[i], this->RowLastPointers[i]);
	  if (TmpIndex >= 0l)
	    {
	      this->MatrixElements[TmpIndex] = x;
	      return;
	    }
	}     
      return;
    }
  if (this->RowPointers[i] >= 0l)
    {
      long TmpIndex = this->FindColumnIndexPosition(j, this->RowPointers[i], this->RowLastPointers[i]);
      if (TmpIndex >= 0l)
	{
	  this->MatrixElements[TmpIndex] = x;
	  return;
	}
      this->IncreaseNbrMatrixElements();
      for (int k = i + 1; k < this->NbrRow; ++k)
	{
	  if (this->RowPointers[k] >= 0l)
	    {
	      ++this->RowPointers[k];
	      ++this->RowLastPointers[k];
	    }
	}
      TmpIndex = this->RowPointers[i];
      while ((TmpIndex <= this->RowLastPointers[i]) && (this->ColumnIndices[TmpIndex] < j))
	++TmpIndex;
      for (long k = this->NbrMatrixElements - 1; k > TmpIndex;  --k)
	{
	  this->MatrixElements[k] = this->MatrixElements[k - 1];
	  this->ColumnIndices[k] = this->ColumnIndices[k - 1];	  
	}
      this->MatrixElements[TmpIndex] = x;
      this->ColumnIndices[TmpIndex] = j;
      ++this->RowLastPointers[i];
      return;
    }
  int TmpIndex1 = i;
  while ((TmpIndex1 >= 0) && (this->RowPointers[TmpIndex1] == -1l))
    --TmpIndex1;
  long TmpIndex = 0l;
  if (TmpIndex1 >= 0)
    TmpIndex = this->RowLastPointers[TmpIndex1] + 1l;

  this->IncreaseNbrMatrixElements();
  for (int k = i + 1; k < this->NbrRow; ++k)
    {
      if (this->RowPointers[k] >= 0l)
	{
	  ++this->RowPointers[k];
	  ++this->RowLastPointers[k];
	}
    }
  for (long k = this->NbrMatrixElements - 1; k > TmpIndex;  --k)
    {
      this->MatrixElements[k] = this->MatrixElements[k - 1];
      this->ColumnIndices[k] = this->ColumnIndices[k - 1];	  
    }
  this->MatrixElements[TmpIndex] = x;
  this->ColumnIndices[TmpIndex] = j;
  this->RowPointers[i] = TmpIndex;
  this->RowLastPointers[i] = TmpIndex;
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void SparseRealMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn) || (this->RowPointers[i] == -1l))
    return;
  long TmpIndex = this->FindColumnIndexPosition(j, this->RowPointers[i], this->RowLastPointers[i]);
  if (TmpIndex == -1l)
    {
      return;
    }
  this->MatrixElements[TmpIndex] += x;
}



// increase the number of matrix elements
//
// nbrElements = number of elements to add

void SparseRealMatrix::IncreaseNbrMatrixElements(long nbrElements)
{
  long TmpNbrMatrixElements = this->NbrMatrixElements + nbrElements;
  if (TmpNbrMatrixElements <= this->MaximumNbrMatrixElements)
    {
      this->NbrMatrixElements = TmpNbrMatrixElements;
      return;
    }
  this->MaximumNbrMatrixElements += this->NbrMatrixElementPacketSize;
  double* TmpMatrixElements = new double[this->MaximumNbrMatrixElements];
  int* TmpColumnIndices = new int[this->MaximumNbrMatrixElements];
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    {
      TmpMatrixElements[i] = this->MatrixElements[i];
      TmpColumnIndices[i] = this->ColumnIndices[i];
    }
  this->NbrMatrixElements = TmpNbrMatrixElements;
  if (this->Flag.Shared() == false)
    {
      delete[] this->MatrixElements;
      delete[] this->ColumnIndices;
    }
  else
    {
      this->Flag.Initialize();
    }
  this->MatrixElements = TmpMatrixElements;
  this->ColumnIndices = TmpColumnIndices;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void SparseRealMatrix::Resize (int nbrRow, int nbrColumn)
{
//   if (this->NbrRow != nbrRow)
//     {
//       for (int i = 0; i < this->NbrColumn; i++)
// 	this->Columns[i].Resize(nbrRow);
//       if (this->TrueNbrRow >= nbrRow)
// 	{
// 	  this->NbrRow = nbrRow;
// 	}
//       else
// 	{
// 	  this->NbrRow = nbrRow;
// 	  this->TrueNbrRow = nbrRow;
// 	}
//     }
//   if (this->TrueNbrColumn >= nbrColumn)
//     {
//       for (int i = this->NbrColumn; i < nbrColumn; i++)
// 	this->Columns[i].Resize(nbrRow);
//       this->NbrColumn = nbrColumn;
//     }
//   else
//     {
//       RealVector* Tmp = new RealVector[nbrColumn];
//       for (int i = 0; i < this->NbrColumn; i++)
// 	Tmp[i] = this->Columns[i];      
//       for (int i = this->NbrColumn; i < nbrColumn; i++)
// 	Tmp[i] = RealVector(nbrRow);
//       if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
// 	{
// 	  delete[] this->Columns;
// 	}
//       this->Columns = Tmp;
//       this->TrueNbrColumn = nbrColumn;
//       this->NbrColumn = nbrColumn;
//       this->Flag = GarbageFlag();
//       this->Flag.Initialize();
//     }
//   return;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void SparseRealMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
//   if (this->NbrRow != nbrRow)
//     {
//       for (int i = 0; i < this->NbrColumn; i++)
// 	this->Columns[i].ResizeAndClean(nbrRow);
//       if (this->TrueNbrRow >= nbrRow)
// 	{
// 	  this->NbrRow = nbrRow;
// 	}
//       else
// 	{
// 	  this->NbrRow = nbrRow;
// 	  this->TrueNbrRow = nbrRow;
// 	}
//     }
//   if (this->TrueNbrColumn >= nbrColumn)
//     {
//       for (int i = this->NbrColumn; i < nbrColumn; i++)
// 	this->Columns[i].ResizeAndClean(nbrRow);
//       this->TrueNbrColumn = nbrColumn;
//     }
//   else
//     {
//       RealVector* Tmp = new RealVector[nbrColumn];
//       for (int i = 0; i < this->NbrColumn; i++)
// 	Tmp[i] = this->Columns[i];      
//       for (int i = this->NbrColumn; i < nbrColumn; i++)
// 	Tmp[i] = RealVector(nbrRow, true);
//       if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
// 	{
// 	  delete[] this->Columns;
// 	}
//       this->Columns = Tmp;
//       this->Flag = GarbageFlag();
//       this->Flag.Initialize();
//       this->TrueNbrColumn = nbrColumn;
//       this->NbrColumn = nbrColumn;
//     }
//   return;
}

// Set all entries in matrix to zero
//

void SparseRealMatrix::ClearMatrix ()
{
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    this->MatrixElements[i] = 0.0;
  return;
}

// add two matrices
//
// matrix1 = first matrix
// matrix2 = second matrix
// return value = sum of the two matrices

SparseRealMatrix operator + (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2)
{
  return SparseRealMatrixLinearCombination(1.0, matrix1, 1.0, matrix2);
}

// difference of two matrices
//
// matrix1 = first matrix
// matrix2 = second matrix
// return value = difference of the two matrices

SparseRealMatrix operator - (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2)
{
  return SparseRealMatrixLinearCombination(1.0, matrix1, -1.0, matrix2);
}

// create the linear combination of two matrices
//
// x1 = prefactor of the first matrix
// matrix1 = first matrix
// x2 = prefactor of the second matrix
// matrix2 = second matrix
// return value = linear combination

SparseRealMatrix SparseRealMatrixLinearCombination(const double& x1, const SparseRealMatrix& matrix1, const double& x2, const SparseRealMatrix& matrix2)
{
  if ((matrix1.NbrRow != matrix2.NbrRow) || (matrix1.NbrColumn != matrix2.NbrColumn))
    return SparseRealMatrix();
  long TmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      if (matrix1.RowPointers[i] < 0l)
	{
	  if (matrix2.RowPointers[i] >= 0l)
	    TmpNbrMatrixElements += matrix2.RowLastPointers[i] - matrix2.RowPointers[i] + 1l;	 	  
	}
      else
	{
	  if (matrix2.RowPointers[i] < 0l)
	    TmpNbrMatrixElements += matrix1.RowLastPointers[i] - matrix1.RowPointers[i] + 1l;	 	  
	  else
	    {
	      long MinPos1 = matrix1.RowPointers[i];
	      long MinPos2 = matrix2.RowPointers[i];
	      long MaxPos1 = matrix1.RowLastPointers[i];
	      long MaxPos2 = matrix2.RowLastPointers[i];
	      while (MinPos1 <= MaxPos1)
		{
		  while ((MinPos1 <= MaxPos1) && (matrix1.ColumnIndices[MinPos1] < matrix2.ColumnIndices[MinPos2]))
		    {
		      ++TmpNbrMatrixElements;
		      ++MinPos1;
		    }
		  if (MinPos1 <= MaxPos1)
		    {
		      if (matrix1.ColumnIndices[MinPos1] == matrix2.ColumnIndices[MinPos2])
			{
			  ++TmpNbrMatrixElements;
			  ++MinPos1;
			  ++MinPos2;
			}
		      while ((MinPos2 <= MaxPos2) && (matrix2.ColumnIndices[MinPos2] < matrix1.ColumnIndices[MinPos1]))
			{
			  ++TmpNbrMatrixElements;
			  ++MinPos2;
			}
		      if ((MinPos2 > MaxPos2) && (MinPos1 <= MaxPos1))
			{
			  TmpNbrMatrixElements += MaxPos1 - MinPos1 + 1l;
			  MinPos1 = MaxPos1 + 1l;
			}
		    }
		}
	      if (MinPos2 <= MaxPos2)
		{
		  TmpNbrMatrixElements += MaxPos2 - MinPos2 + 1l;
		}		  
	    }
	}
    }  

  SparseRealMatrix TmpMatrix (matrix1.NbrRow, matrix1.NbrColumn, TmpNbrMatrixElements);
  TmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      if (matrix1.RowPointers[i] < 0l)
	{
	  if (matrix2.RowPointers[i] >= 0l)
	    {
	      TmpMatrix.RowPointers[i] = TmpNbrMatrixElements;
	      long MinPos2 = matrix2.RowPointers[i];
	      long MaxPos2 = matrix2.RowLastPointers[i];
	      for (; MinPos2 <= MaxPos2; ++MinPos2)
		{
		  TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x2 * matrix2.MatrixElements[MinPos2];
		  TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix2.ColumnIndices[MinPos2];
		  ++TmpNbrMatrixElements;
		}
	      TmpMatrix.RowLastPointers[i] = TmpNbrMatrixElements - 1l;
	    }
	  else
	    {
	      TmpMatrix.RowLastPointers[i] = -1l;
	      TmpMatrix.RowPointers[i] = -1l;
	    }
	}
      else
	{
	  if (matrix2.RowPointers[i] < 0l)
	    {
	      TmpMatrix.RowPointers[i] = TmpNbrMatrixElements;
	      long MinPos1 = matrix1.RowPointers[i];
	      long MaxPos1 = matrix1.RowLastPointers[i];
	      for (; MinPos1 <= MaxPos1; ++MinPos1)
		{
		  TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x1 * matrix1.MatrixElements[MinPos1];
		  TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix1.ColumnIndices[MinPos1];
		  ++TmpNbrMatrixElements;
		}
	      TmpMatrix.RowLastPointers[i] = TmpNbrMatrixElements - 1l;
	    }
	  else
	    {
	      long MinPos1 = matrix1.RowPointers[i];
	      long MinPos2 = matrix2.RowPointers[i];
	      long MaxPos1 = matrix1.RowLastPointers[i];
	      long MaxPos2 = matrix2.RowLastPointers[i];
	      TmpMatrix.RowPointers[i] = TmpNbrMatrixElements;
	      while (MinPos1 <= MaxPos1)
		{
		  while ((MinPos1 <= MaxPos1) && (matrix1.ColumnIndices[MinPos1] < matrix2.ColumnIndices[MinPos2]))
		    {
		      TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x1 * matrix1.MatrixElements[MinPos1];
		      TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix1.ColumnIndices[MinPos1];
		      ++TmpNbrMatrixElements;
		      ++MinPos1;
		    }
		  if (MinPos1 <= MaxPos1)
		    {
		      if (matrix1.ColumnIndices[MinPos1] == matrix2.ColumnIndices[MinPos2])
			{
			  TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x1 * matrix1.MatrixElements[MinPos1] + x2 * matrix2.MatrixElements[MinPos2];
			  TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix1.ColumnIndices[MinPos1];
			  ++TmpNbrMatrixElements;
			  ++MinPos1;
			  ++MinPos2;
			}
		      while ((MinPos2 <= MaxPos2) && (matrix2.ColumnIndices[MinPos2] < matrix1.ColumnIndices[MinPos1]))
			{
			  TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x2 * matrix2.MatrixElements[MinPos2];
			  TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix2.ColumnIndices[MinPos2];
			  ++TmpNbrMatrixElements;
			  ++MinPos2;
			}		  
		      if (MinPos2 > MaxPos2)		  
			{
			  while (MinPos1 <= MaxPos1)
			    {
			      TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x1 * matrix1.MatrixElements[MinPos1];
			      TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix1.ColumnIndices[MinPos1];
			      ++TmpNbrMatrixElements;
			      ++MinPos1;
			    }
			}
		    }
		}
	      while (MinPos2 <= MaxPos2)
		{
		  TmpMatrix.MatrixElements[TmpNbrMatrixElements] = x2 * matrix2.MatrixElements[MinPos2];
		  TmpMatrix.ColumnIndices[TmpNbrMatrixElements] = matrix2.ColumnIndices[MinPos2];
		  ++TmpNbrMatrixElements;
		  ++MinPos2;
		}		  
	      TmpMatrix.RowLastPointers[i] = TmpNbrMatrixElements - 1l;
	    }
	}
    }      
  return TmpMatrix;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::operator *= (double x) 
{
  for (long i = 0; i < this->NbrMatrixElements; ++i)
    this->MatrixElements[i] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::operator /= (double x)
{
  x = 1.0 / x;;
  for (long i = 0; i < this->NbrMatrixElements; ++i)
    this->MatrixElements[i] *= x;
  return *this;
}

// multiply a matrix to the right by another matrix
//
// matrix = matrix used as multiplicator
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::Multiply (const SparseRealMatrix& matrix)
{
  if (matrix.NbrRow != this->NbrColumn)
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return *this; 
    }
  double* TmpMatrixElements = new double[this->NbrMatrixElements * matrix.NbrMatrixElements];
  int* TmpColumnIndices = new int[this->NbrMatrixElements * matrix.NbrMatrixElements];
  long TmpNbrMatrixElements = 0l;
  long PreviousTmpNbrMatrixElements = 0l;
  double* TmpElements = new double [matrix.NbrColumn];
  for (int i = 0; i < matrix.NbrColumn; ++i)
    {
      TmpElements[i] = 0.0;
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = this->ColumnIndices[MinPos];
	      long MinPos2 = matrix.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = this->MatrixElements[MinPos];
		  long MaxPos2 = matrix.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      TmpElements[matrix.ColumnIndices[MinPos2]] += Tmp * matrix.MatrixElements[MinPos2];
		    }      
		}
	    }	 
   
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < matrix.NbrColumn; ++j)
	    if (TmpElements[j] != 0)
	      {
		TmpMatrixElements[TmpNbrMatrixElements] = TmpElements[j];
		TmpColumnIndices[TmpNbrMatrixElements] = j;
		TmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	  this->RowPointers[i] = PreviousTmpNbrMatrixElements;
	  this->RowLastPointers[i] = TmpNbrMatrixElements - 1;
	}
    }
  delete[] TmpElements;
  delete[] this->MatrixElements;
  delete[] this->ColumnIndices;
  this->NbrMatrixElements = TmpNbrMatrixElements;
  this->MatrixElements = new double[this->NbrMatrixElements];
  this->ColumnIndices = new int[this->NbrMatrixElements];
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    {
      this->MatrixElements[i] = TmpMatrixElements[i];
      this->ColumnIndices[i] = TmpColumnIndices[i];
    }
  delete[] TmpMatrixElements;
  delete[] TmpColumnIndices;
  return *this;
}

// multiply a matrix to the right by another matrix, providing all the required temporary arrays
//
// matrix = matrix used as multiplicator
// tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::Multiply (const SparseRealMatrix& matrix, 
					      double* tmpMatrixElements, int* tmpColumnIndices, double* tmpElements)
{
  if (matrix.NbrRow != this->NbrColumn)
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return *this; 
    }
  long TmpNbrMatrixElements = 0l;
  long PreviousTmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix.NbrColumn; ++i)
    {
      tmpElements[i] = 0.0;
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = this->ColumnIndices[MinPos];
	      long MinPos2 = matrix.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = this->MatrixElements[MinPos];
		  long MaxPos2 = matrix.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      tmpElements[matrix.ColumnIndices[MinPos2]] += Tmp * matrix.MatrixElements[MinPos2];
		    }      
		}
	    }	 
   
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < matrix.NbrColumn; ++j)
	    if (tmpElements[j] != 0)
	      {
		tmpMatrixElements[TmpNbrMatrixElements] = tmpElements[j];
		tmpColumnIndices[TmpNbrMatrixElements] = j;
		tmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	  this->RowPointers[i] = PreviousTmpNbrMatrixElements;
	  this->RowLastPointers[i] = TmpNbrMatrixElements - 1;
	}
    }
  delete[] this->MatrixElements;
  delete[] this->ColumnIndices;
  this->NbrMatrixElements = TmpNbrMatrixElements;
  this->MatrixElements = new double[this->NbrMatrixElements];
  this->ColumnIndices = new int[this->NbrMatrixElements];
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    {
      this->MatrixElements[i] = tmpMatrixElements[i];
      this->ColumnIndices[i] = tmpColumnIndices[i];
    }
  return *this;
}


// multiply a matrix to the right by another matrix, providing all the required temporary arrays, extend their capacity if needed
//
// matrix = matrix used as multiplicator
// tmpMatrixElements = reference on the temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
// tmpColumnIndices = reference on the temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
// nbrElements = reference ont the number of elements in tmpMatrixElements and tmpColumnIndices
// tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
// return value = reference on current matrix

SparseRealMatrix& SparseRealMatrix::Multiply (const SparseRealMatrix& matrix, 
					      double*& tmpMatrixElements, int*& tmpColumnIndices, 
					      long& nbrElements, double* tmpElements)
{
  if (matrix.NbrRow != this->NbrColumn)
    {
      cout << "error, cannot multiply the two matrices" << endl;
      return *this; 
    }
  long TmpNbrMatrixElements = 0l;
  long PreviousTmpNbrMatrixElements = 0l;
  for (int i = 0; i < matrix.NbrColumn; ++i)
    {
      tmpElements[i] = 0.0;
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      int TmpIndex = this->ColumnIndices[MinPos];
	      long MinPos2 = matrix.RowPointers[TmpIndex];
	      if (MinPos2 >= 0)
		{
		  double Tmp = this->MatrixElements[MinPos];
		  long MaxPos2 = matrix.RowLastPointers[TmpIndex];
		  for (; MinPos2 <= MaxPos2; ++MinPos2)
		    {
		      tmpElements[matrix.ColumnIndices[MinPos2]] += Tmp * matrix.MatrixElements[MinPos2];
		    }      
		}
	    }	 
   
	  PreviousTmpNbrMatrixElements = TmpNbrMatrixElements;
	  for (int j = 0; j < matrix.NbrColumn; ++j)
	    if (tmpElements[j] != 0)
	      {
		tmpMatrixElements[TmpNbrMatrixElements] = tmpElements[j];
		tmpColumnIndices[TmpNbrMatrixElements] = j;
		tmpElements[j] = 0.0;
		++TmpNbrMatrixElements;
	      }	  
	  this->RowPointers[i] = PreviousTmpNbrMatrixElements;
	  this->RowLastPointers[i] = TmpNbrMatrixElements - 1;
	}
    }
  delete[] this->MatrixElements;
  delete[] this->ColumnIndices;
  this->NbrMatrixElements = TmpNbrMatrixElements;
  this->MatrixElements = new double[this->NbrMatrixElements];
  this->ColumnIndices = new int[this->NbrMatrixElements];
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    {
      this->MatrixElements[i] = tmpMatrixElements[i];
      this->ColumnIndices[i] = tmpColumnIndices[i];
    }
  return *this;
}

// compute the tensor product of two sparse matrices (matrix1 x matrix2), and store the result in a sparse matrix
//
// matrix1 = reference on the left matrix
// matrix2 = reference on the right matrix
// return value = tensor product

SparseRealMatrix TensorProduct (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2)
{
  SparseRealMatrix TmpMatrix (matrix1.NbrRow * matrix2.NbrRow, 
				 matrix1.NbrColumn * matrix2.NbrColumn, 
				 matrix1.NbrMatrixElements * matrix2.NbrMatrixElements);
  for (int i = 0; i < TmpMatrix.NbrRow; ++i)
    {
      TmpMatrix.RowPointers[i] = -1;
      TmpMatrix.RowLastPointers[i] = -1;
    }
  long TmpPosition = 0l;
  for (int i = 0; i < matrix1.NbrRow; ++i)
    {
      if (matrix1.RowPointers[i] >= 0l)
	{
	  long MinM1 = matrix1.RowPointers[i];
	  long MaxM1 = matrix1.RowLastPointers[i];
	  for (int j = 0; j < matrix2.NbrRow; ++j)
	    {
	      if (matrix2.RowPointers[j] >= 0l)
		{
		  int RowIndex = i *  matrix2.NbrRow + j;
		  TmpMatrix.RowPointers[RowIndex] = TmpPosition;
		  long MinM2 = matrix2.RowPointers[j];
		  long MaxM2 = matrix2.RowLastPointers[j];
		  for (long k1 = MinM1; k1 <= MaxM1; ++k1)
		    {
		      int Shift =  matrix1.ColumnIndices[k1] * matrix2.NbrColumn;
		      double Tmp =  matrix1.MatrixElements[k1];;
		      for (long k2 = MinM2; k2 <= MaxM2; ++k2)
			{
			  TmpMatrix.MatrixElements[TmpPosition] =  Tmp * matrix2.MatrixElements[k2];
			  TmpMatrix.ColumnIndices[TmpPosition] =  matrix2.ColumnIndices[k2] + Shift;			
			  ++TmpPosition;
			}
		    }
		  TmpMatrix.RowLastPointers[RowIndex] = TmpPosition - 1l;
		}
	    }
	}
    }
  return TmpMatrix;  
}

// compute the transpose of the current matrix
//
// return value = hermitian transposed matrix

SparseRealMatrix SparseRealMatrix::Transpose ()
{
  SparseRealMatrix TmpMatrix (this->NbrColumn, this->NbrRow, this->NbrMatrixElements);
  for (int i = 0; i < TmpMatrix.NbrRow; ++i)
    {
      TmpMatrix.RowPointers[i] = 0l;
    }
  for (long i = 0l; i < this->NbrMatrixElements; ++i)
    ++TmpMatrix.RowPointers[this->ColumnIndices[i]];
  int PreviousIndex = 0;
  while ((PreviousIndex < TmpMatrix.NbrRow) && (TmpMatrix.RowPointers[PreviousIndex] == 0l))
    {
      TmpMatrix.RowPointers[PreviousIndex] = -1l;
      TmpMatrix.RowLastPointers[PreviousIndex] = -1l;
      ++PreviousIndex;        
    }
  int Index = PreviousIndex + 1;
  TmpMatrix.RowLastPointers[PreviousIndex] = TmpMatrix.RowPointers[PreviousIndex] - 1l;
  TmpMatrix.RowPointers[PreviousIndex] = 0l;
  while (Index < TmpMatrix.NbrRow)
    {
      while ((Index < TmpMatrix.NbrRow) && (TmpMatrix.RowPointers[Index] == 0l))
	{
	  TmpMatrix.RowPointers[Index] = -1l;
	  TmpMatrix.RowLastPointers[Index] = -1l;
	  ++Index;
	}
      if (Index < TmpMatrix.NbrRow)
	{
	  TmpMatrix.RowLastPointers[Index] = TmpMatrix.RowPointers[Index] + TmpMatrix.RowLastPointers[PreviousIndex];;
	  TmpMatrix.RowPointers[Index] = TmpMatrix.RowLastPointers[PreviousIndex] + 1l;
	  PreviousIndex = Index;	  
	  ++Index;
	}
    }
  for (int i = 0; i < this->NbrRow; ++i)
    {
      long MinPos =  this->RowPointers[i];
      if (MinPos >= 0l)
	{
	  long MaxPos = this->RowLastPointers[i];
	  for (; MinPos <= MaxPos; ++MinPos)
	    {
	      long& TmpPointer = TmpMatrix.RowPointers[this->ColumnIndices[MinPos]];
	      TmpMatrix.MatrixElements[TmpPointer] = this->MatrixElements[MinPos];
	      TmpMatrix.ColumnIndices[TmpPointer] = i;
	      ++TmpPointer;
	    }
	}
    }
  Index = 0;
  while ((Index < TmpMatrix.NbrRow) && (TmpMatrix.RowPointers[Index] < 0l))
    {
      ++Index;        
    }
  long PreviousPointer =  TmpMatrix.RowPointers[Index];
  TmpMatrix.RowPointers[Index] = 0l;
  ++Index;
   while (Index < TmpMatrix.NbrRow)
    {
      while ((Index < TmpMatrix.NbrRow) && (TmpMatrix.RowPointers[Index] < 0l))
	{
	  ++Index;
	}
      if (Index < TmpMatrix.NbrRow)
	{
	  long Tmp = TmpMatrix.RowPointers[Index];
	  TmpMatrix.RowPointers[Index] = PreviousPointer;
	  PreviousPointer = Tmp;
	  ++Index;
	}
    }
	
  return TmpMatrix;
}

// compute the number of non-zero matrix elements (zero having strictly zero square norm)
//
// return value = number of non-zero matrix elements

long SparseRealMatrix::ComputeNbrNonZeroMatrixElements()
{
  long NbrNonZero = 0l;
  for (long i = 0; i < this->NbrMatrixElements; ++i)
    if (this->MatrixElements[i] != 0.0)
      ++NbrNonZero;
  return NbrNonZero;
}

// compute the total amount of memory needed to store the sparse matrix
//
// return value = amount of memory (in bytes)

unsigned long SparseRealMatrix::GetAllocatedMemory()
{
  return (((sizeof(double)+ sizeof(int)) * this->NbrMatrixElements)
	  + ((2ul * sizeof(int)) * this->NbrRow));
}

// evaluate the real part of the matrix trace
//
// return value = real part of the matrix trace 

double SparseRealMatrix::Tr ()
{
  double Trace = 0.0;
  long TmpIndex;
  for (long i = 0; i < this->NbrRow; ++i)
    if (this->RowPointers[i] >= 0l)
      {
	TmpIndex = this->FindColumnIndexPosition(i, this->RowPointers[i], this->RowLastPointers[i]);
	if (TmpIndex >= 0l)
	  Trace += this->MatrixElements[TmpIndex];
      }
  return Trace;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const SparseRealMatrix& P) 
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      if (P.RowPointers[i] == -1l)
	{
	  for (int j = 0; j < P.NbrColumn; ++j)
	    Str << "0    ";
	}
      else
	{
	  long MinPos = P.RowPointers[i];
	  long MaxPos = P.RowLastPointers[i];
	  for (int j = 0; j < P.NbrColumn; ++j)
	    {
	      if ((MinPos <= MaxPos) && (P.ColumnIndices[MinPos] == j))
		{
		  Str << P.MatrixElements[MinPos] << "    ";
		  ++MinPos;
		}
	      else
		Str << "0    ";
	    }
	}
      Str << endl;
    }
  return Str;
}

// output the matrix in a sparse display (column formatted output)
//
// str = reference on output stream
// return value = reference on output stream

ostream& SparseRealMatrix::PrintNonZero (ostream& str) 
{
  for (long i = 0; i < this->NbrRow; ++i)
    if (this->RowPointers[i] >= 0l)
      {
	long MinPos = this->RowPointers[i];
	long MaxPos = this->RowLastPointers[i];
	while (MinPos <= MaxPos)
	  {
	    str << i << " " << this->ColumnIndices[MinPos] << " " << this->MatrixElements[MinPos] << endl;
	    ++MinPos;
	  }
      }
  return str;
}

#ifdef USE_OUTPUT

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// P = matrix to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const SparseRealMatrix& P) 
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; ++i)
    {
      Str << "{";
      if (P.RowPointers[0] == -1l)
	{
	  Str << "0,";
	}
      else
	{
	  long MinPos = P.RowPointers[i];
	  long MaxPos = P.RowLastPointers[i];
	  for (int j = 0; j < (P.NbrColumn - 1); ++j)
	    {
	      if ((MinPos <= MaxPos) && (P.ColumnIndices[MinPos] == j))
		{
		  Str << P.MatrixElements[MinPos] << ",";
		  ++MinPos;
		}
	      else
		Str << "0,";
	    }
	  if ((MinPos <= MaxPos) && (P.ColumnIndices[MinPos] == (P.NbrColumn - 1)))
	    {
	      Str << P.MatrixElements[MinPos];
	      ++MinPos;
	    }
	  else
	    Str << "0";	  
	  Str << "}";
	}
      Str << "}";
      if (i != (P.NbrRow - 1))
	Str << ",";
    }
  Str << "}";
  return Str;
}

#endif


