////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of complex matrix with sparse storage                //
//                                                                            //
//                        last modification : 04/10/2012                      //
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


#ifndef SPARSECOMPLEXMATRIX_H
#define SPARSECOMPLEXMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Vector/ComplexVector.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class SparseComplexMatrix : public Matrix
{

  friend class ComplexMatrix;
  friend class RealVector;
  friend class ComplexVector;

 protected:

  // number of non-zero mtarix elements
  long NbrMatrixElements; 

  // array that contains the matrix elements
  Complex* MatrixElements; 

  // array that give the column index of each matrix element
  int* ColumnIndices;

  // array that gives the entry point of each row in MatrixElements
  long* RowPointers;
  // array that gives the entry point last element of each row in MatrixElements
  long* RowLastPointers;
  
  // garbage collector flag
  GarbageFlag Flag;

 public:

  // default constructor
  //
  SparseComplexMatrix();

  // constructor for a sparse matrix without any specific struture but a given number of non-zero matrix elements
  //
  // nbrRow = number of rows
  // nbrColumn = number of columns
  // nbrMatrixElements = number of non-zero matrix elements
  // zero = true if matrix elements have to be set to zero
  SparseComplexMatrix(int nbrRow, int nbrColumn, long nbrMatrixElements, bool zero = false);

  // copy constructor from a complex matrix
  //
  // M = matrix to copy
  SparseComplexMatrix(const SparseComplexMatrix& M);

  // copy constructor (duplicating all datas)
  //
  // M = matrix to copy
  SparseComplexMatrix(Matrix& M);

  // destructor
  //
  ~SparseComplexMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  SparseComplexMatrix& operator = (const SparseComplexMatrix& M);

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  Matrix* Clone ();  

  // copy a matrix into another (duplicating data)
  //
  // matrix = matrix to copy
  // return value = reference on current matrix
  SparseComplexMatrix& Copy (SparseComplexMatrix& matrix);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, double x);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, const Complex& x);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // real = new real value for matrix element
  // imag = new imaginary value for matrix element
  void SetMatrixElement(int i, int j, double real, double imag);

  // get a matrix element (real part if complex)
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, double& x) const;

  // get a matrix element
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, Complex& x) const;

  // add a value to a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  void AddToMatrixElement(int i, int j, double x);

  // add a value  a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  void AddToMatrixElement(int i, int j, const Complex& x);

  // Resize matrix
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  void Resize (int nbrRow, int nbrColumn);

  // Resize matrix and set to zero all elements that have been added
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  void ResizeAndClean (int nbrRow, int nbrColumn);

  // Set all entries in matrix to zero
  //
  void ClearMatrix ();

  // add two matrices
  //
  // matrix1 = first matrix
  // matrix2 = second matrix
  // return value = sum of the two matrices
  friend SparseComplexMatrix operator + (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2);

  // difference of two matrices
  //
  // matrix1 = first matrix
  // matrix2 = second matrix
  // return value = difference of the two matrices
  friend SparseComplexMatrix operator - (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2);

  // create the linear combination of two matrices
  //
  // x1 = prefactor of the first matrix
  // matrix1 = first matrix
  // x2 = prefactor of the second matrix
  // matrix2 = second matrix
   // return value = linear combination
  friend SparseComplexMatrix SparseComplexMatrixLinearCombination(const Complex& x1, const SparseComplexMatrix& matrix1, const Complex& x2, const SparseComplexMatrix& matrix2);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  SparseComplexMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  SparseComplexMatrix& operator /= (double x);
  
  // multiply a matrix to the right by another matrix
  //
  // matrix = matrix used as multiplicator
  // return value = reference on current matrix
  SparseComplexMatrix& Multiply (const SparseComplexMatrix& matrix);

  // multiply a matrix to the right by another matrix, providing all the required temporary arrays
  //
  // matrix = matrix used as multiplicator
  // tmpMatrixElements = temporary array of complex numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of complex numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  SparseComplexMatrix& Multiply (const SparseComplexMatrix& matrix, Complex* tmpMatrixElements, 
				 int* tmpColumnIndices, Complex* tmpElements);

  // multiply a matrix to the right by another matrix, providing all the required temporary arrays, extend their capacity if needed
  //
  // matrix = matrix used as multiplicator
  // tmpMatrixElements = reference on the temporary array of complex numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = reference on the temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // nbrElements = reference ont the number of elements in tmpMatrixElements and tmpColumnIndices
  // tmpElements = temporary array of complex numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  SparseComplexMatrix& Multiply (const SparseComplexMatrix& matrix, Complex*& tmpMatrixElements, int*& tmpColumnIndices, 
				 long& nbrElements, Complex* tmpElements);

  // compute the number of non-zero matrix elements (zero having strictly zero square norm)
  //
  // return value = number of non-zero matrix elements
  long ComputeNbrNonZeroMatrixElements();

  // compute the total amount of memory needed to store the sparse matrix
  //
  // return value = amount of memory (in bytes)
  unsigned long GetAllocatedMemory();

  // evaluate the real part of the matrix trace
  //
  // return value = real part of the matrix trace 
  double Tr ();

  // evaluate the matrix trace
  //
  // return value = matrix trace 
  Complex ComplexTr ();

  // compute the tensor product of two sparse matrices (matrix1 x matrix2), and store the result in a sparse matrix
  //
  // matrix1 = reference on the left matrix
  // matrix2 = reference on the right matrix
  // return value = tensor product
  friend SparseComplexMatrix TensorProduct (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2);

  // compute the hermitian transpose of the current matrix
  //
  // return value = hermitian transposed matrix
  SparseComplexMatrix HermitianTranspose ();

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const SparseComplexMatrix& P);

#ifdef USE_OUTPUT

 // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const SparseComplexMatrix& P);

#endif

 public:

  // find the position of a given column index in a row
  //
  // index = column index
  // minPosition = minimum position of the row in the compressed row storage
  // maxPosition = maximum position of the row in the compressed row storage
  // return value = position of the column index (-1 if it does not exist)
  long FindColumnIndexPosition(int index, long minPosition, long maxPosition) const;

  // get the matrix element
  // i = position

  Complex GetMatrixElement(int i);

  // get the column index
  // i = position

  int GetColumnIndex(int i);

  // get the nbr of matrix elements
  // i = position

  long GetNbrMatrixElements();

  //returns the array with indices of rows

  void GetRowIndices(int* RowIndices);

};

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void SparseComplexMatrix::GetMatrixElement(int i, int j, double& x) const
{
  if (this->RowPointers[i] == -1l)
    {
      x = 0.0;
      return;
    }
  long TmpIndex = this->FindColumnIndexPosition(j, this->RowPointers[i], this->RowLastPointers[i]);
  if (TmpIndex == -1l)
    {
      x = 0.0;
      return;      
    }
  x = this->MatrixElements[TmpIndex].Re;
  return;
}

// get the matrix element
// i = position

inline Complex SparseComplexMatrix::GetMatrixElement(int i)
{
  return this->MatrixElements[i];
}

// get the column index
// i = position

inline int SparseComplexMatrix::GetColumnIndex(int i)
{
  return this->ColumnIndices[i];
}

// get the nbr of matrix elements
// i = position

inline long SparseComplexMatrix::GetNbrMatrixElements()
{
  return this->NbrMatrixElements;
}


// get a matrix element
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void SparseComplexMatrix::GetMatrixElement(int i, int j, Complex& x) const
{
  if (this->RowPointers[i] == -1l)
    {
      x = 0.0;
      return;
    }
  long TmpIndex = this->FindColumnIndexPosition(j, this->RowPointers[i], this->RowLastPointers[i]);
  if (TmpIndex == -1l)
    {
      x = 0.0;
      return;      
    }
  x = this->MatrixElements[TmpIndex];
  return;
}

// find the position of a given column index in a row
//
// index = column index
// minPosition = minimum position of the row in the compressed row storage
// maxPosition = maximum position of the row in the compressed row storage
// return value = position of the column index (-1 if it does not exist)

inline long SparseComplexMatrix::FindColumnIndexPosition(int index, long minPosition, long maxPosition) const
{
  while (minPosition <= maxPosition)
    {
      if (this->ColumnIndices[minPosition] == index)
	return minPosition;
      ++minPosition;
    }
  return -1l;
}

#endif
