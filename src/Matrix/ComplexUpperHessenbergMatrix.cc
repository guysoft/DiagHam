////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of complex upper hessenberg matrix                 //
//                                                                            //
//                        last modification : 27/11/2012                      //
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


#include "Matrix/ComplexUpperHessenbergMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "GeneralTools/ListIterator.h"
#include "MathTools/Complex.h"

#include <stdlib.h>
#include <fstream>


using std::endl;


#ifdef HAVE_LAPACK

// binding to the LAPACK function
//
// extern "C" void FORTRAN_NAME(zgeev)(const char* jobVL, const char* jobVR, const int* nbrColumn, const double* matrix, const int* leadingDimension,
// 				    const double* eigenvaluesRealPart, const double* eigenvaluesImaginaryPart, 
// 				    const double* eigenvectorLeftMatrix, const int* leadingDimensionEigenvectorLeftMatrix,
// 				    const double* eigenvectorRightMatrix, const int* leadingDimensionEigenvectorRightMatrix,
// 				    const double* workingArea, const int* workingAreaSize, const int* information);

#endif


// default constructor
//

ComplexUpperHessenbergMatrix::ComplexUpperHessenbergMatrix() 
{
  this->UpperOffDiagonalElements = 0;
  this->DiagonalElements = 0;
  this->LowerDiagonalElements = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hessenberg | Matrix::Upper;
  this->Dummy = 0.0;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

ComplexUpperHessenbergMatrix::ComplexUpperHessenbergMatrix(int dimension, bool zero) 
{
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hessenberg | Matrix::Upper;
  this->DiagonalElements = new Complex [this->NbrRow];
  this->LowerDiagonalElements = new Complex [this->NbrRow];
  long TmpNbrOffDiagonalElements = (((long) this->NbrRow) * (((long) this->NbrRow) - 1l)) / 2l;
  this->UpperOffDiagonalElements = new Complex [TmpNbrOffDiagonalElements];
  if (zero == true)
    {
      long pos = 0;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  this->DiagonalElements[i] = 0.0;
	  this->LowerDiagonalElements[i] = 0.0;
	  for (int j = i + 1; j < this->NbrRow; j++)
	    {
	      this->UpperOffDiagonalElements[pos] = 0.0;
	      pos++;
	    }
	}
    }
  this->Dummy = 0.0;
}

// constructor from matrix elements (without duplicating datas)
//
// diagonalElements = diagonal elements
// offDiagonalElements = upper off-diagonal elements
// lowerDiagonalElements = lower diagonal elements
// dimension = matrix dimension

ComplexUpperHessenbergMatrix::ComplexUpperHessenbergMatrix(Complex* diagonalElements, Complex* offDiagonalElements, 
							   Complex* lowerDiagonalElements, int dimension) 
{
  this->DiagonalElements = diagonalElements;
  this->LowerDiagonalElements = lowerDiagonalElements;
  this->UpperOffDiagonalElements = offDiagonalElements;
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hessenberg | Matrix::Upper;
  this->Dummy = 0.0;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

ComplexUpperHessenbergMatrix::ComplexUpperHessenbergMatrix(const ComplexUpperHessenbergMatrix& M) 
{
  this->DiagonalElements = M.DiagonalElements;
  this->LowerDiagonalElements = M.LowerDiagonalElements;
  this->UpperOffDiagonalElements = M.UpperOffDiagonalElements;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hessenberg | Matrix::Upper;
  this->Dummy = 0.0;
}

// constructor from a matrix, copying the data and discarding all elements below the lower diagonal
//
// M = reference on the matrix

ComplexUpperHessenbergMatrix::ComplexUpperHessenbergMatrix(Matrix& M) 
{
  this->Flag.Initialize();  
  if (M.GetNbrColumn() >= M.GetNbrRow())
    {
      this->NbrRow = M.GetNbrRow();
     this->NbrColumn = M.GetNbrRow();
     }
  else
    {    
      this->NbrRow = M.GetNbrColumn();
      this->NbrColumn = M.GetNbrColumn();
    }
  this->TrueNbrRow = this->TrueNbrRow;
  this->TrueNbrColumn = this->TrueNbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hessenberg | Matrix::Upper;
  this->Dummy = 0.0;
  this->DiagonalElements = new Complex [this->NbrRow];
  this->LowerDiagonalElements = new Complex [this->NbrRow];
  long TmpNbrOffDiagonalElements = (((long) this->NbrRow) * (((long) this->NbrRow) - 1l)) / 2l;
  this->UpperOffDiagonalElements = new Complex [TmpNbrOffDiagonalElements];
  Complex Tmp;
  long pos = 0;
  for (int i = 0; i < this->NbrRow; i++)
    {
      M.GetMatrixElement(i, i, Tmp);
      this->DiagonalElements[i] = Tmp;
      if (i > 0)
	M.GetMatrixElement(i, i - 1, Tmp);
      this->LowerDiagonalElements[i] = Tmp;
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  M.GetMatrixElement(i, j, Tmp);
 	  this->UpperOffDiagonalElements[pos] = Tmp;
	  ++pos;
	}
    }  
}

// destructor
//

ComplexUpperHessenbergMatrix::~ComplexUpperHessenbergMatrix() 
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      if (this->DiagonalElements != 0)
	{
	  delete[] this->DiagonalElements;
	}
      if (this->LowerDiagonalElements != 0)
	{
	  delete[] this->LowerDiagonalElements;
	}
      if (this->UpperOffDiagonalElements != 0)
	{
	  delete[] this->UpperOffDiagonalElements;
	}
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

ComplexUpperHessenbergMatrix& ComplexUpperHessenbergMatrix::operator = (const ComplexUpperHessenbergMatrix& M) 
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      if (this->DiagonalElements != 0)
	{
	  delete[] this->DiagonalElements;
	}
      if (this->LowerDiagonalElements != 0)
	{
	  delete[] this->LowerDiagonalElements;
	}
      if (this->UpperOffDiagonalElements != 0)
	{
	  delete[] this->UpperOffDiagonalElements;
	}
    }
  this->DiagonalElements = M.DiagonalElements;
  this->LowerDiagonalElements = M.LowerDiagonalElements;
  this->UpperOffDiagonalElements = M.UpperOffDiagonalElements;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::ComplexElements | Matrix::Hessenberg | Matrix::Upper;
  this->Dummy = 0.0;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* ComplexUpperHessenbergMatrix::Clone ()
{
  return ((Matrix*) new ComplexUpperHessenbergMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexUpperHessenbergMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (i < j)
	{
	  this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] = x;
	}
      else
	{
	  if ((j + 1) == i)
	    {
	      this->LowerDiagonalElements[i] = x;
	    }
	}
    }      
  return;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexUpperHessenbergMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if (i == j)
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (i < j)
	{
	  this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] = x;
	}
      else
	{
	  if ((j + 1) == i)
	    {
	      this->LowerDiagonalElements[i] = x;
	    }
	}
    }
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void ComplexUpperHessenbergMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (i < j)
	{
	  this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] += x;
	}
      else
	{
	  if (i == (j + 1))
	    {
	      this->LowerDiagonalElements[i] += x;
	    }
	}
    }      
  return;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void ComplexUpperHessenbergMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn))
    return;
  if (i == j)
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (i < j)
	{
	  this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] += x;
	}
      else
	{
	  if ((j + 1) == i)
	    {
	      this->LowerDiagonalElements[i] += x;
	    }
	}
    }      
  return;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexUpperHessenbergMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  Complex* TmpDiagonalElements = new Complex [nbrRow];
  Complex* TmpLowerDiagonalElements = new Complex [nbrRow];
  long TmpNbrOffDiagonalElements = (((long) nbrRow) * (((long) nbrRow) - 1l)) / 2l;
  Complex* TmpUpperOffDiagonalElements = new Complex [TmpNbrOffDiagonalElements];
  int MinNbrRow = nbrRow;
  if (nbrRow > this->NbrRow)
    MinNbrRow = this->NbrRow;
  for (int i = 0; i < MinNbrRow; ++i)
    {
      TmpDiagonalElements[i] = this->DiagonalElements[i];
      TmpLowerDiagonalElements[i] = this->LowerDiagonalElements[i];
      for (int j = i + 1; j < MinNbrRow; ++j)
	{
	  long Index = i + (j * (j - 1l)) / 2l;
	  TmpUpperOffDiagonalElements[Index] = this->UpperOffDiagonalElements[Index];
	}
    }
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      if (this->DiagonalElements != 0)
	{
	  delete[] this->DiagonalElements;
	}
      if (this->LowerDiagonalElements != 0)
	{
	  delete[] this->LowerDiagonalElements;
	}
      if (this->UpperOffDiagonalElements != 0)
	{
	  delete[] this->UpperOffDiagonalElements;
	}
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->DiagonalElements = TmpDiagonalElements;
  this->LowerDiagonalElements = TmpLowerDiagonalElements;
  this->UpperOffDiagonalElements = TmpUpperOffDiagonalElements;
  this->Flag.Initialize();
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexUpperHessenbergMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  Complex* TmpDiagonalElements = new Complex [nbrRow];
  Complex* TmpLowerDiagonalElements = new Complex [nbrRow];
  long TmpNbrOffDiagonalElements = (((long) nbrRow) * (((long) nbrRow) - 1l)) / 2l;
  Complex* TmpUpperOffDiagonalElements = new Complex [TmpNbrOffDiagonalElements];
  int MinNbrRow = nbrRow;
  if (nbrRow > this->NbrRow)
    MinNbrRow = this->NbrRow;
  for (int i = 0; i < MinNbrRow; ++i)
    {
      TmpDiagonalElements[i] = this->DiagonalElements[i];
      TmpLowerDiagonalElements[i] = this->LowerDiagonalElements[i];
      for (int j = i + 1; j < MinNbrRow; ++j)
	{
	  long Index = i + (j * (j - 1l)) / 2l;
	  TmpUpperOffDiagonalElements[Index] = this->UpperOffDiagonalElements[Index];
	}
      for (int j = MinNbrRow; j < nbrRow; ++j)
	{
	  long Index = i + (j * (j - 1l)) / 2l;
	  TmpUpperOffDiagonalElements[Index] = 0.0;
	}
    }
  for (int i = this->NbrRow; i < nbrRow; ++i)
    {
      TmpDiagonalElements[i] = 0.0;
      TmpLowerDiagonalElements[i] = 0.0;
      for (int j = i + 1; j < nbrRow; ++j)
	{
	  TmpUpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] = 0.0;
	}
    }
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      if (this->DiagonalElements != 0)
	{
	  delete[] this->DiagonalElements;
	}
      if (this->LowerDiagonalElements != 0)
	{
	  delete[] this->LowerDiagonalElements;
	}
      if (this->UpperOffDiagonalElements != 0)
	{
	  delete[] this->UpperOffDiagonalElements;
	}
    }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->DiagonalElements = TmpDiagonalElements;
  this->LowerDiagonalElements = TmpLowerDiagonalElements;
  this->UpperOffDiagonalElements = TmpUpperOffDiagonalElements;
  this->Flag.Initialize();
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

ComplexUpperHessenbergMatrix operator + (const ComplexUpperHessenbergMatrix& M1, const ComplexUpperHessenbergMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexUpperHessenbergMatrix();
  Complex* TmpDiagonalElements = new Complex [M1.NbrRow];
  Complex* TmpLowerDiagonalElements = new Complex [M1.NbrRow];
  long TmpNbrOffDiagonalElements = (((long) M1.NbrRow) * (((long) M1.NbrRow) - 1l)) / 2l;
  Complex* TmpUpperOffDiagonalElements = new Complex [TmpNbrOffDiagonalElements];
  long pos = 0;
  for (int i = 0; i < M1.NbrRow; i++)
    {
      TmpDiagonalElements[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
      TmpLowerDiagonalElements[i] = M1.LowerDiagonalElements[i] + M2.LowerDiagonalElements[i];
      for (int j = i + 1; j < M1.NbrRow; j++)
	{
	  TmpUpperOffDiagonalElements[pos] = M1.UpperOffDiagonalElements[pos] + M2.UpperOffDiagonalElements[pos];
	  ++pos;
	}
    }
  return ComplexUpperHessenbergMatrix(TmpDiagonalElements, TmpUpperOffDiagonalElements, TmpLowerDiagonalElements, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexUpperHessenbergMatrix operator - (const ComplexUpperHessenbergMatrix& M1, const ComplexUpperHessenbergMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexUpperHessenbergMatrix();
  Complex* TmpDiagonalElements = new Complex [M1.NbrRow];
  Complex* TmpLowerDiagonalElements = new Complex [M1.NbrRow];
  long TmpNbrOffDiagonalElements = (((long) M1.NbrRow) * (((long) M1.NbrRow) - 1l)) / 2l;
  Complex* TmpUpperOffDiagonalElements = new Complex [TmpNbrOffDiagonalElements];
  long pos = 0;
  for (int i = 0; i < M1.NbrRow; i++)
    {
      TmpDiagonalElements[i] = M1.DiagonalElements[i] - M2.DiagonalElements[i];
      TmpLowerDiagonalElements[i] = M1.LowerDiagonalElements[i]  - M2.LowerDiagonalElements[i];
      for (int j = i + 1; j < M1.NbrRow; j++)
	{
	  TmpUpperOffDiagonalElements[pos] = M1.UpperOffDiagonalElements[pos] - M2.UpperOffDiagonalElements[pos];
	  ++pos;
	}
    }
  return ComplexUpperHessenbergMatrix(TmpDiagonalElements, TmpUpperOffDiagonalElements, TmpLowerDiagonalElements, M1.NbrRow);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexUpperHessenbergMatrix operator * (const ComplexUpperHessenbergMatrix& M, double x) 
{
  Complex* TmpDiagonalElements = new Complex [M.NbrRow];
  Complex* TmpLowerDiagonalElements = new Complex [M.NbrRow];
  long TmpNbrOffDiagonalElements = (((long) M.NbrRow) * (((long) M.NbrRow) - 1l)) / 2l;
  Complex* TmpUpperOffDiagonalElements = new Complex [TmpNbrOffDiagonalElements];
  long pos = 0;
  for (int i = 0; i < M.NbrRow; i++)
    {
      TmpDiagonalElements[i] = M.DiagonalElements[i] * x;
      TmpLowerDiagonalElements[i] = M.LowerDiagonalElements[i] * x;
      for (int j = i + 1; j < M.NbrRow; j++)
	{
	  TmpUpperOffDiagonalElements[pos] = M.UpperOffDiagonalElements[pos] * x;
	  ++pos;
	}
    }
  return ComplexUpperHessenbergMatrix(TmpDiagonalElements, TmpUpperOffDiagonalElements, TmpLowerDiagonalElements, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexUpperHessenbergMatrix operator * (double x, const ComplexUpperHessenbergMatrix& M) 
{
  return (M * x);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

ComplexUpperHessenbergMatrix operator / (const ComplexUpperHessenbergMatrix& M, double x) 
{
  x = 1.0 / x;
  return (M * x);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexUpperHessenbergMatrix& ComplexUpperHessenbergMatrix::operator += (const ComplexUpperHessenbergMatrix& M) 
{
  if (this->NbrRow == 0)
    return *this;
  long pos = 0;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
      this->LowerDiagonalElements[i] += M.LowerDiagonalElements[i];
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  this->UpperOffDiagonalElements[pos] += M.UpperOffDiagonalElements[pos];
	  ++pos;
	}
    }
  return *this;
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to add to current matrix
// return value = reference on current matrix

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexUpperHessenbergMatrix& ComplexUpperHessenbergMatrix::operator -= (const ComplexUpperHessenbergMatrix& M) 
{
  if (this->NbrRow == 0)
    return *this;
  long pos = 0;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
      this->LowerDiagonalElements[i] -= M.LowerDiagonalElements[i];
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  this->UpperOffDiagonalElements[pos] -= M.UpperOffDiagonalElements[pos];
	  ++pos;
	}
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexUpperHessenbergMatrix& ComplexUpperHessenbergMatrix::operator *= (double x) 
{
  long pos = 0;
  for (int i = 0; i < this->NbrRow; i++)
    {
      this->DiagonalElements[i] *= x;
      this->LowerDiagonalElements[i] *= x;
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  this->UpperOffDiagonalElements[pos] *= x;
	  ++pos;
	}
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexUpperHessenbergMatrix& ComplexUpperHessenbergMatrix::operator /= (double x)
{
  x = 1.0 / x;
  (*this) *= x;
  return *this;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ComplexUpperHessenbergMatrix::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  Complex x = 0.0;
  if ((V1.GetVectorDimension() != this->NbrRow) || (V2.GetVectorDimension() != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      Complex x2 = this->DiagonalElements[i] * V2[i];
      x2 += this->LowerDiagonalElements[i] * V2[i - 1];
      for (int j = i + 1; j < this->NbrColumn; ++j)
	{
	  x2 +=  this->UpperOffDiagonalElements[i + (j * (j - 1l)) / 2l] * V2[j];
	}
      x += Conj(V1[i]) * x2;
    }
  return x;
}

// evaluate matrix trace
//
// return value = matrix trace 

double ComplexUpperHessenbergMatrix::Tr () 
{
  return 0.0;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double ComplexUpperHessenbergMatrix::Det () 
{
  return 0.0;
}

// Diagonalize a real matrix using the LAPACK library
//
// M = reference on complex diagonal matrix where result has to be stored
// leftFlag = compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// return value = reference on complex diagonal matrix

ComplexDiagonalMatrix& ComplexUpperHessenbergMatrix::LapackDiagonalize (ComplexDiagonalMatrix& M, bool leftFlag)
{
#ifdef HAVE_LAPACK
//   int Information = 0;
//   int WorkingAreaSize = -1;
//   char JobVL;
//   char JobVR;
//   if (leftFlag == true)
//     {
//       JobVL = 'V';
//       JobVR = 'N';
//     }
//   else
//     {
//       JobVL = 'N';
//       JobVR = 'V';
//     }
//   double* TmpMatrix = new double [this->NbrColumn * this->NbrRow];
//   long TotalIndex = 0l;
//   for (int j = 0; j < this->NbrColumn; ++j)
//     {
//       for (int i = 0; i < this->NbrRow; ++i)
// 	{
// 	  TmpMatrix[TotalIndex] = this->Columns[j][i];
// 	  ++TotalIndex;
// 	}
//     }
//   double* TmpEigenvalueReal = new double[this->NbrColumn];
//   double* TmpEigenvalueImaginary = new double[this->NbrColumn];
//   int TmpLeadingDimension = 1;
//   double* Dummy = 0;
//   double TmpWorkingArea;
//   FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
// 		      TmpEigenvalueReal, TmpEigenvalueImaginary,
// 		      Dummy, &TmpLeadingDimension, Dummy, &TmpLeadingDimension,
// 		      &TmpWorkingArea, &WorkingAreaSize, &Information);
//   WorkingAreaSize = (int) TmpWorkingArea;
//   double* WorkingArea = new double [WorkingAreaSize];
//   FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
// 		      TmpEigenvalueReal, TmpEigenvalueImaginary,
// 		      Dummy, &TmpLeadingDimension, Dummy, &TmpLeadingDimension,
// 		      WorkingArea, &WorkingAreaSize, &Information);
//   for (int i = 0; i < this->NbrRow; ++i)
//     {
//       M[i].Re = TmpEigenvalueReal[i];
//       M[i].Im = TmpEigenvalueImaginary[i];
//     }
//   delete[] WorkingArea;
//   delete[] TmpEigenvalueReal;
//   delete[] TmpEigenvalueImaginary;
//   delete[] TmpMatrix;
#endif
  return M;
}

// Diagonalize a real matrix and evaluate the left eigenstates using the LAPACK library
//
// M = reference on complex diagonal matrix where result has to be stored
// Q = matrix where transformation matrix has to be stored
// leftFlag = compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// return value = reference on complex diagonal matrix

ComplexDiagonalMatrix& ComplexUpperHessenbergMatrix::LapackDiagonalize (ComplexDiagonalMatrix& M, ComplexMatrix& Q, bool leftFlag)
{
#ifdef HAVE_LAPACK
//   int Information = 0;
//   int WorkingAreaSize = -1;
//   char JobVL;
//   char JobVR;
//   if (leftFlag == true)
//     {
//       JobVL = 'V';
//       JobVR = 'N';
//     }
//   else
//     {
//       JobVL = 'N';
//       JobVR = 'V';
//     }
//   double* TmpMatrix = new double [this->NbrColumn * this->NbrRow];
//   double* TmpLeftEigenstates = new double [this->NbrColumn * this->NbrRow];
//   long TotalIndex = 0l;
//   for (int j = 0; j < this->NbrColumn; ++j)
//     {
//       for (int i = 0; i < this->NbrRow; ++i)
// 	{
// 	  TmpMatrix[TotalIndex] = this->Columns[j][i];
// 	  ++TotalIndex;
// 	}
//     }
//   double* TmpEigenvalueReal = new double[this->NbrColumn];
//   double* TmpEigenvalueImaginary = new double[this->NbrColumn];
//   int TmpLeadingLeftDimension;
//   int TmpLeadingRightDimension;
//   if (leftFlag == true)
//     {
//       TmpLeadingLeftDimension = this->NbrColumn;
//       TmpLeadingRightDimension = 1;
//     }
//   else
//     {
//       TmpLeadingRightDimension = this->NbrColumn;
//       TmpLeadingLeftDimension = 1;
//     }
//   double* Dummy = 0;
//   double TmpWorkingArea;
//   if (leftFlag == true)
//     {
//       FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
// 			  TmpEigenvalueReal, TmpEigenvalueImaginary,
// 			  TmpLeftEigenstates, &TmpLeadingLeftDimension, Dummy, &TmpLeadingRightDimension, 
// 			  &TmpWorkingArea, &WorkingAreaSize, &Information);
//     }
//   else
//     {
//       FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
// 			  TmpEigenvalueReal, TmpEigenvalueImaginary,
// 			  Dummy, &TmpLeadingLeftDimension, TmpLeftEigenstates, &TmpLeadingRightDimension, 
// 			  &TmpWorkingArea, &WorkingAreaSize, &Information);
//     }
//   WorkingAreaSize = (int) TmpWorkingArea;
//   double* WorkingArea = new double [WorkingAreaSize];
//   if (leftFlag == true)
//     {
//       FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
// 			  TmpEigenvalueReal, TmpEigenvalueImaginary,
// 			  TmpLeftEigenstates, &TmpLeadingLeftDimension, Dummy, &TmpLeadingRightDimension, 
// 			  WorkingArea, &WorkingAreaSize, &Information);
//     }
//   else
//     {
//       FORTRAN_NAME(dgeev)(&JobVL, &JobVR, &this->NbrRow, TmpMatrix, &this->NbrRow, 
// 			  TmpEigenvalueReal, TmpEigenvalueImaginary,
// 			  Dummy, &TmpLeadingLeftDimension, TmpLeftEigenstates, &TmpLeadingRightDimension, 
// 			  WorkingArea, &WorkingAreaSize, &Information);
//     }
//   for (int i = 0; i < this->NbrRow; ++i)
//     {
//       M[i].Re = TmpEigenvalueReal[i];
//       M[i].Im = TmpEigenvalueImaginary[i];
//     }
//   TotalIndex = 0l;
//   for (int i = 0; i < this->NbrRow;)
//     {
//       if ((i == (this->NbrRow - 1)) || (TmpEigenvalueImaginary[i] != -TmpEigenvalueImaginary[i + 1]) 
// 	  || (TmpEigenvalueImaginary[i] == 0.0))
// 	{
// 	  Complex Tmp;
// 	  for (int j = 0; j < this->NbrRow; ++j)
// 	    {
// 	      Tmp.Re = TmpLeftEigenstates[TotalIndex];
// 	      Q.SetMatrixElement(j, i, Tmp);
// 	      ++TotalIndex;
// 	    }
// 	  ++i;
// 	}
//       else
// 	{
// 	  Complex Tmp;
// 	  for (int j = 0; j < this->NbrRow; ++j)
// 	    {
// 	      Tmp.Re = TmpLeftEigenstates[TotalIndex];
// 	      Q.SetMatrixElement(j, i, Tmp);
// 	      Q.SetMatrixElement(j, i + 1, Tmp);
// 	      ++TotalIndex;
// 	    }
// 	  for (int j = 0; j < this->NbrRow; ++j)
// 	    {
// 	      Q.GetMatrixElement(j, i, Tmp);
// 	      Tmp.Im = TmpLeftEigenstates[TotalIndex];
// 	      Q.SetMatrixElement(j, i, Tmp);
// 	      Tmp.Im = -TmpLeftEigenstates[TotalIndex];
// 	      Q.SetMatrixElement(j, i + 1, Tmp);
// 	      ++TotalIndex;
// 	    }
// 	  i += 2;
// 	}
//     }
//   delete[] TmpLeftEigenstates;
//   delete[] WorkingArea;
//   delete[] TmpEigenvalueReal;
//   delete[] TmpEigenvalueImaginary;
//   delete[] TmpMatrix;
#endif
  return M;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const ComplexUpperHessenbergMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      for (int j = 0; j < (i - 1); j ++)
	{
	  Str << 0.0 << "    ";
	}
      if (i > 0)
	{
	  Str << P.LowerDiagonalElements[i];
	}
      Str << P.DiagonalElements[i] << "    ";
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << P.UpperOffDiagonalElements[i + ((j - 1l) * j) / 2l] << "    ";
	}
      Str << endl;
    }
  return Str;
}

#ifdef USE_OUTPUT

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// P = matrix to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexUpperHessenbergMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; ++i)
    {
      Str << "{";
      int pos = i - 1;
      for (int j = 0; j < (i - 1); ++j)
	{
	  Str << "0,";
	}
      if (i > 0)
	{
	  Str << P.LowerDiagonalElements[i];
	}
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << "," << P.UpperOffDiagonalElements[i + ((j - 1l) * j) / 2l];
	}
      if (i != (P.NbrRow - 1))
	Str << "},";
      else
	Str << "}";
    }
  Str << "}";
  return Str;
}

#endif
