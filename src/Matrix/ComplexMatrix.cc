////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                             class of real matrix                           //
//                                                                            //
//                        last modification : 05/02/2001                      //
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


#include "Matrix/ComplexMatrix.h"
#include "Vector/ComplexVector.h"


using std::endl;
using std::cout;


// default constructor
//

ComplexMatrix::ComplexMatrix() 
{
  this->Columns = 0;
  this->ColumnGarbageFlag = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;  
  this->MatrixType = Matrix::ComplexElements;
}

// constructor for an empty matrix
//
// nbrRow = number of rows
// nbrColumn = number of columns
// zero = tue if matrix elements have to be set to zero

ComplexMatrix::ComplexMatrix(int nbrRow, int nbrColumn, bool zero)
{
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->NbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Columns = new ComplexVector [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] = ComplexVector (this->NbrRow, zero);
  this->MatrixType = Matrix::ComplexElements;
}

// constructor from matrix elements (without duplicating datas)
//
// columns = pointer an array of vector
// nbrColumn = number of columns

ComplexMatrix::ComplexMatrix(ComplexVector* columns, int nbrColumn) 
{
  this->Columns = columns;
  this->ColumnGarbageFlag = new int;
  *(this->ColumnGarbageFlag) = 1;
  this->NbrRow = columns[0].GetVectorDimension();
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;  
  this->MatrixType = Matrix::ComplexElements;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

ComplexMatrix::ComplexMatrix(const ComplexMatrix& M) 
{
  this->Columns = M.Columns;
  this->ColumnGarbageFlag = M.ColumnGarbageFlag;
  (*(this->ColumnGarbageFlag))++;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;  
  this->MatrixType = Matrix::ComplexElements;
}

// copy constructor (duplicating all datas)
//
// M = matrix to copy

ComplexMatrix::ComplexMatrix(Matrix& M)
{
  if ((M.GetNbrRow() == 0) || (M.GetNbrColumn() == 0))
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = 0;
      this->TrueNbrColumn = 0;
      this->MatrixType = Matrix::ComplexElements;
    }
  else
    {
      this->ColumnGarbageFlag = new int;
      *(this->ColumnGarbageFlag) = 1;
      this->NbrColumn = M.GetNbrColumn();
      this->NbrRow = M.GetNbrRow();
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;
      this->Columns = new ComplexVector [this->NbrColumn];
      Complex Tmp;
      for (int i = 0; i < this->NbrColumn; i++)
	{
	  this->Columns[i] = ComplexVector (this->NbrRow);
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      M.GetMatrixElement(j, i, Tmp);
	      this->Columns[i].Re(j) = Tmp.Re;
	      this->Columns[i].Im(j) = Tmp.Im;
	    }
	}
      this->MatrixType = Matrix::ComplexElements;
    }
}

// destructor
//

ComplexMatrix::~ComplexMatrix() 
{
  if (this->ColumnGarbageFlag != 0)
    {
      if ((*(this->ColumnGarbageFlag)) == 1)
	{
	  delete[] this->Columns;
	  delete this->ColumnGarbageFlag;
	}
      else
	(*(this->ColumnGarbageFlag))--;
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

ComplexMatrix& ComplexMatrix::operator = (const ComplexMatrix& M) 
{
  if (this->ColumnGarbageFlag != 0)
    {
      if ((*(this->ColumnGarbageFlag)) == 1)
	{
	  delete[] this->Columns;
	  delete this->ColumnGarbageFlag;
	}
      else
	(*(this->ColumnGarbageFlag))--;
    }
  if (M.ColumnGarbageFlag != 0)
    {
      this->Columns = M.Columns;
      this->ColumnGarbageFlag = M.ColumnGarbageFlag;
      (*(this->ColumnGarbageFlag))++;
      this->NbrRow = M.NbrRow;
      this->NbrColumn = M.NbrColumn;
      this->TrueNbrRow = M.TrueNbrRow;
      this->TrueNbrColumn = M.TrueNbrColumn;  
      this->MatrixType = Matrix::ComplexElements;
    }
  else
    {
      this->Columns = 0;
      this->ColumnGarbageFlag = 0;
      this->NbrRow = 0;
      this->NbrColumn = 0;
      this->TrueNbrRow = this->NbrRow;
      this->TrueNbrColumn = this->NbrColumn;  
      this->MatrixType = Matrix::ComplexElements;
    }
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* ComplexMatrix::Clone ()
{
  return ((Matrix*) new ComplexMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].RealComponents[i] = x;
  this->Columns[j].ImaginaryComponents[i] = 0.0;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].RealComponents[i] = x.Re;
  this->Columns[j].ImaginaryComponents[i] = x.Im;
}


// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void ComplexMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].RealComponents[i] += x;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void ComplexMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i > this->NbrRow) || (j > this->NbrColumn))
    return;
  this->Columns[j].RealComponents[i] += x.Re;
  this->Columns[j].ImaginaryComponents[i] += x.Im;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (this->NbrRow != nbrRow)
    {
      for (int i = 0; i < this->NbrColumn; i++)
	this->Columns[i].Resize(nbrRow);
      if (this->TrueNbrRow >= nbrRow)
	{
	  this->NbrRow = nbrRow;
	}
      else
	{
	  this->NbrRow = nbrRow;
	  this->TrueNbrRow = nbrRow;
	}
    }
  if (this->TrueNbrColumn >= nbrColumn)
    {
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	this->Columns[i].Resize(nbrRow);
      this->NbrColumn = nbrColumn;
    }
  else
    {
      ComplexVector* Tmp = new ComplexVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = ComplexVector(nbrRow);
      delete[] this->Columns;
      this->Columns = Tmp;
      this->TrueNbrColumn = nbrColumn;
      this->NbrColumn = nbrColumn;
    }
  return;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (this->NbrRow != nbrRow)
    {
      for (int i = 0; i < this->NbrColumn; i++)
	this->Columns[i].ResizeAndClean(nbrRow);
      if (this->TrueNbrRow >= nbrRow)
	{
	  this->NbrRow = nbrRow;
	}
      else
	{
	  this->NbrRow = nbrRow;
	  this->TrueNbrRow = nbrRow;
	}
    }
  if (this->TrueNbrColumn >= nbrColumn)
    {
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	this->Columns[i].ResizeAndClean(nbrRow);
      this->TrueNbrColumn = nbrColumn;
    }
  else
    {
      ComplexVector* Tmp = new ComplexVector[nbrColumn];
      for (int i = 0; i < this->NbrColumn; i++)
	Tmp[i] = this->Columns[i];      
      for (int i = this->NbrColumn; i < nbrColumn; i++)
	Tmp[i] = ComplexVector(nbrRow, true);
      delete[] this->Columns;
      this->Columns = Tmp;
      this->TrueNbrColumn = nbrColumn;
      this->NbrColumn = nbrColumn;
    }
  return;
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

ComplexMatrix operator + (const ComplexMatrix& M1, const ComplexMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] + M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j] + M2.Columns[i].ImaginaryComponents[j];
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// add two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

ComplexMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); ++j)
	{
	  TmpColumns[i].RealComponents[j] = M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M2.Columns[i].ImaginaryComponents[j];
	}
      if (i > 0)
	{
	  TmpColumns[i].RealComponents[j] = M1.UpperDiagonalElements[i - 1] + M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M2.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      TmpColumns[i].RealComponents[j] = M1.DiagonalElements[i] + M2.Columns[i].RealComponents[j];
      TmpColumns[i].ImaginaryComponents[j] = M2.Columns[i].ImaginaryComponents[j];
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].RealComponents[j] = M1.UpperDiagonalElements[i + 1] + M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M2.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M2.Columns[i].RealComponents[j];	
	  TmpColumns[i].ImaginaryComponents[j] = M2.Columns[i].ImaginaryComponents[j];	
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = left matrix
// M2 = right matrix
// return value = sum of the two matrices

ComplexMatrix operator + (const ComplexMatrix& M1, 
			  const RealTriDiagonalSymmetricMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); ++j)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	}
      if (i > 0)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] + M2.UpperDiagonalElements[i - 1];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] + M2.DiagonalElements[i];
      TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] + M2.UpperDiagonalElements[i + 1];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j];	
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];	
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexMatrix operator - (const ComplexMatrix& M1, const ComplexMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] - M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j] - M2.Columns[i].ImaginaryComponents[j];	  
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices where the left one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); ++j)
	{
	  TmpColumns[i].RealComponents[j] = -M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = -M2.Columns[i].ImaginaryComponents[j];
	}
      if (i > 0)
	{
	  TmpColumns[i].RealComponents[j] = M1.UpperDiagonalElements[i - 1] - M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = -M2.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      TmpColumns[i].RealComponents[j] = M1.DiagonalElements[i] - M2.Columns[i].RealComponents[j];
      TmpColumns[i].ImaginaryComponents[j] = -M2.Columns[i].ImaginaryComponents[j];
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].RealComponents[j] = M1.UpperDiagonalElements[i + 1] - M2.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = -M2.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].RealComponents[j] = -M2.Columns[i].RealComponents[j];	
	  TmpColumns[i].ImaginaryComponents[j] = -M2.Columns[i].ImaginaryComponents[j];	
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexMatrix operator - (const ComplexMatrix& M1, 
			  const RealTriDiagonalSymmetricMatrix& M2)
{
  if ((M1.NbrColumn != M2.NbrColumn) || (M1.NbrRow != M2.NbrRow))
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M1.NbrColumn];
  for (int i = 0; i < M1.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      int j = 0;
      for (; j < (i - 1); ++j)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	}
      if (i > 0)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] - M2.UpperDiagonalElements[i - 1];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] - M2.DiagonalElements[i];
      TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
      ++j;
      if (i < (M1.NbrColumn - 1))
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j] - M2.UpperDiagonalElements[i + 1];
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];
	  ++j;
	}
      ++j;
      for (; j < M1.NbrColumn; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M1.Columns[i].RealComponents[j];	
	  TmpColumns[i].ImaginaryComponents[j] = M1.Columns[i].ImaginaryComponents[j];	
	}
    }
  return ComplexMatrix(TmpColumns, M1.NbrColumn);
}

// multiply two matrices
//
// M1 = first matrix
// M2 = matrix to multiply to M1
// return value = product of the two matrices

ComplexMatrix operator * (const ComplexMatrix& M1, const ComplexMatrix& M2)
{
  if (M1.NbrColumn != M2.NbrRow)
    return ComplexMatrix();
  ComplexVector* TmpColumns = new ComplexVector [M2.NbrColumn];
  for (int i = 0; i < M2.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M1.NbrRow);
      for (int j = 0; j < M1.NbrRow; ++j)
	{
	  TmpColumns[i].RealComponents[j] = 0.0;
	  TmpColumns[i].ImaginaryComponents[j + 1] = 0.0;
	  for (int k = 0; k < M2.NbrRow; ++k)	
	    {
	      TmpColumns[i].RealComponents[j] += (M1.Columns[k].RealComponents[j] * M2.Columns[i].RealComponents[k] - 
						  M1.Columns[k].ImaginaryComponents[j] * M2.Columns[i].ImaginaryComponents[k]);
	      TmpColumns[i].ImaginaryComponents[j] += (M1.Columns[k].RealComponents[j] * M2.Columns[i].ImaginaryComponents[k] + 
						       M1.Columns[k].ImaginaryComponents[j] * M2.Columns[i].RealComponents[k]);
	    }
	}
    }
  return ComplexMatrix(TmpColumns, M2.NbrColumn);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexMatrix operator * (const ComplexMatrix& M, double x) 
{
  ComplexVector* TmpColumns = new ComplexVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M.NbrRow);
      for (int j = 0; j < M.NbrRow; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M.Columns[i].RealComponents[j] * x;
	  TmpColumns[i].ImaginaryComponents[j] = M.Columns[i].ImaginaryComponents[j] * x;
	}
    }
  return ComplexMatrix(TmpColumns, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexMatrix operator * (double x, const ComplexMatrix& M) 
{
  ComplexVector* TmpColumns = new ComplexVector [M.NbrColumn];
  for (int i = 0; i < M.NbrColumn; ++i)
    {
      TmpColumns[i] = ComplexVector (M.NbrRow);
      for (int j = 0; j < M.NbrRow; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M.Columns[i].RealComponents[j] * x;
	  TmpColumns[i].ImaginaryComponents[j] = M.Columns[i].ImaginaryComponents[j] * x;
	}
    }
  return ComplexMatrix(TmpColumns, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

ComplexMatrix operator / (const ComplexMatrix& M, double x) 
{
  ComplexVector* TmpColumns = new ComplexVector [M.NbrColumn];
  x = 1.0 / x;
  for (int i = 0; i < M.NbrColumn; i++)
    {
      TmpColumns[i] = ComplexVector (M.NbrRow);
      for (int j = 0; j < M.NbrRow; ++j)
	{
	  TmpColumns[i].RealComponents[j] = M.Columns[i].RealComponents[j] * x;
	  TmpColumns[i].ImaginaryComponents[j] = M.Columns[i].ImaginaryComponents[j] * x;
	}
    }
  return ComplexMatrix(TmpColumns, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator += (const ComplexMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] += M.Columns[i];
  return *this;
}

// add two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator += (const RealTriDiagonalSymmetricMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow) || (this->ColumnGarbageFlag == 0))
    return *this;  
  this->Columns[0].RealComponents[0] += M.DiagonalElements[0];
  for (int i = 1; i < this->NbrColumn; i++)
    {
      this->Columns[i].RealComponents[i] += M.DiagonalElements[i];
      this->Columns[i].RealComponents[i - 1] += M.UpperDiagonalElements[i - 1];
      this->Columns[i - 1].RealComponents[i] += M.UpperDiagonalElements[i - 1];
    }
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator -= (const ComplexMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow))
    return *this;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] -= M.Columns[i];
  return *this;
}

// substract two matrices where the right one is a real tridiagonal symmetric matrix
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator -= (const RealTriDiagonalSymmetricMatrix& M) 
{
  if ((this->NbrColumn != M.NbrColumn) || (this->NbrRow != M.NbrRow) || (this->ColumnGarbageFlag == 0))
    return *this;  
  this->Columns[0].RealComponents[0] -= M.DiagonalElements[0];
  for (int i = 1; i < this->NbrColumn; i++)
    {
      this->Columns[i].RealComponents[i] -= M.DiagonalElements[i];
      this->Columns[i].RealComponents[i - 1] -= M.UpperDiagonalElements[i - 1];
      this->Columns[i - 1].RealComponents[i] -= M.UpperDiagonalElements[i - 1];
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator *= (double x) 
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] *= x;
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::operator /= (double x)
{
  x = 1.0 / x;;
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i] *= x;
  return *this;
}

// normalize matrix column vectors
//
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::NormalizeColumns ()
{
  for (int i = 0; i < this->NbrColumn; i++)
    this->Columns[i].Normalize();
  return *this;
}

// orthonormalize matrix column vectors
//
// return value = reference on current matrix

ComplexMatrix& ComplexMatrix::OrthoNormalizeColumns ()
{
  Complex* tmp = new Complex [this->NbrColumn];
  for (int i = 0; i < this->NbrColumn; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  tmp[j] = this->Columns[i] * this->Columns[j];
	}
      for (int j = 0; j < i; j++)
	{
	  this->Columns[i].AddLinearCombination(tmp[j], this->Columns[j]);
	}
      this->Columns[i].Normalize();
    }      
  delete[] tmp;
  return *this;
}
// evaluate matrix determinant (skrewing up matrix elements)
//
// return value = matrix determinant 

Complex ComplexMatrix::Determinant () 
{
  if (this->NbrColumn != this->NbrRow)
    return 0.0;
  Complex TmpDet (1.0);
  int ReducedNbrRow = this->NbrRow - 1;
  Complex Pivot;
  Complex Factor;
  int PivotPos = 0;
  double PivotNorm;
  for (int k = 0; k < ReducedNbrRow; ++k)
    {
      Pivot.Re = this->Columns[k].Re(k);
      Pivot.Im = this->Columns[k].Im(k);
      PivotNorm = (Pivot.Re * Pivot.Re) + (Pivot.Im * Pivot.Im);
      PivotPos = k + 1;
      while ((PivotPos < this->NbrRow) && 
	     (((this->Columns[PivotPos].Re(k) * this->Columns[PivotPos].Re(k)) + (this->Columns[PivotPos].Im(k) * this->Columns[PivotPos].Im(k)))< PivotNorm))
	{
	  ++PivotPos;
	}
      if (PivotPos == this->NbrRow)
	{
	  if (PivotNorm == 0.0)
	    return Complex(0.0);
	}
      else
	{
	  Pivot.Re = this->Columns[PivotPos].Re(k);
	  Pivot.Im = this->Columns[PivotPos].Im(k);
	  ComplexVector TmpColumn3(this->Columns[k]);
	  this->Columns[k] = this->Columns[PivotPos];
	  this->Columns[PivotPos] = TmpColumn3;
	  TmpDet *= -1.0;
	}
      if (PivotNorm == 0.0)
	return Complex(0.0);
      TmpDet *= Pivot;
      Pivot = 1.0 / Pivot;       
      for (int i = k + 1; i < this->NbrRow; ++i)
	{
	  ComplexVector& TmpColumn = this->Columns[i];
	  ComplexVector& TmpColumn2 = this->Columns[k];
	  Factor.Re = ((Pivot.Re * TmpColumn.Re(k)) - (Pivot.Im * TmpColumn.Im(k)));
	  Factor.Im = ((Pivot.Im * TmpColumn.Re(k)) + (Pivot.Re * TmpColumn.Im(k)));
	  for (int j = k + 1; j < this->NbrRow; ++j)
	    {
	      TmpColumn.Re(j) -= ((TmpColumn2.Re(j) * Factor.Re) - (TmpColumn2.Im(j) * Factor.Im)); 
	      TmpColumn.Im(j) -= ((TmpColumn2.Re(j) * Factor.Im) + (TmpColumn2.Im(j) * Factor.Re)); 
	    }
	}
    } 
  Pivot.Re = this->Columns[ReducedNbrRow].Re(ReducedNbrRow);
  Pivot.Im = this->Columns[ReducedNbrRow].Im(ReducedNbrRow);
  TmpDet *= Pivot;
  return TmpDet;
}


// evaluate permanent associated to the (square) matrix using Ryser algorithm
//
// return value = permanent associated to the matrix
                                                                                                                                          
Complex ComplexMatrix::Permanent()
{
  if (this->NbrColumn != this->NbrRow)
    return 0.0;
  Complex Perm;
  double Sign = 1.0;
  if ((this->NbrColumn & 1) == 0)
    Sign = -1.0;
  Complex* Tmp = new Complex [this->NbrColumn];
  Complex Tmp2;
  int Lim = 1 << this->NbrColumn;
  for (int i = 0; i < this->NbrColumn; ++i)
    Tmp[i] = 0.0;
  int GrayCode = 0;
  int ChangedBit;
  int Index;
  for (int k = 1; k < Lim; ++k)
    {
      ChangedBit = (k ^ (k >> 1)) ^ GrayCode;
      GrayCode = k ^ (k >> 1);
      if ((GrayCode & ChangedBit) == 0)
	{
	  Index = 0;
	  while (ChangedBit != 1)
	    {
	      ChangedBit >>= 1;
	      ++Index;
	    }
	  ComplexVector& TmpColumn = this->Columns[Index];
	  for (int i = 0; i < this->NbrColumn; ++i)
	    {
	      Tmp[i].Re -= TmpColumn.RealComponents[i];
	      Tmp[i].Im -= TmpColumn.ImaginaryComponents[i];
	    }
	}
      else
	{
	  Index = 0;
	  while (ChangedBit != 1)
	    {
	      ChangedBit >>= 1;
	      ++Index;
	    }
	  ComplexVector& TmpColumn = this->Columns[Index];
	  for (int i = 0; i < this->NbrColumn; ++i)
	    {
	      Tmp[i].Re += TmpColumn.RealComponents[i];
	      Tmp[i].Im += TmpColumn.ImaginaryComponents[i];
	    }
	}
      Tmp2 = Tmp[0];
      for (int i = 1; i < this->NbrColumn; ++i)
        Tmp2 *= Tmp[i];
      Perm += Sign * Tmp2;
      Sign *= -1.0;
    }
  delete[] Tmp;
  return Perm;
}

// evaluate minor develomment of permanent associated to the (square) matrix using Ryser algorithm
//
// column = index of the column from which permnanent will developped
// minors = reference on an array where minors will be stored
                                                                                                                                          
void ComplexMatrix::PermanentMinorDevelopment(int column, Complex*& minors)
{
  if (this->NbrColumn != this->NbrRow)
    return;
  int ReducedNbrColumn = this->NbrColumn - 1;
  Complex* Tmp = new Complex  [ReducedNbrColumn];
  Complex Tmp2;
  int Lim = 1 << ReducedNbrColumn;
  for (int l = 0; l < this->NbrColumn; ++l)
    {
      minors[l].Re = 0.0;
      minors[l].Im = 0.0;
      double Sign = 1.0;
      if ((ReducedNbrColumn & 1) == 0)
	Sign = -1.0;
      for (int i = 0; i < ReducedNbrColumn; ++i)
	Tmp[i] = 0.0;
      int GrayCode = 0;
      int ChangedBit;
      int Index;
      for (int k = 1; k < Lim; ++k)
	{
	  ChangedBit = (k ^ (k >> 1)) ^ GrayCode;
	  GrayCode = k ^ (k >> 1);
	  if ((GrayCode & ChangedBit) == 0)
	    {
	      Index = 0;
	      while (ChangedBit != 1)
		{
		  ChangedBit >>= 1;
		  ++Index;
		}
	      if (Index >= column)
		++Index;
	      int i = 0;
	      ComplexVector& TmpColumn = this->Columns[Index];
	      for (; i < l; ++i)
		{
		  Tmp[i].Re -= TmpColumn.RealComponents[i];
		  Tmp[i].Im -= TmpColumn.ImaginaryComponents[i];
		}
	      for (; i < ReducedNbrColumn; ++i)
		{
		  Tmp[i].Re -= TmpColumn.RealComponents[i + 1];
		  Tmp[i].Im -= TmpColumn.ImaginaryComponents[i + 1];
		}
	    }
	  else
	    {
	      Index = 0;
	      while (ChangedBit != 1)
		{
		  ChangedBit >>= 1;
		  ++Index;
		}
	      if (Index >= column)
		++Index;
	      int i = 0;
	      ComplexVector& TmpColumn = this->Columns[Index];
	      for (; i < l; ++i)
		{
		  Tmp[i].Re += TmpColumn.RealComponents[i];
		  Tmp[i].Im += TmpColumn.ImaginaryComponents[i];
		}
	      for (; i < ReducedNbrColumn; ++i)
		{
		  Tmp[i].Re += TmpColumn.RealComponents[i + 1];
		  Tmp[i].Im += TmpColumn.ImaginaryComponents[i + 1];
		}
	    }
 	  Tmp2 = Tmp[0];
	  for (int i = 1; i < ReducedNbrColumn; ++i)
	    Tmp2 *= Tmp[i];
	  minors[l] += Sign * Tmp2;
	  Sign *= -1.0;
	}
    }
  delete[] Tmp;
}

// evaluate permanent associated to the (square) matrix using Ryser algorithm using precalculation array (faster)
//
// changeBit = array indicating which bit is changed at the i-th iteration of the Gray code
// changeBitSign = array with 0 if the changed bit is from 1 to 0, +1 either
// return value = permanent associated to the matrix
                                                                                                                                          
Complex ComplexMatrix::FastPermanent(int* changeBit, int* changeBitSign)
{
  Complex Perm;
  double Sign = 1.0;
  if ((this->NbrColumn & 1) == 0)
    Sign = -1.0;
  Complex* Tmp = new Complex [this->NbrColumn];
  Complex Tmp2;
  int Lim = 1 << this->NbrColumn;
  for (int i = 0; i < this->NbrColumn; ++i)
    Tmp[i] = 0.0;
  for (int k = 1; k < Lim; ++k)
    {
      ComplexVector& TmpColumn = this->Columns[changeBit[k]];
      if (changeBitSign[k] == 0)
	{
	  for (int i = 0; i < this->NbrColumn; ++i)
	    {
	      Tmp[i].Re -= TmpColumn.RealComponents[i];
	      Tmp[i].Im -= TmpColumn.ImaginaryComponents[i];	  
	    }
	}
      else
	{
	  for (int i = 0; i < this->NbrColumn; ++i)
	    {
	      Tmp[i].Re += TmpColumn.RealComponents[i];
	      Tmp[i].Im += TmpColumn.ImaginaryComponents[i];	  
	    }
	}
      Tmp2 = Tmp[0];
      for (int i = 1; i < this->NbrColumn; ++i)
        Tmp2 *= Tmp[i];
      Perm += Sign * Tmp2;
      Sign *= -1.0;
    }
  delete[] Tmp;
  return Perm;
}


// evaluate minor develomment of permanent associated to the (square) matrix using Ryser algorithm and precalculation array (faster)
//
// changeBit = array indicating which bit is changed at the i-th iteration of the Gray code
// changeBitSign = array with -1 if the changed bit is from 1 to 0, +1 either
// column = index of the column from which permnanent will developped
// minors = reference on an array where minors will be stored

void ComplexMatrix::FastPermanentMinorDevelopment(int* changeBit, int* changeBitSign, int column, Complex*& minors)
{
  int ReducedNbrColumn = this->NbrColumn - 1;
  Complex* Tmp = new Complex [ReducedNbrColumn];
  Complex Tmp2;
  ComplexVector* TmpColumn;
  int Lim = 1 << ReducedNbrColumn;
  for (int l = 0; l < this->NbrColumn; ++l)
    {
      minors[l].Re = 0.0;
      minors[l].Im = 0.0;
      double Sign = 1.0;
      if ((ReducedNbrColumn & 1) == 0)
	Sign = -1.0;
      for (int i = 0; i < ReducedNbrColumn; ++i)
	Tmp[i] = 0.0;
      for (int k = 1; k < Lim; ++k)
	{
	  if (changeBit[k] < column)
	    {
	      TmpColumn = &(this->Columns[changeBit[k]]);
	    }
	  else
	    {
	      TmpColumn = &(this->Columns[changeBit[k] + 1]);
	    }
	  if (changeBitSign[k] == 0)
	    {
	      int i = 0;
	      for (; i < l; ++i)
		{
		  Tmp[i].Re -= TmpColumn->RealComponents[i];
		  Tmp[i].Im -= TmpColumn->ImaginaryComponents[i];	  
		}
	      for (; i < ReducedNbrColumn; ++i)
		{
		  Tmp[i].Re -= TmpColumn->RealComponents[i + 1];
		  Tmp[i].Im -= TmpColumn->ImaginaryComponents[i + 1];	  
		}
	    }
	  else
	    {
	      int i = 0;
	      for (; i < l; ++i)
		{
		  Tmp[i].Re += TmpColumn->RealComponents[i];
		  Tmp[i].Im += TmpColumn->ImaginaryComponents[i];	  
		}
	      for (; i < ReducedNbrColumn; ++i)
		{
		  Tmp[i].Re += TmpColumn->RealComponents[i + 1];
		  Tmp[i].Im += TmpColumn->ImaginaryComponents[i + 1];	  
		}
	    }
	  Tmp2 = Tmp[0];
	  for (int i = 1; i < ReducedNbrColumn; ++i)
	    Tmp2 *= Tmp[i];
	  minors[l] += Sign * Tmp2;
	  Sign *= -1.0;
	}
    }
  delete[] Tmp;
}
  
// evaluate precalculation array  neede for the fast permanent calculation
//
// changeBit = reference on the array indicating which bit is changed at the i-th iteration of the Gray code
// changeBitSign = reference on array with 0 if the changed bit is from 1 to 0, +1 either
// minor = flag that indicated if precalculation will be used for minor development

void ComplexMatrix::EvaluateFastPermanentPrecalculationArray(int*& changeBit, int*& changeBitSign, bool minor)
{
  int Lim ;
  if (minor == false)
    {
      Lim = 1 << this->NbrColumn;
    }
  else
    {
      Lim = 1 << (this->NbrColumn - 1);
    }
  int GrayCode = 0;
  int ChangedBit;
  int Index;
  changeBit = new int [Lim];
  changeBitSign = new int [Lim];
  for (int k = 1; k < Lim; ++k)
    {
      ChangedBit = (k ^ (k >> 1)) ^ GrayCode;
      GrayCode = k ^ (k >> 1);
      if ((GrayCode & ChangedBit) == 0)
	{
	  changeBitSign[k] = 0;
	}
      else
	{
	  changeBitSign[k] = 1;
	}
      Index = 0;
      while (ChangedBit != 1)
	{
	  ChangedBit >>= 1;
	  ++Index;
	}
      changeBit[k] = Index;
   }
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const ComplexMatrix& P) 
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      for (int j = 0; j < (P.NbrColumn - 1); j ++)
	{
	  Str << P.Columns[j].RealComponents[i];      
	  if (P.Columns[j].ImaginaryComponents[i] < 0.0)
	    Str << P.Columns[j].ImaginaryComponents[i] << "i    ";
	  else
	    if (P.Columns[j].ImaginaryComponents[i] != 0.0)
	      Str << "+" << P.Columns[j].ImaginaryComponents[i] << "i    ";
	    else
	      Str << "    ";
	}
      Str << P.Columns[P.NbrColumn - 1].RealComponents[i];      
      if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] < 0.0)
	Str << P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] << "i";
      else
	if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] != 0.0)
	  Str << "+" << P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] << "i";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexMatrix& P) 
{
  Str << "{";
  for (int i = 0; i < (P.NbrRow - 1); ++i)
    {
      Str << "{";
      for (int j = 0; j < (P.NbrColumn - 1); ++j)
	{
	  Str << P.Columns[j].RealComponents[i];      
	  if (P.Columns[j].ImaginaryComponents[i] < 0.0)
	    Str << P.Columns[j].ImaginaryComponents[i] << "I,";
	  else
	    if (P.Columns[j].ImaginaryComponents[i] != 0.0)
	      Str << "+" << P.Columns[j].ImaginaryComponents[i] << "I,";
	    else
	      Str << ",";
	}
      Str << P.Columns[P.NbrColumn - 1].RealComponents[i];      
      if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] < 0.0)
	Str << P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] << "I";
      else
	if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] != 0.0)
	  Str << "+" << P.Columns[P.NbrColumn - 1].ImaginaryComponents[i] << "I";
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (P.NbrColumn - 1); ++j)
    {
      Str << P.Columns[j].RealComponents[P.NbrRow - 1];      
      if (P.Columns[j].ImaginaryComponents[P.NbrRow - 1] < 0.0)
	Str << P.Columns[j].ImaginaryComponents[P.NbrRow - 1] << "I,";
      else
	if (P.Columns[j].ImaginaryComponents[P.NbrRow - 1] != 0.0)
	  Str << "+" << P.Columns[j].ImaginaryComponents[P.NbrRow - 1] << "I,";
	else
	  Str << ",";
    }
  Str << P.Columns[P.NbrColumn - 1].RealComponents[P.NbrRow - 1];      
  if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[P.NbrRow - 1] < 0.0)
    Str << P.Columns[P.NbrColumn - 1].ImaginaryComponents[P.NbrRow - 1] << "I";
  else
    if (P.Columns[P.NbrColumn - 1].ImaginaryComponents[P.NbrRow - 1] != 0.0)
      Str << "+" << P.Columns[P.NbrColumn - 1].ImaginaryComponents[P.NbrRow - 1] << "I";
  Str << "}}";
  return Str;
}

#endif
