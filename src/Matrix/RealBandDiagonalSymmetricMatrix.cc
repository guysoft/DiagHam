////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of real band-diagonal symmetric matrix                //
//                                                                            //
//                        last modification : 16/03/2005                      //
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


#include "Matrix/RealBandDiagonalSymmetricMatrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"

#include <math.h>


using std::endl;


// default constructor
//

RealBandDiagonalSymmetricMatrix::RealBandDiagonalSymmetricMatrix() 
{
  this->DiagonalElements = 0;
  this->UpperOffDiagonalElements = 0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->NbrBands = 0; 
  this->TrueNbrBands = this->NbrBands;
  this->TrueNbrRow = 0;
  this->TrueNbrColumn = 0;
  this->MatrixType = Matrix::RealElements | Matrix::BandDiagonal | Matrix::Symmetric;
  this->Dummy = 0.0;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

RealBandDiagonalSymmetricMatrix::RealBandDiagonalSymmetricMatrix(int dimension, int nbrBands, bool zero)
{
  this->DiagonalElements = new double [dimension];
  this->UpperOffDiagonalElements = new double* [this->NbrBands];
  this->NbrBands = nbrBands;
  this->TrueNbrBands = this->NbrBands;
  for (int i = 0; i < this->NbrBands; ++i)
    this->UpperOffDiagonalElements[i] = new double [dimension];
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = dimension;
  this->TrueNbrColumn = dimension;
  this->MatrixType = Matrix::RealElements | Matrix::BandDiagonal | Matrix::Symmetric;
  if (zero == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	this->DiagonalElements[i] = 0.0;
      double* Tmp;
      for (int j = 0; j < this->NbrBands; ++j)
	{
	  Tmp = this->UpperOffDiagonalElements[j];
	  for (int i = 0; i < this->NbrRow; i++)
	    Tmp[i] = 0.0;
	}
    }
  this->Dummy = 0.0;
}

// constructor from matrix elements (without duplicating datas)
//
// diagonal = pointer to diagonal element array
// upperOffDiagonal = pointer to the array which contains upper off-diagonal elements (first index is used as row index)
// dimension = matrix dimension
// nbrBands = number of bands in the upper part of the matrix

RealBandDiagonalSymmetricMatrix::RealBandDiagonalSymmetricMatrix(double* diagonal, double** upperDiagonal, int dimension, int nbrBands) 
{
  this->DiagonalElements = diagonal;
  this->UpperOffDiagonalElements = upperDiagonal;
  this->NbrBands = nbrBands;
  this->TrueNbrBands = this->NbrBands;
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = dimension;
  this->TrueNbrColumn = dimension;
  this->MatrixType = Matrix::RealElements | Matrix::BandDiagonal | Matrix::Symmetric;
  this->Dummy = 0.0;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

RealBandDiagonalSymmetricMatrix::RealBandDiagonalSymmetricMatrix(const RealBandDiagonalSymmetricMatrix& M) 
{  
  this->DiagonalElements = M.DiagonalElements;
  this->UpperOffDiagonalElements = M.UpperOffDiagonalElements;
  this->NbrBands = M.NbrBands;
  this->TrueNbrBands = M.NbrBands;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::BandDiagonal | Matrix::Symmetric;
  this->Dummy = 0.0;
}

// destructor
//

RealBandDiagonalSymmetricMatrix::~RealBandDiagonalSymmetricMatrix() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->TrueNbrBands; ++i)
	delete[] this->UpperOffDiagonalElements[i];
      delete[] this->UpperOffDiagonalElements;
      delete[] this->DiagonalElements;
    }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

RealBandDiagonalSymmetricMatrix& RealBandDiagonalSymmetricMatrix::operator = (const RealBandDiagonalSymmetricMatrix& M) 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->TrueNbrBands; ++i)
	delete[] this->UpperOffDiagonalElements[i];
      delete[] this->UpperOffDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = M.DiagonalElements;
  this->UpperOffDiagonalElements = M.UpperOffDiagonalElements;
  this->NbrBands = M.NbrBands;
  this->TrueNbrBands = M.NbrBands;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->MatrixType = Matrix::RealElements | Matrix::BandDiagonal | Matrix::Symmetric;
  this->Dummy = 0.0;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* RealBandDiagonalSymmetricMatrix::Clone ()
{
  return ((Matrix*) new RealBandDiagonalSymmetricMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealBandDiagonalSymmetricMatrix::SetMatrixElement(int i, int j, double x)
{
  if ((i == j) && (i < this->NbrRow))
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      this->UpperOffDiagonalElements[j - 1][i] = x;
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      this->UpperOffDiagonalElements[i - 1][j] = x;
	    }	
	}
    }    
}

// return refernce on real part of a given matrix element
//
// i = line position
// j = column position
// return value = reference on real part

double& RealBandDiagonalSymmetricMatrix::operator () (int i, int j)
{
  if ((i == j) && (i < this->NbrRow))
    {
      return this->DiagonalElements[i];
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      return this->UpperOffDiagonalElements[j - 1][i];
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      return this->UpperOffDiagonalElements[i - 1][j];
	    }	
	}
    }    
  return this->Dummy;
}

// get a matrix element 
// 
// i = Row number
// j = Column number
// return value = matrix element M_(i,j)

double RealBandDiagonalSymmetricMatrix::GetElement(int i, int j)
{
  if ((i == j) && (i < this->NbrRow))
    {
      return this->DiagonalElements[i];
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      return this->UpperOffDiagonalElements[j - 1][i];
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      return this->UpperOffDiagonalElements[i - 1][j];
	    }	
	}
    }    
  return 0.0;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void RealBandDiagonalSymmetricMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
}

// access to i-th diagonal element
// 
// i = position 
// return value = reference on i-th diagonal element

double& RealBandDiagonalSymmetricMatrix::DiagonalElement(int i)
{
  return this->DiagonalElements[i];
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void RealBandDiagonalSymmetricMatrix::AddToMatrixElement(int i, int j, double x)
{
  if ((i == j) && (i < this->NbrRow))
    {
      this->DiagonalElements[i] = x;
    }
  else
    {
      if (j > i)
	{
	  j -= i;
	  if ((j <= this->NbrBands) && (i < (this->NbrRow - 1)))
	    {
	      this->UpperOffDiagonalElements[j - 1][i] += x;
	    }
	}
      else
	{
	  i -= j;
	  if ((i <= this->NbrBands) && (j < (this->NbrRow - 1)))
	    {
	      this->UpperOffDiagonalElements[i - 1][j] += x;
	    }	
	}
    }    
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void RealBandDiagonalSymmetricMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealBandDiagonalSymmetricMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      return;
    }
  double* TmpDiag = new double [nbrRow];
  double** TmpUpperDiag = new double* [this->TrueNbrBands];
  for (int i = 0; i < this->TrueNbrBands; ++i)
    TmpUpperDiag[i] = new double [nbrRow];
  if (this->Flag.Used() == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	TmpDiag[i] = this->DiagonalElements[i];
      for (int j = 0; j < this->NbrBands; ++j)
	for (int i = 0; i < this->NbrRow; ++i)
	  TmpUpperDiag[j][i] = this->UpperOffDiagonalElements[j][i];
    }
   if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int j = 0; j < this->TrueNbrBands; ++j)
	delete[] this->UpperOffDiagonalElements[j];
      delete[] this->UpperOffDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = TmpDiag;
  this->UpperOffDiagonalElements = TmpUpperDiag;
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = nbrRow;
  this->TrueNbrColumn = nbrColumn;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
  return;
}

// Resize matrix and change the number of bands
//
// nbrRow = new number of rows
// nbrColumn = new number of columns
// nbrBands = new number of bands

void RealBandDiagonalSymmetricMatrix::Resize (int nbrRow, int nbrColumn, int nbrBands)
{
  if (nbrRow != nbrColumn)
    return;
  if ((nbrRow <= this->TrueNbrRow) && (nbrBands <= this->TrueNbrBands))
    {
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      this->NbrBands = nbrBands;
      return;
    }
  double* TmpDiag = new double [nbrRow];
  double** TmpUpperDiag = new double* [nbrBands];
  for (int i = 0; i < nbrBands; ++i)
    TmpUpperDiag[i] = new double [nbrRow];
  if (this->Flag.Used() == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	TmpDiag[i] = this->DiagonalElements[i];
      for (int j = 0; j < this->NbrBands; ++j)
	for (int i = 0; i < this->NbrRow; ++i)
	  TmpUpperDiag[j][i] = this->UpperOffDiagonalElements[j][i];
    }
   if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int j = 0; j < this->TrueNbrBands; ++j)
	delete[] this->UpperOffDiagonalElements[j];
      delete[] this->UpperOffDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = TmpDiag;
  this->UpperOffDiagonalElements = TmpUpperDiag;
  this->NbrBands = nbrBands;
  this->TrueNbrBands = nbrBands;
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = nbrRow;
  this->TrueNbrColumn = nbrColumn;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
  return;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void RealBandDiagonalSymmetricMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      for (int i = this->NbrRow; i < nbrRow; i++)
	{
	  this->DiagonalElements[i] = 0.0;
	}
      for (int j = 0; j < this->NbrBands; ++j)
	for (int i = this->NbrRow; i < nbrRow; ++i)
	  this->UpperOffDiagonalElements[j][i] = 0.0;	  
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      return;
    }
  double* TmpDiag = new double [nbrRow];
  double** TmpUpperDiag = new double* [this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    TmpUpperDiag[i] = new double [nbrRow];
  if (this->Flag.Used() == true)
    {
      int i = 0;
      for (; i < this->NbrRow; i++)
	TmpDiag[i] = this->DiagonalElements[i];
      for (; i < nbrRow; i++)
	TmpDiag[i] = 0.0;
      for (int j = 0; j < this->NbrBands; ++j)
	for (i = 0; i < this->NbrRow; ++i)
	  TmpUpperDiag[j][i] = this->UpperOffDiagonalElements[j][i];
      for (int j = 0; j < this->NbrBands; ++j)
	for (i = nbrRow; i < this->NbrRow; ++i)
	  TmpUpperDiag[j][i] = 0.0;
    }
   if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int j = 0; j < this->NbrBands; ++j)
	delete[] this->UpperOffDiagonalElements[j];
      delete[] this->UpperOffDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = TmpDiag;
  this->UpperOffDiagonalElements = TmpUpperDiag;
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = nbrRow;
  this->TrueNbrColumn = nbrColumn;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
  return;
}

// copy matrix
//
// M = matrix to copy
// return value = refence on current matrix

RealBandDiagonalSymmetricMatrix& RealBandDiagonalSymmetricMatrix::Copy (RealBandDiagonalSymmetricMatrix& M)
{
  if (this->NbrBands != M.NbrBands)
    this->Resize(M.NbrRow, M.NbrColumn, M.NbrBands);
  else
    if (this->NbrRow != M.NbrRow)
      this->Resize(M.NbrRow, M.NbrColumn);
  for (int i = 0; i < M.NbrColumn; i++)
    this->DiagonalElements[i] = M.DiagonalElements[i];
  double* Tmp1;
  double* Tmp2;
  for (int j = 0; j < this->NbrBands; ++j)
    {
      Tmp1 = this->UpperOffDiagonalElements[j];
      Tmp2 = M.UpperOffDiagonalElements[j];
      for (int i = 0; i < this->NbrRow; i++)
	Tmp1[i] = Tmp2[i];
    }
  return *this;
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

RealBandDiagonalSymmetricMatrix operator + (const RealBandDiagonalSymmetricMatrix& M1, const RealBandDiagonalSymmetricMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return RealBandDiagonalSymmetricMatrix();
  int MaxNbrBands = M1.NbrBands;
  int MinNbrBands = M2.NbrBands;
  if (M1.NbrBands < M2.NbrBands)
    {
      MaxNbrBands = M2.NbrBands;
      MinNbrBands = M1.NbrBands;
    }
  double* Diagonal = new double [M1.NbrRow];
  double** UpperDiagonal = new double* [MaxNbrBands];
  for (int i = 0; i < MaxNbrBands; ++i)
    UpperDiagonal[i] = new double [M1.NbrRow];
  for (int i = 0; i < M1.NbrRow; i++)
    Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
  double* Tmp1;
  double* Tmp2;
  double* Tmp3;
  for (int j = 0; j < MinNbrBands; ++j)
    {      
      Tmp1 = M1.UpperOffDiagonalElements[j];
      Tmp2 = M2.UpperOffDiagonalElements[j];
      Tmp3 = UpperDiagonal[j];
      for (int i = 0; i < M1.NbrRow; i++)
	Tmp3[i] = Tmp1[i] + Tmp2[i];
    }
  for (int j = MinNbrBands; j < MaxNbrBands; ++j)
    {      
      if (M1.NbrBands == MaxNbrBands)
	Tmp1 = M1.UpperOffDiagonalElements[j];
      else
	Tmp1 = M2.UpperOffDiagonalElements[j];
      Tmp3 = UpperDiagonal[j];
      for (int i = 0; i < M1.NbrRow; i++)
	Tmp3[i] = Tmp1[i];
    }
  return RealBandDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M1.NbrRow, MaxNbrBands);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

RealBandDiagonalSymmetricMatrix operator - (const RealBandDiagonalSymmetricMatrix& M1, const RealBandDiagonalSymmetricMatrix& M2) 
{
  if (M1.NbrRow != M2.NbrRow)
    return RealBandDiagonalSymmetricMatrix();
  int MaxNbrBands = M1.NbrBands;
  int MinNbrBands = M2.NbrBands;
  if (M1.NbrBands < M2.NbrBands)
    {
      MaxNbrBands = M2.NbrBands;
      MinNbrBands = M1.NbrBands;
    }
  double* Diagonal = new double [M1.NbrRow];
  double** UpperDiagonal = new double* [MaxNbrBands];
  for (int i = 0; i < MaxNbrBands; ++i)
    UpperDiagonal[i] = new double [M1.NbrRow];
  for (int i = 0; i < M1.NbrRow; i++)
    Diagonal[i] = M1.DiagonalElements[i] - M2.DiagonalElements[i];
  double* Tmp1;
  double* Tmp2;
  double* Tmp3;
  for (int j = 0; j < MinNbrBands; ++j)
    {      
      Tmp1 = M1.UpperOffDiagonalElements[j];
      Tmp2 = M2.UpperOffDiagonalElements[j];
      Tmp3 = UpperDiagonal[j];
      for (int i = 0; i < M1.NbrRow; i++)
	Tmp3[i] = Tmp1[i] + Tmp2[i];
    }
  if (M1.NbrBands == MaxNbrBands)
    for (int j = MinNbrBands; j < MaxNbrBands; ++j)
      {      
	Tmp1 = M1.UpperOffDiagonalElements[j];
	Tmp3 = UpperDiagonal[j];
	for (int i = 0; i < M1.NbrRow; i++)
	  Tmp3[i] = Tmp1[i];
      }
  else
    for (int j = MinNbrBands; j < MaxNbrBands; ++j)
      {      
	Tmp1 = M2.UpperOffDiagonalElements[j];
	Tmp3 = UpperDiagonal[j];
	for (int i = 0; i < M1.NbrRow; i++)
	  Tmp3[i] = -Tmp1[i];
      }
  return RealBandDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M1.NbrRow, MaxNbrBands);
}

// multiply a matrix with a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealBandDiagonalSymmetricMatrix operator * (const RealBandDiagonalSymmetricMatrix& M, double x) 
{
  double* Diagonal = new double [M.NbrRow];
  double** UpperDiagonal = new double* [M.NbrBands];
  for (int i = 0; i < M.NbrRow; ++i)
    Diagonal[i] = M.DiagonalElements[i] * x;      
  for (int i = 0; i < M.NbrBands; ++i)
    {
      UpperDiagonal[i] = new double [M.NbrRow];
      for (int j = 0; j < M.NbrRow; ++j)
	UpperDiagonal[i][j] = M.UpperOffDiagonalElements[i][j] * x;      	
    }
  return RealBandDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M.NbrRow, M.NbrBands);
}

// multiply a matrix with a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealBandDiagonalSymmetricMatrix operator * (double x, const RealBandDiagonalSymmetricMatrix& M) 
{
  double* Diagonal = new double [M.NbrRow];
  double** UpperDiagonal = new double* [M.NbrBands];
  for (int i = 0; i < M.NbrRow; ++i)
    Diagonal[i] = M.DiagonalElements[i] * x;      
  for (int i = 0; i < M.NbrBands; ++i)
    {
      UpperDiagonal[i] = new double [M.NbrRow];
      for (int j = 0; j < M.NbrRow; ++j)
	UpperDiagonal[i][j] = M.UpperOffDiagonalElements[i][j] * x;      	
    }
  return RealBandDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M.NbrRow, M.NbrBands);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

RealBandDiagonalSymmetricMatrix operator / (const RealBandDiagonalSymmetricMatrix& M, double x) 
{
  x = 1.0 / x;
  double* Diagonal = new double [M.NbrRow];
  double** UpperDiagonal = new double* [M.NbrBands];
  for (int i = 0; i < M.NbrRow; ++i)
    Diagonal[i] = M.DiagonalElements[i] * x;      
  for (int i = 0; i < M.NbrBands; ++i)
    {
      UpperDiagonal[i] = new double [M.NbrRow];
      for (int j = 0; j < M.NbrRow; ++j)
	UpperDiagonal[i][j] = M.UpperOffDiagonalElements[i][j] * x;      	
    }
  return RealBandDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M.NbrRow, M.NbrBands);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

RealBandDiagonalSymmetricMatrix& RealBandDiagonalSymmetricMatrix::operator += (const RealBandDiagonalSymmetricMatrix& M) 
{
  if ((this->NbrRow != M.NbrRow) || (this->NbrBands < M.NbrBands))
    return *this;
  for (int i = 0; i < this->NbrRow; ++i)
    this->DiagonalElements[i] += M.DiagonalElements[i];
  for (int j = 0; j < M.NbrBands; ++j)
    for (int i = 0; i < this->NbrRow; ++i)
      this->UpperOffDiagonalElements[j][i] += M.UpperOffDiagonalElements[j][i];      
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

RealBandDiagonalSymmetricMatrix& RealBandDiagonalSymmetricMatrix::operator -= (const RealBandDiagonalSymmetricMatrix& M) 
{
  if ((this->NbrRow != M.NbrRow) || (this->NbrBands < M.NbrBands))
    return *this;
  for (int i = 0; i < this->NbrRow; ++i)
    this->DiagonalElements[i] -= M.DiagonalElements[i];
  for (int j = 0; j < M.NbrBands; ++j)
    for (int i = 0; i < this->NbrRow; ++i)
      this->UpperOffDiagonalElements[j][i] -= M.UpperOffDiagonalElements[j][i];      
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealBandDiagonalSymmetricMatrix& RealBandDiagonalSymmetricMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  for (int i = 0; i < this->NbrRow; i++)
    this->DiagonalElements[i] *= x;
  for (int j = 0; j < this->NbrBands; ++j)
    for (int i = 0; i < this->NbrRow; i++)    
      this->UpperOffDiagonalElements[i][j] *= x;      
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

RealBandDiagonalSymmetricMatrix& RealBandDiagonalSymmetricMatrix::operator /= (double x) 
{  
  if (this->NbrRow == 0)
    return *this;
  x = 1.0 / x;
  for (int i = 0; i < this->NbrRow; i++)
    this->DiagonalElements[i] *= x;
  for (int j = 0; j < this->NbrBands; ++j)
    for (int i = 0; i < this->NbrRow; i++)    
      this->UpperOffDiagonalElements[i][j] *= x;      
  return *this;
}

// evaluate matrix trace
//
// return value = matrix trace 

double RealBandDiagonalSymmetricMatrix::Tr () 
{
  if (this->NbrRow == 0)
    return 0.0;
  double x = this->DiagonalElements[0];
  for (int i = 1; i < this->NbrRow; i++)
    {
      x += this->DiagonalElements[i];
    }
  return x;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double RealBandDiagonalSymmetricMatrix::Det () 
{
  return 1.0;
}

// Tridiagonalize a real band symmetric matrix using Rutishauer-Schwarz (modifying current matrix)
//
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// err = absolute error on matrix element
// return value = reference on real tridiagonal symmetric matrix

RealTriDiagonalSymmetricMatrix& RealBandDiagonalSymmetricMatrix::Tridiagonalize (RealTriDiagonalSymmetricMatrix& M, double err)
{
  if (M.NbrRow != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  int ReducedNbrRow = this->NbrRow - 2;
  double TmpNorm;
  double Cosinus;
  double Sinus;
  double Tmp;
  double Tmp2;
  double FillInElement;
  int MinJ;
  int Pos;
  int Pos2;
  int Max;
  int GivenColumnPosition;
  int GivenRowPosition;
  double SquareErr = err * err;
  if (SquareErr < MACHINE_PRECISION)
    SquareErr = MACHINE_PRECISION;

  for (int i = 0; i < ReducedNbrRow; ++i)
    {
      MinJ = ReducedNbrRow - i;
      if (MinJ >= this->NbrBands)
	MinJ = this->NbrBands - 1;
      for (int j = MinJ; j >= 1; --j)
	{
	  GivenRowPosition = i;
	  GivenColumnPosition = j;
	  FillInElement = this->UpperOffDiagonalElements[GivenColumnPosition][GivenRowPosition]; 
	  Cosinus = this->UpperOffDiagonalElements[GivenColumnPosition - 1][GivenRowPosition];
	  Sinus = -FillInElement;
	  TmpNorm =  ((Cosinus * Cosinus) + (Sinus * Sinus));
	  while ((GivenRowPosition < this->NbrRow) && (((Sinus * Sinus) > (TmpNorm * SquareErr)) && (TmpNorm > SquareErr)))
	    {
	      // zeroing outmost element of the i-th line using Given rotation and apllying Given rotations to chase out the fill-in element produced by the previous Given rotation
	      TmpNorm = sqrt (TmpNorm);
	      this->UpperOffDiagonalElements[GivenColumnPosition - 1][GivenRowPosition] = TmpNorm;
	      TmpNorm = 1.0 / TmpNorm;
	      Cosinus *= TmpNorm;
	      Sinus *= TmpNorm;
	      Pos = GivenRowPosition + 1;
	      Pos2 = GivenColumnPosition - 1;
	      while (Pos2 > 0)
		{
		  Tmp = this->UpperOffDiagonalElements[Pos2][Pos];
		  this->UpperOffDiagonalElements[Pos2][Pos] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2][Pos] += Sinus * this->UpperOffDiagonalElements[Pos2 - 1][Pos];
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos] -= Sinus * Tmp;
		  --Pos2;
		  ++Pos;
		}
	      
	      Tmp = this->UpperOffDiagonalElements[0][Pos];
	      Tmp2 = this->DiagonalElements[Pos];
	      this->DiagonalElements[Pos] *= Cosinus * Cosinus;
	      this->DiagonalElements[Pos] += Sinus * ((Sinus * this->DiagonalElements[Pos + 1]) - (2.0 * Cosinus * Tmp));
	      this->UpperOffDiagonalElements[0][Pos] *= (Cosinus * Cosinus - Sinus * Sinus);
	      this->UpperOffDiagonalElements[0][Pos] += Cosinus * Sinus * (Tmp - this->DiagonalElements[Pos + 1]);
	      this->DiagonalElements[Pos + 1] *= Cosinus * Cosinus;
	      this->DiagonalElements[Pos + 1] += Sinus * ((Sinus * Tmp2) + (2.0 * Cosinus * Tmp));
	      
	      ++Pos;
	      Pos2 = 1;
	      Max = ReducedNbrRow - GivenRowPosition + 1;
	      if (Max > this->NbrBands)
		Max = this->NbrBands;
	      while (Pos2 < Max)
		{
		  Tmp = this->UpperOffDiagonalElements[Pos2][Pos];
		  this->UpperOffDiagonalElements[Pos2][Pos] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2][Pos] -= Sinus * this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1];
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1] += Sinus * Tmp;
		  ++Pos2;
		}
	      
	      
	      FillInElement = -this->UpperOffDiagonalElements[Max - 1][Pos] * Sinus;
	      this->UpperOffDiagonalElements[Max - 1][Pos] *= Cosinus;


	      GivenRowPosition += GivenColumnPosition;
	      GivenColumnPosition = this->NbrBands - 1;
	      Cosinus = this->UpperOffDiagonalElements[GivenColumnPosition - 1][GivenRowPosition];
	      Sinus = -FillInElement;
	      TmpNorm =  ((Cosinus * Cosinus) + (Sinus * Sinus));
	    }
	}
    }

  for (int i = 0; i < this->NbrRow; ++i)
    M.DiagonalElements[i] = this->DiagonalElements[i];
  ++ReducedNbrRow;
  double* TmpColumn = this->UpperOffDiagonalElements[0];
  for (int i = 0; i < ReducedNbrRow; ++i)
    M.UpperDiagonalElements[i] = TmpColumn[i];

  return M;
}

// Tridiagonalize a real band symmetric matrix using Rutishauer-Schwarz and evaluate transformation matrix  (modifying current matrix)
//
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// err = absolute error on matrix element
// Q = matrix where transformation matrix has to be stored
// return value = reference on real tridiagonal symmetric matrix

RealTriDiagonalSymmetricMatrix& RealBandDiagonalSymmetricMatrix::Tridiagonalize (RealTriDiagonalSymmetricMatrix& M, double err, RealMatrix& Q)
{
  if (M.NbrRow != this->NbrRow)
    M.Resize(this->NbrRow, this->NbrColumn);
  int ReducedNbrRow = this->NbrRow - 2;
  double TmpNorm;
  double Cosinus;
  double Sinus;
  double Tmp;
  double Tmp2;
  double FillInElement;
  int MinJ;
  int Pos;
  int Pos2;
  int Max;
  int GivenColumnPosition;
  int GivenRowPosition;
  double SquareErr = err * err;
  if (SquareErr < MACHINE_PRECISION)
    SquareErr = MACHINE_PRECISION;

  for (int i = 0; i < ReducedNbrRow; ++i)
    {
      MinJ = ReducedNbrRow - i;
      if (MinJ >= this->NbrBands)
	MinJ = this->NbrBands - 1;
      for (int j = MinJ; j >= 1; --j)
	{
	  GivenRowPosition = i;
	  GivenColumnPosition = j;
	  FillInElement = this->UpperOffDiagonalElements[GivenColumnPosition][GivenRowPosition]; 
	  Cosinus = this->UpperOffDiagonalElements[GivenColumnPosition - 1][GivenRowPosition];
	  Sinus = -FillInElement;
	  TmpNorm =  ((Cosinus * Cosinus) + (Sinus * Sinus));
	  while ((GivenRowPosition < this->NbrRow) && (((Sinus * Sinus) > (TmpNorm * SquareErr)) && (TmpNorm > SquareErr)))
	    {
	      // zeroing outmost element of the i-th line using Given rotation and apllying Given rotations to chase out the fill-in element produced by the previous Given rotation
	      TmpNorm = sqrt (TmpNorm);
	      this->UpperOffDiagonalElements[GivenColumnPosition - 1][GivenRowPosition] = TmpNorm;
	      TmpNorm = 1.0 / TmpNorm;
	      Cosinus *= TmpNorm;
	      Sinus *= TmpNorm;
	      Pos = GivenRowPosition + 1;
	      Pos2 = GivenColumnPosition - 1;
	      while (Pos2 > 0)
		{
		  Tmp = this->UpperOffDiagonalElements[Pos2][Pos];
		  this->UpperOffDiagonalElements[Pos2][Pos] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2][Pos] += Sinus * this->UpperOffDiagonalElements[Pos2 - 1][Pos];
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos] -= Sinus * Tmp;
		  --Pos2;
		  ++Pos;
		}
	      
	      Tmp = this->UpperOffDiagonalElements[0][Pos];
	      Tmp2 = this->DiagonalElements[Pos];
	      this->DiagonalElements[Pos] *= Cosinus * Cosinus;
	      this->DiagonalElements[Pos] += Sinus * ((Sinus * this->DiagonalElements[Pos + 1]) - (2.0 * Cosinus * Tmp));
	      this->UpperOffDiagonalElements[0][Pos] *= (Cosinus * Cosinus - Sinus * Sinus);
	      this->UpperOffDiagonalElements[0][Pos] += Cosinus * Sinus * (Tmp - this->DiagonalElements[Pos + 1]);
	      this->DiagonalElements[Pos + 1] *= Cosinus * Cosinus;
	      this->DiagonalElements[Pos + 1] += Sinus * ((Sinus * Tmp2) + (2.0 * Cosinus * Tmp));
	      
	      ++Pos;
	      Pos2 = 1;
	      Max = ReducedNbrRow - GivenRowPosition + 1;
	      if (Max > this->NbrBands)
		Max = this->NbrBands;
	      while (Pos2 < Max)
		{
		  Tmp = this->UpperOffDiagonalElements[Pos2][Pos];
		  this->UpperOffDiagonalElements[Pos2][Pos] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2][Pos] -= Sinus * this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1];
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1] *= Cosinus;
		  this->UpperOffDiagonalElements[Pos2 - 1][Pos + 1] += Sinus * Tmp;
		  ++Pos2;
		}
	      
	      
	      RealVector& Vector1 =  Q[GivenRowPosition + GivenColumnPosition];
	      RealVector& Vector2 =  Q[GivenRowPosition + GivenColumnPosition + 1];
	      for (int k = 0; k < this->NbrRow; ++k)
		{
		  Tmp = Vector2[k];
		  Vector2[k] *= Cosinus;
		  Vector2[k] += Sinus * Vector1[k];
		  Vector1[k] *= Cosinus;
		  Vector1[k] += Sinus * Tmp;
		}
	      
	      FillInElement = -this->UpperOffDiagonalElements[Max - 1][Pos] * Sinus;
	      this->UpperOffDiagonalElements[Max - 1][Pos] *= Cosinus;


	      GivenRowPosition += GivenColumnPosition;
	      GivenColumnPosition = this->NbrBands - 1;
	      Cosinus = this->UpperOffDiagonalElements[GivenColumnPosition - 1][GivenRowPosition];
	      Sinus = -FillInElement;
	      TmpNorm =  ((Cosinus * Cosinus) + (Sinus * Sinus));
	    }
	}
    }

  for (int i = 0; i < this->NbrRow; ++i)
    M.DiagonalElements[i] = this->DiagonalElements[i];
  ++ReducedNbrRow;
  double* TmpColumn = this->UpperOffDiagonalElements[0];
  for (int i = 0; i < ReducedNbrRow; ++i)
    M.UpperDiagonalElements[i] = TmpColumn[i];

  return M;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const RealBandDiagonalSymmetricMatrix& P)
{
  for (int i = 0; i < P.NbrRow; ++i)
    {
      int j = 0;
      for (; j < (i - P.NbrBands); ++j)
	Str << "0    ";
      for (; j < i; ++j)
	{
	  Str << P.UpperOffDiagonalElements[i - j - 1][j] << "    ";
	}
      Str << P.DiagonalElements[i] << "    ";
      ++j;
      for (; ((j < (P.NbrBands + i)) && (j < P.NbrColumn)); ++j)
	{
	  Str << P.UpperOffDiagonalElements[j - i - 1][i] << "    ";
	}
      for (; j < P.NbrColumn; j++)
	Str << "0    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, const RealBandDiagonalSymmetricMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; i++)
    {
      Str << "{";
      int j = 0;
      for (; j < (i - P.NbrBands); ++j)
	Str << "0,";
      for (; j < i; ++j)
	Str << P.UpperOffDiagonalElements[i - j - 1][j] << ",";
      Str << P.DiagonalElements[i] << ",";
      j++;

      for (; ((j < (P.NbrBands + i)) && (j < P.NbrColumn)); ++j)
	{
	  Str << P.UpperOffDiagonalElements[j - i - 1][i];
	  if (j != (P.NbrColumn - 1))
	    Str << ",";
	}
      for (; j < (P.NbrColumn - 1); j++)
	Str << "0,";
      if (j == (P.NbrColumn - 1))
	Str << "0";
      Str << "}";
      if (i != (P.NbrRow - 1))
	Str << ",";
    }
  Str << "}";
  return Str;
}

#endif
