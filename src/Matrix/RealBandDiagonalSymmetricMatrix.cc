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
#include "Matrix/RealUpperTriangularMatrix.h"
#include "Matrix/RealLowerTriangularMatrix.h"
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
  for (int i = 0; i < this->NbrBands; ++i)
    this->UpperOffDiagonalElements[i] = new double [dimension]
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = dimension;
  this->TrueNbrColumn = dimension;
  this->MatrixType = Matrix::RealElements | Matrix::BandDiagonal | Matrix::Symmetric;
  if (zero == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	{
	  this->DiagonalElements[i] = 0.0;
	  Tmp = this->UpperOffDiagonalElements[i];
	}
      double* Tmp;
      for (int j = 0; j < this->NbrBands; ++i)
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
      for (int i = 0; i < this->NbrBands; ++i)
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
      for (int i = 0; i < this->NbrBands; ++i)
	delete[] this->UpperOffDiagonalElements[i];
      delete[] this->UpperOffDiagonalElements;
      delete[] this->DiagonalElements;
    }
  this->DiagonalElements = M.DiagonalElements;
  this->UpperOffDiagonalElements = M.UpperOffDiagonalElements;
  this->NbrBands = M.NbrBands;
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
      j -= i;
      if ((j == 1) && (i < (this->NbrRow - 1)))
	{
	  this->UpperOffDiagonalElements[i] = x;
	}
      else
	if ((j == -1) && (i < this->NbrRow))
	  {
	    this->UpperOffDiagonalElements[i - 1] = x;
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
	  if ((j <= this->NbrBand) && (i < (this->NbrRow - 1)))
	    {
	      return this->UpperOffDiagonalElements[j - 1][i];
	    }
	}
      else
	{
	  i -= j;
	  if ((i > 0) && (i <= this->NbrBand) && (j > 0))
	    {
	      return this->UpperOffDiagonalElements[i - 1][j - 1];
	    }	
	}
    }    
  return this->Dummy;
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
      this->DiagonalElements[i] += x;
    }
  else
    {
      j -= i;
      if ((j == 1) && (i < (this->NbrRow - 1)))
	{
	  this->UpperOffDiagonalElements[i] += x;
	}
      else
	if ((j == -1) && (i < this->NbrRow))
	  {
	    this->UpperOffDiagonalElements[i - 1] += x;
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
  double* TmpUpperDiag = new double [nbrRow];
  if (this->Flag.Used() == true)
    {
      for (int i = 0; i < this->NbrRow; i++)
	{
	  TmpDiag[i] = this->DiagonalElements[i];
	  TmpUpperDiag[i] = this->UpperOffDiagonalElements[i]; 
	}
    }
   if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
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
	  this->UpperOffDiagonalElements[i] = 0.0;
	}
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      return;
    }
  double* TmpDiag = new double [nbrRow];
  double* TmpUpperDiag = new double [nbrRow];
  if (this->Flag.Used() == true)
    {
      int i = 0;
      for (; i < this->NbrRow; i++)
	{
	  TmpDiag[i] = this->DiagonalElements[i];
	  TmpUpperDiag[i] = this->UpperOffDiagonalElements[i]; 
	}
      for (; i < nbrRow; i++)
	{
	  TmpDiag[i] = 0.0;
	  TmpUpperDiag[i] = 0.0; 
	}
    }
   if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
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
  if (this->NbrRow != M.NbrRow)
    this->Resize(M.NbrRow, M.NbrColumn);
  for (int i = 0; i < M.NbrColumn; i++)
    {
      this->DiagonalElements[i] = M.DiagonalElements[i];
    }
  double* Tmp1;
  double* Tmp2;
  for (int j = 0; j < this->NbrBands; ++i)
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
  int ReducedNbr = M1.NbrRow - 1;
  double* Diagonal = new double [M1.NbrRow];
  double** UpperDiagonal = new double [M1.NbrBands];
  for (int i = 0; i < M1.NbrRow; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] + M2.DiagonalElements[i];
      UpperDiagonal[i] = M1.UpperOffDiagonalElements[i] + M2.UpperOffDiagonalElements[i];      
    }
  double* Tmp1;
  double* Tmp2;
  for (int j = 0; j < this->NbrBands; ++i)
    {
      
      Tmp1 = M1.UpperOffDiagonalElements[j];
      Tmp2 = M2.UpperOffDiagonalElements[j];
      for (int i = 0; i < M1.NbrRow; i++)
	Tmp1[i] = Tmp2[i];
    }
  Diagonal[ReducedNbr] = M1.DiagonalElements[ReducedNbr] + M2.DiagonalElements[ReducedNbr];
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
  int ReducedNbr = M1.NbrRow - 1;
  double* Diagonal = new double [M1.NbrRow];
  double* UpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; i++)
    {
      Diagonal[i] = M1.DiagonalElements[i] - M2.DiagonalElements[i];
      UpperDiagonal[i] = M1.UpperOffDiagonalElements[i] - M2.UpperOffDiagonalElements[i];      
    }
  Diagonal[ReducedNbr] = M1.DiagonalElements[ReducedNbr] - M2.DiagonalElements[ReducedNbr];
  return RealBandDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M1.NbrRow);
}

// multiply a matrix with a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealBandDiagonalSymmetricMatrix operator * (const RealBandDiagonalSymmetricMatrix& M, double x) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* Diagonal = new double [M.NbrRow];
  double* UpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
      UpperDiagonal[i] = M.UpperOffDiagonalElements[i] * x;      
    }
  Diagonal[ReducedNbr] = M.DiagonalElements[ReducedNbr] * x;
  return RealBandDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M.NbrRow);
}

// multiply a matrix with a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

RealBandDiagonalSymmetricMatrix operator * (double x, const RealBandDiagonalSymmetricMatrix& M) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* Diagonal = new double [M.NbrRow];
  double* UpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] * x;
      UpperDiagonal[i] = M.UpperOffDiagonalElements[i] * x;      
    }
  Diagonal[ReducedNbr] = M.DiagonalElements[ReducedNbr] * x;
  return RealBandDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

RealBandDiagonalSymmetricMatrix operator / (const RealBandDiagonalSymmetricMatrix& M, double x) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* Diagonal = new double [M.NbrRow];
  double* UpperDiagonal = new double [ReducedNbr];
  for (int i = 0; i < ReducedNbr; i++)
    {
      Diagonal[i] = M.DiagonalElements[i] / x;
      UpperDiagonal[i] = M.UpperOffDiagonalElements[i] / x;      
    }
  Diagonal[ReducedNbr] = M.DiagonalElements[ReducedNbr] / x;
  return RealBandDiagonalSymmetricMatrix(Diagonal, UpperDiagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

RealBandDiagonalSymmetricMatrix& RealBandDiagonalSymmetricMatrix::operator += (const RealBandDiagonalSymmetricMatrix& M) 
{
  if (this->NbrRow != M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < ReducedNbr; i++)
    {
      this->DiagonalElements[i] += M.DiagonalElements[i];
      this->UpperOffDiagonalElements[i] += M.UpperOffDiagonalElements[i];      
    }
  this->DiagonalElements[ReducedNbr] += DiagonalElements[ReducedNbr];
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

RealBandDiagonalSymmetricMatrix& RealBandDiagonalSymmetricMatrix::operator -= (const RealBandDiagonalSymmetricMatrix& M) 
{
  if (this->NbrRow != M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  for (int i = 0; i < ReducedNbr; i++)
    {
      this->DiagonalElements[i] -= M.DiagonalElements[i];
      this->UpperOffDiagonalElements[i] -= M.UpperOffDiagonalElements[i];      
    }
  this->DiagonalElements[ReducedNbr] -= DiagonalElements[ReducedNbr];
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
  int ReducedNbr = this->NbrRow - 1;
  for (int i = 0; i < ReducedNbr; i++)
    {
      this->DiagonalElements[i] *= x;
      this->UpperOffDiagonalElements[i] *= x;      
    }
  this->DiagonalElements[ReducedNbr] *= x;
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
  int ReducedNbr = this->NbrRow - 1;
  for (int i = 0; i < ReducedNbr; i++)
    {
      this->DiagonalElements[i] /= x;
      this->UpperOffDiagonalElements[i] /= x;      
    }
  this->DiagonalElements[ReducedNbr] /= x;
  return *this;
}

// get a matrix element 
// 
// i = Row number
// j = Column number
// return value = matrix element M_(i,j)

double RealBandDiagonalSymmetricMatrix::GetElement(int i, int j)
{
  if (i == j)
    return this->DiagonalElements[i];
  if (j < i)
    {
      int tmp = i;
      i = j;
      j = tmp;
    }
  if ((j - i) == 1)
    return this->UpperOffDiagonalElements[i];
  return 0.0;
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
  if (this->NbrRow == 0)
    return 0.0;
  double d0 = this->DiagonalElements[0];
  if (this->NbrRow == 1)
    return d0;
  double d1 = d0 * this->DiagonalElements[0] - this->UpperOffDiagonalElements[0] * this->UpperOffDiagonalElements[0];
  if (this->NbrRow == 2)
    return d0;
  double d;
  for (int i = 2; i < this->NbrRow; i++)
    {
      d = this->DiagonalElements[i] * d1 - this->UpperOffDiagonalElements[i - 1] * this->UpperOffDiagonalElements[i - 1] * d0;
      d0 = d1;
      d1 = d;
    }
  return d1;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const RealBandDiagonalSymmetricMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      int j = 0;
      for (; j < (i-1); j++)
	Str << "0    ";
      if (i > 0)
	{
	  Str << P.UpperOffDiagonalElements[(i - 1)] << "    ";
	  j++;
	}
      Str << P.DiagonalElements[i] << "    ";
      j++;
      if (i < (P.NbrRow -1))
	{
	  Str << P.UpperOffDiagonalElements[i] << "    ";
	  j++;
	}
      for (; j < P.NbrColumn; j++)
	Str << "0    ";
      Str << endl;
    }
/*  for (int i = 0; i < (P.NbrColumn - 1); i++)
    if (i > 0)
      {
	Str << P.DiagonalElements[i] << "    " << P.UpperOffDiagonalElements[i] << endl;
      }
  Str << P.DiagonalElements[P.NbrColumn - 1] << endl;*/
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
  for (int i = 0; i < (P.NbrRow - 1); i++)
    {
      Str << "{";
      int j = 0;
      for (; j < (i-1); j++)
	Str << "0,";
      if (i > 0)
	{
	  Str << P.UpperOffDiagonalElements[(i - 1)] << ",";
	  j++;
	}
      Str << P.DiagonalElements[i] << ",";
      j++;
      if (i < (P.NbrRow -1))
	{
	  Str << P.UpperOffDiagonalElements[i];
	  if (j != (P.NbrColumn - 1))
	    Str << ",";
	  j++;
	}
      for (; j < (P.NbrColumn - 1); j++)
	Str << "0,";
      if (j == (P.NbrColumn - 1))
	Str << "0";
      Str << "},";
    }
  Str << "{";
  int j = 0;
  for (; j < (P.NbrRow-2); j++)
    Str << "0,";
  if (P.NbrRow > 0)
    {
      Str << P.UpperOffDiagonalElements[P.NbrRow - 2] << ",";
    }
  Str << P.DiagonalElements[P.NbrRow - 1];
  Str << "}}";
  return Str;
}

#endif
