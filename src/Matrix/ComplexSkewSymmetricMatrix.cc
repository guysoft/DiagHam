////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of complex skew symmetric matrix                  //
//                                                                            //
//                        last modification : 19/08/2004                      //
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


#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "GeneralTools/ListIterator.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include <stdlib.h>


using std::endl;


// default constructor
//

ComplexSkewSymmetricMatrix::ComplexSkewSymmetricMatrix() 
{
  this->Dummy = 0.0;
  this->RealOffDiagonalElements = 0;
  this->ImaginaryOffDiagonalElements =  0;
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = 0;
  this->MatrixType = Matrix::ComplexElements | Matrix::Antisymmetric;
}

// constructor for an empty matrix
//
// dimension = matrix dimension
// zero = true if matrix has to be filled with zeros

ComplexSkewSymmetricMatrix::ComplexSkewSymmetricMatrix(int dimension, bool zero) 
{
  this->Dummy = 0.0;
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Antisymmetric;
  if (this->NbrRow > 1)
    {
      this->RealOffDiagonalElements = new double [(this->NbrRow * (this->NbrRow - 1)) / 2];
      this->ImaginaryOffDiagonalElements = new double [(this->NbrRow * (this->NbrRow - 1)) / 2];
    }
  else
    {    
      this->RealOffDiagonalElements = new double [1];
      this->ImaginaryOffDiagonalElements = new double [1];
    }
  if (zero == true)
    {
      int pos = 0;
      for (int i = 0; i < this->NbrRow; i++)
	{
	  for (int j = i + 1; j < this->NbrRow; j++)
	    {
	      this->RealOffDiagonalElements[pos] = 0.0;
	      this->ImaginaryOffDiagonalElements[pos] = 0.0;
	      pos++;
	    }
	}
    }
}

// constructor from matrix elements (without duplicating datas)
//
// realUpperDiagonal = pointer to real part of the upper-diagonal elements
// imaginaryUpperDiagonal = pointer to imaginary part of the upper-diagonal elements
// dimension = matrix dimension

ComplexSkewSymmetricMatrix::ComplexSkewSymmetricMatrix(double* realUpperDiagonal, double* imaginaryUpperDiagonal, int dimension) 
{
  this->Dummy = 0.0;
  this->RealOffDiagonalElements = realUpperDiagonal;  
  this->ImaginaryOffDiagonalElements = imaginaryUpperDiagonal;  
  this->Flag.Initialize();
  this->NbrRow = dimension;
  this->NbrColumn = dimension;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Antisymmetric;
}

// copy constructor (without duplicating datas)
//
// M = matrix to copy

ComplexSkewSymmetricMatrix::ComplexSkewSymmetricMatrix(const ComplexSkewSymmetricMatrix& M) 
{
  this->Dummy = 0.0;
  this->RealOffDiagonalElements = M.RealOffDiagonalElements;  
  this->ImaginaryOffDiagonalElements = M.ImaginaryOffDiagonalElements;  
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = (this->TrueNbrRow - this->NbrRow);
  this->MatrixType = Matrix::ComplexElements | Matrix::Antisymmetric;
}

// destructor
//

ComplexSkewSymmetricMatrix::~ComplexSkewSymmetricMatrix() 
{
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->Flag.Used() == true) && (this->Flag.Shared() == false))
      {
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
}

// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::operator = (const ComplexSkewSymmetricMatrix& M) 
{
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
  this->RealOffDiagonalElements = M.RealOffDiagonalElements;
  this->ImaginaryOffDiagonalElements = M.ImaginaryOffDiagonalElements;
  this->Flag = M.Flag;
  this->NbrRow = M.NbrRow;
  this->NbrColumn = M.NbrColumn;
  this->TrueNbrRow = M.TrueNbrRow;
  this->TrueNbrColumn = M.TrueNbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  return *this;
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* ComplexSkewSymmetricMatrix::Clone ()
{
  return ((Matrix*) new ComplexSkewSymmetricMatrix (*this));
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void ComplexSkewSymmetricMatrix::SetMatrixElement(int i, int j, double x)
{
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void ComplexSkewSymmetricMatrix::SetMatrixElement(int i, int j, const Complex& x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i == j))
    return;
  if (i > j)
    {
      i -= (j * (j + 1)) / 2 - j * (this->NbrRow + this->Increment - 1) + 1;
      this->RealOffDiagonalElements[i] = -x.Re;
      this->ImaginaryOffDiagonalElements[i] = -x.Im;
    }
  else
    {
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      this->RealOffDiagonalElements[j] = x.Re;
      this->ImaginaryOffDiagonalElements[j] = x.Im;
    }
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void ComplexSkewSymmetricMatrix::AddToMatrixElement(int i, int j, double x)
{
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element
void ComplexSkewSymmetricMatrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  if ((i >= this->NbrRow) || (j >= this->NbrColumn) || (i == j))
    return;
  if (i > j)
    {
      i -= (j * (j + 1)) / 2 - j * (this->NbrRow + this->Increment - 1) + 1;
      this->RealOffDiagonalElements[i] -= x.Re;
      this->ImaginaryOffDiagonalElements[i] -= x.Im;
    }
  else
    {
      j -= (i * (i + 1)) / 2 - i * (this->NbrRow + this->Increment - 1) + 1;
      this->RealOffDiagonalElements[j] += x.Re;
      this->ImaginaryOffDiagonalElements[j] += x.Im;
    }
}

// get reference of a given matrix element supposing i < j
//
// i = line position
// j = column position
// return value = reference om matrix elememt

double& ComplexSkewSymmetricMatrix::operator () (int i, int j)
{
  return this->Dummy;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexSkewSymmetricMatrix::Resize (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      this->Increment = this->TrueNbrRow - this->NbrRow;
      return;
    }
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpRealOffDiag = new double [Tot];
  double* TmpImaginaryOffDiag = new double [Tot];
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  TmpRealOffDiag[k] = this->RealOffDiagonalElements[l];
	  TmpImaginaryOffDiag[k] = this->ImaginaryOffDiagonalElements[l];
	  ++l;
	  ++k;
	}
      l += this->Increment;
      for (int j = this->NbrRow; j < nbrRow; j++)
	{
	  TmpRealOffDiag[k] = 0.0;
	  TmpImaginaryOffDiag[k] = 0.0;
	  ++k;
	}      
    }
  for (int i = this->NbrRow * (this->NbrRow - 1); i < Tot; i++)
    {
      TmpRealOffDiag[i] = 0.0;
      TmpImaginaryOffDiag[i] = 0.0;
    }
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  this->RealOffDiagonalElements = TmpRealOffDiag;
  this->ImaginaryOffDiagonalElements = TmpImaginaryOffDiag;
  this->Flag.Initialize();
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void ComplexSkewSymmetricMatrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow != nbrColumn)
    return;
  if (nbrRow <= this->TrueNbrRow)
    {
      if (this->NbrRow < nbrRow)
	{
	  int Tot = (nbrRow * (nbrRow - 1));
	  int k = (this->NbrRow - 1);
	  for (int i = 0; i < (this->NbrRow - 1); i++)
	    {
	      for (int j = this->NbrRow; j < nbrRow; j++)
		{
		  this->RealOffDiagonalElements[k] = 0.0;
		  this->ImaginaryOffDiagonalElements[k] = 0.0;
		  ++k;
		}
	      k += (this->NbrRow - 2 - i);
	    }
	  for (int i = this->NbrRow * (this->NbrRow - 1); i < Tot; i++)
	    {
	      this->RealOffDiagonalElements[i] = 0.0;
	      this->ImaginaryOffDiagonalElements[i] = 0.0;
	    }
	}
      this->NbrRow = nbrRow;
      this->NbrColumn = nbrColumn;
      this->Increment = (this->TrueNbrRow - this->NbrRow);
      return;
    }
  int Tot = (nbrRow * (nbrRow - 1)) / 2;
  double* TmpRealOffDiag = new double [Tot];
  double* TmpImaginaryOffDiag = new double [Tot];
  int k = 0;
  int l = 0;
  for (int i = 0; i < (this->NbrRow - 1); i++)
    {
      for (int j = i + 1; j < this->NbrRow; j++)
	{
	  TmpRealOffDiag[k] = this->RealOffDiagonalElements[l];
	  TmpImaginaryOffDiag[k] = this->ImaginaryOffDiagonalElements[l];
	  ++l;
	  ++k;
	}
      l += this->Increment;
      for (int j = this->NbrRow; j < nbrRow; j++)
	{
	  TmpRealOffDiag[k] = 0.0;
	  TmpImaginaryOffDiag[k] = 0.0;
	  ++k;
	}      
    }
  for (int i = this->NbrRow * (this->NbrRow - 1); i < Tot; i++)
    {
      TmpRealOffDiag[i] = 0.0;
      TmpImaginaryOffDiag[i] = 0.0;
    }
  if ((this->RealOffDiagonalElements != 0) && (this->ImaginaryOffDiagonalElements != 0) 
      && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealOffDiagonalElements;
	delete[] this->ImaginaryOffDiagonalElements;
      }
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
  this->TrueNbrRow = this->NbrRow;
  this->TrueNbrColumn = this->NbrColumn;
  this->Increment = this->TrueNbrRow - this->NbrRow;
  this->RealOffDiagonalElements = TmpRealOffDiag;
  this->ImaginaryOffDiagonalElements = TmpImaginaryOffDiag;
  this->Flag.Initialize();
}

// add two matrices
//
// M1 = first matrix
// M2 = second matrix
// return value = sum of the two matrices

ComplexSkewSymmetricMatrix operator + (const ComplexSkewSymmetricMatrix& M1, const ComplexSkewSymmetricMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexSkewSymmetricMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* RealOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  double* ImaginaryOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  int k = 0;
  int l1 = 0;
  int l2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  RealOffDiagonal[k] = M1.RealOffDiagonalElements[l1] + M2.RealOffDiagonalElements[l2];      
	  ImaginaryOffDiagonal[k++] = M1.ImaginaryOffDiagonalElements[l1++] + M2.ImaginaryOffDiagonalElements[l2++];      
	}
      l1 += M2.Increment;
      l2 += M2.Increment;
    }
  return ComplexSkewSymmetricMatrix(RealOffDiagonal, ImaginaryOffDiagonal, M1.NbrRow);
}

// substract two matrices
//
// M1 = first matrix
// M2 = matrix to substract to M1
// return value = difference of the two matrices

ComplexSkewSymmetricMatrix operator - (const ComplexSkewSymmetricMatrix& M1, const ComplexSkewSymmetricMatrix& M2)
{
  if (M1.NbrRow != M2.NbrRow)
    return ComplexSkewSymmetricMatrix();
  int ReducedNbr = M1.NbrRow - 1;
  double* RealOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  double* ImaginaryOffDiagonal = new double [M1.NbrRow * ReducedNbr];
  int k = 0;
  int l1 = 0;
  int l2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = 0; j < i; j++)
	{
	  RealOffDiagonal[k] = M1.RealOffDiagonalElements[l1] - M2.RealOffDiagonalElements[l2];      
	  ImaginaryOffDiagonal[k++] = M1.ImaginaryOffDiagonalElements[l1++] - M2.ImaginaryOffDiagonalElements[l2++];      
	}
      l1 += M2.Increment;
      l2 += M2.Increment;
    }
  return ComplexSkewSymmetricMatrix(RealOffDiagonal, ImaginaryOffDiagonal, M1.NbrRow);
}

// multiply a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexSkewSymmetricMatrix operator * (const ComplexSkewSymmetricMatrix& M, double x) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* RealOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  double* ImaginaryOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  RealOffDiagonal[k] = M.RealOffDiagonalElements[k2] * x;
	  ImaginaryOffDiagonal[k++] = M.ImaginaryOffDiagonalElements[k2++] * x;
	}
      k2 += M.Increment;
    }
  return ComplexSkewSymmetricMatrix(RealOffDiagonal, ImaginaryOffDiagonal, M.NbrRow);
}

// multiply a matrix by a real number (left multiplication)
//
// M = source matrix
// x = real number to use
// return value = product result

ComplexSkewSymmetricMatrix operator * (double x, const ComplexSkewSymmetricMatrix& M) 
{
  int ReducedNbr = M.NbrRow - 1;
  double* RealOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  double* ImaginaryOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  RealOffDiagonal[k] = M.RealOffDiagonalElements[k2] * x;
	  ImaginaryOffDiagonal[k++] = M.ImaginaryOffDiagonalElements[k2++] * x;
	}
      k2 += M.Increment;
    }
  return ComplexSkewSymmetricMatrix(RealOffDiagonal, ImaginaryOffDiagonal, M.NbrRow);
}

// divide a matrix by a real number (right multiplication)
//
// M = source matrix
// x = real number to use
// return value = division result

ComplexSkewSymmetricMatrix operator / (const ComplexSkewSymmetricMatrix& M, double x) 
{
  x = 1.0 / x;
  int ReducedNbr = M.NbrRow - 1;
  double* RealOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  double* ImaginaryOffDiagonal = new double [(M.NbrRow * ReducedNbr) / 2];
  int k = 0;
  int k2 = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  RealOffDiagonal[k] = M.RealOffDiagonalElements[k2] * x;
	  ImaginaryOffDiagonal[k++] = M.ImaginaryOffDiagonalElements[k2++] * x;
	}
      k2 += M.Increment;
    }
  return ComplexSkewSymmetricMatrix(RealOffDiagonal, ImaginaryOffDiagonal, M.NbrRow);
}

// add two matrices
//
// M = matrix to add to current matrix
// return value = reference on current matrix

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::operator += (const ComplexSkewSymmetricMatrix& M) 
{
  if (this->NbrRow < M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  int k = 0;
  int k2 = 0;  
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] += M.RealOffDiagonalElements[k2];
	  this->ImaginaryOffDiagonalElements[k++] += M.ImaginaryOffDiagonalElements[k2++];
	}
      k += this->Increment;
      k2 += M.Increment;
    }
  return *this;
}

// substract two matrices
//
// M = matrix to substract to current matrix
// return value = reference on current matrix

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::operator -= (const ComplexSkewSymmetricMatrix& M) 
{
  if (this->NbrRow < M.NbrRow)
    return *this;
  int ReducedNbr = M.NbrRow - 1;
  int k = 0;
  int k2 = 0;  
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] -= M.RealOffDiagonalElements[k2];
	  this->ImaginaryOffDiagonalElements[k++] -= M.ImaginaryOffDiagonalElements[k2++];
	}
      k += this->Increment;
      k2 += M.Increment;
    }
  return *this;
}

// multiply a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::operator *= (double x) 
{
  if (this->NbrRow == 0)
    return *this;
  int ReducedNbr = this->NbrRow - 1;
  int k = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] *= x;
	  this->ImaginaryOffDiagonalElements[k++] *= x;
	}
      k += this->Increment;
    }
  return *this;
}

// divide a matrix by a real number
//
// x = real number to use
// return value = reference on current matrix

ComplexSkewSymmetricMatrix& ComplexSkewSymmetricMatrix::operator /= (double x)
{
  if (this->NbrRow == 0)
    return *this;
  x = 1.0 / x;
  int ReducedNbr = this->NbrRow - 1;
  int k = 0;
  for (int i = 0; i < ReducedNbr; i++)
    {
      for (int j = i; j < ReducedNbr; j++)
	{
	  this->RealOffDiagonalElements[k] *= x;
	  this->ImaginaryOffDiagonalElements[k++] *= x;
	}
      k += this->Increment;
    }
  return *this;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ComplexSkewSymmetricMatrix::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  Complex x;
  if ((V1.Dimension != this->NbrRow) || (V2.Dimension != this->NbrColumn))
    return x;
  for (int i = 0; i < this->NbrRow ; i++)
    {
      Complex x2;
      int l = (i - 1);
      for (int k = 0; k < i; k++)
	{
	  x2.Re -= (this->RealOffDiagonalElements[l] * V2.RealComponents[k] - 
		    this->ImaginaryOffDiagonalElements[l] * V2.ImaginaryComponents[k]);
	  x2.Im -= (this->ImaginaryOffDiagonalElements[l] * V2.RealComponents[k] + 
		    this->RealOffDiagonalElements[l] * V2.ImaginaryComponents[k]);
	  l += (this->NbrColumn - 2 - k) + this->Increment;
	}
      l++;
      for (int k = i + 1; k < this->NbrColumn; k++)
	{
	  x2.Re += (this->RealOffDiagonalElements[l] * V2.RealComponents[k] -
		    this->ImaginaryOffDiagonalElements[l] * V2.ImaginaryComponents[k]);
	  x2.Im += (this->ImaginaryOffDiagonalElements[l] * V2.RealComponents[k] + 
		    this->RealOffDiagonalElements[l] * V2.ImaginaryComponents[k]);
	  ++l;
	}      
      x.Re += (x2.Re * V1.RealComponents[i] + x2.Im * V1.ImaginaryComponents[i]);
      x.Re += (x2.Im * V1.RealComponents[i] - x2.Re * V1.ImaginaryComponents[i]);
    }
  return x;
}

// conjugate an hermitian matrix with an unitary matrix (Uh M U)
//
// UnitaryM = unitary matrix to use
// return value = conjugated matrix

Matrix* ComplexSkewSymmetricMatrix::Conjugate(ComplexMatrix& UnitaryM)
{
  if (UnitaryM.NbrRow != this->NbrColumn)
    return 0;
  int NbrOffDiag = (UnitaryM.NbrColumn * (UnitaryM.NbrColumn - 1)) / 2;
  double* TmpRealOffDiag = new double [NbrOffDiag];
  double* TmpImaginaryOffDiag = new double [NbrOffDiag];
  int i2 = 0;
  int ReducedNbrColumn = UnitaryM.NbrColumn - 1;
  int Inc = this->NbrColumn - 3 + this->Increment;
  Complex tmp1;
  int k;
  int l;
  for (int i = 0; i < ReducedNbrColumn; i++)
    {
      for (int m = i + 1; m < UnitaryM.NbrColumn; m++)
	{    
	  TmpRealOffDiag[i2] = 0.0;
	  TmpImaginaryOffDiag[i2] = 0.0;
	  for (int j = 0; j < this->NbrColumn; j++)
	    {
	      tmp1 = 0.0;
	      k = 0;
	      l = (j - 1);
	      for (; k < j; k++)
		{
		  tmp1.Re -= (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k] - 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k]);
		  tmp1.Im -= (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k] + 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k]);
		  l += Inc - k;
		}
	      l++;
	      k++;
	      for (; k < this->NbrColumn; k++)
		{
		  tmp1.Re += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k] - 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k]);
		  tmp1.Im += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k] + 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k]);
		  ++l;
		}
	      TmpRealOffDiag[i2] += tmp1.Re * UnitaryM.Columns[i].RealComponents[j] + tmp1.Im * UnitaryM.Columns[i].ImaginaryComponents[j];
	      TmpImaginaryOffDiag[i2] += tmp1.Im * UnitaryM.Columns[i].RealComponents[j] - tmp1.Re * UnitaryM.Columns[i].ImaginaryComponents[j] ;
	    }
	  ++i2;
	}    
    }
  return new ComplexSkewSymmetricMatrix(TmpRealOffDiag, TmpImaginaryOffDiag, UnitaryM.NbrColumn);
}

// conjugate a block of the matrix with an unitary matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// sourcePosition = index of the row where the block to conjugate starts
// destinationPosition = index of the row where the conjugated block has to be stored
// matrix = matrix where result has to be stored

void ComplexSkewSymmetricMatrix::Conjugate(ComplexMatrix& UnitaryM, int sourcePosition, int destinationPosition,
					   ComplexSkewSymmetricMatrix& matrix)
{
  if (((UnitaryM.NbrRow + sourcePosition) > this->NbrColumn) || 
      ((UnitaryM.NbrColumn + destinationPosition) > matrix.NbrColumn))
    return;
  int i2 = (destinationPosition - (destinationPosition * (destinationPosition + 1)) / 2 +
	    destinationPosition * (matrix.NbrRow + matrix.Increment - 1));
  int ReducedNbrColumn = UnitaryM.NbrColumn - 1;
  for (int i = 0; i < ReducedNbrColumn; i++)
    {
      for (int m = i + 1; m < UnitaryM.NbrColumn; m++)
	{    
	  matrix.RealOffDiagonalElements[i2] = 0.0;
	  matrix.ImaginaryOffDiagonalElements[i2] = 0.0;
	  for (int j = 0; j < UnitaryM.NbrRow; j++)
	    {
	      Complex tmp1;
	      int k = 0;
	      int l = ((j + sourcePosition) - 1 - (sourcePosition * (sourcePosition + 1)) / 2 +
		       sourcePosition * (this->NbrRow + this->Increment - 1));
	      for (; k < j; k++)
		{
		  tmp1.Re -= (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k] + 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k]);
		  tmp1.Im -= (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k] - 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k]);
		  l += (this->NbrColumn - 2 - k - sourcePosition) + this->Increment;
		}
	      l++;
	      k++;
	      for (; k < UnitaryM.NbrRow; k++)
		{
		  tmp1.Re += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k] - 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k]);
		  tmp1.Im += (this->RealOffDiagonalElements[l] * UnitaryM.Columns[m].ImaginaryComponents[k] + 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryM.Columns[m].RealComponents[k]);
		  ++l;
		}
	      matrix.RealOffDiagonalElements[i2] += tmp1.Re * UnitaryM.Columns[i].RealComponents[j] + tmp1.Im * UnitaryM.Columns[i].ImaginaryComponents[j];
	      matrix.RealOffDiagonalElements[i2] += tmp1.Im * UnitaryM.Columns[i].RealComponents[j] - tmp1.Re * UnitaryM.Columns[i].ImaginaryComponents[j];
	    }
	  i2++;
	}
      i2 += matrix.NbrColumn - destinationPosition - UnitaryM.NbrColumn + matrix.Increment;
    }
}

// conjugate a block of the matrix (in the upper diagonal part) with two matrix matrix (Vt M U)
//
// UnitaryMl = unitary matrix to use at the left hand side
// UnitaryMr = unitary matrix to use at the right hand side
// sourceRowIndex = index of the row where the block to conjugate starts
// sourceColumnIndex = index of the column where the block to conjugate starts
// destinationRowIndex = index of the row where the conjugated block has to be stored
// destinationColumnIndex = index of the column where the conjugated block has to be stored
// matrix = matrix where result has to be stored

void ComplexSkewSymmetricMatrix::Conjugate(ComplexMatrix& UnitaryMl, ComplexMatrix& UnitaryMr, int sourceRowIndex, 
					   int sourceColumnIndex, int destinationRowIndex,
					   int destinationColumnIndex, ComplexSkewSymmetricMatrix& matrix)
{
  if (((UnitaryMr.NbrRow + sourceColumnIndex) > this->NbrColumn) || 
      ((UnitaryMl.NbrRow + sourceRowIndex) > this->NbrRow) || 
      ((UnitaryMr.NbrColumn + destinationColumnIndex) > matrix.NbrColumn) || 
      ((UnitaryMl.NbrColumn + destinationRowIndex) > matrix.NbrRow))
    return;
  for (int i = 0; i < UnitaryMl.NbrColumn; i++)
    {
      int i2 = (destinationColumnIndex - 1 - ((i + destinationRowIndex) * ((i + destinationRowIndex) + 1)) / 2 +
		(i + destinationRowIndex) * (matrix.NbrRow + matrix.Increment - 1));
      for (int m = 0; m < UnitaryMr.NbrColumn; m++)
	{    
	  matrix.RealOffDiagonalElements[i2] = 0.0;
	  matrix.ImaginaryOffDiagonalElements[i2] = 0.0;
	  for (int j = 0; j < UnitaryMl.NbrRow; j++)
	    {
	      Complex tmp1;
	      int l = (sourceColumnIndex - 1 - 
		       ((j + sourceRowIndex) * ((sourceRowIndex + j) + 1)) / 2 +
		       (sourceRowIndex + j) * (this->NbrRow + this->Increment - 1));
	      for (int k = 0; k < UnitaryMr.NbrRow; k++)
		{
		  tmp1.Re += (this->RealOffDiagonalElements[l] * UnitaryMr.Columns[m].RealComponents[k] - 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryMr.Columns[m].ImaginaryComponents[k]);
		  tmp1.Im += (this->RealOffDiagonalElements[l] * UnitaryMr.Columns[m].ImaginaryComponents[k] + 
			      this->ImaginaryOffDiagonalElements[l] * UnitaryMr.Columns[m].RealComponents[k]);
		  ++l;
		}
	      matrix.RealOffDiagonalElements[i2] += tmp1.Re * UnitaryMl.Columns[i].RealComponents[j] + tmp1.Im * UnitaryMl.Columns[i].ImaginaryComponents[j];
	      matrix.ImaginaryOffDiagonalElements[i2] += tmp1.Im * UnitaryMl.Columns[i].RealComponents[j] - tmp1.Re * UnitaryMl.Columns[i].ImaginaryComponents[j] ;
	    }
	  i2++;
	}
    }
  return;
}

// evaluate matrix trace
//
// return value = matrix trace 

double ComplexSkewSymmetricMatrix::Tr () 
{
  return 0.0;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double ComplexSkewSymmetricMatrix::Det () 
{
  return 1.0;
}

// evaluate matrix pfaffian
//
// return value = matrix pfaffian 

Complex ComplexSkewSymmetricMatrix::Pfaffian()
{
  Complex Tmp;
  if (this->NbrColumn == 2)
    {
      Tmp.Re = this->RealOffDiagonalElements[0];
      Tmp.Im = this->ImaginaryOffDiagonalElements[0];
      return Tmp;
    }
  if (this->NbrColumn == 4)
    {
      Tmp.Re = (this->RealOffDiagonalElements[0] * this->RealOffDiagonalElements[5 + (2 * this->Increment)]
		 - this->ImaginaryOffDiagonalElements[0] * this->ImaginaryOffDiagonalElements[5 + (2 * this->Increment)]);
      Tmp.Im = (this->ImaginaryOffDiagonalElements[0] * this->RealOffDiagonalElements[5 + (2 * this->Increment)]
		+ this->RealOffDiagonalElements[0] * this->ImaginaryOffDiagonalElements[5 + (2 * this->Increment)]);
      Tmp.Re += (this->RealOffDiagonalElements[1] * this->RealOffDiagonalElements[4 + this->Increment]
		 - this->ImaginaryOffDiagonalElements[1] * this->ImaginaryOffDiagonalElements[4 + this->Increment]);
      Tmp.Im += (this->ImaginaryOffDiagonalElements[1] * this->RealOffDiagonalElements[4 + this->Increment]
		 + this->RealOffDiagonalElements[1] * this->ImaginaryOffDiagonalElements[4 + this->Increment]);
      Tmp.Re += (this->RealOffDiagonalElements[2] * this->RealOffDiagonalElements[3 + this->Increment]
		 - this->ImaginaryOffDiagonalElements[2] * this->ImaginaryOffDiagonalElements[3 + this->Increment]);
      Tmp.Im += (this->ImaginaryOffDiagonalElements[2] * this->RealOffDiagonalElements[3 + this->Increment]
		 + this->RealOffDiagonalElements[2] * this->ImaginaryOffDiagonalElements[3 + this->Increment]);
      return Tmp;
    }
  return Tmp;
}

// Output Stream overload
//
// Str = reference on output stream
// P = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& Str, const ComplexSkewSymmetricMatrix& P)
{
  for (int i = 0; i < P.NbrRow; i++)
    {
      int pos = (i - 1);
      for (int j = 0; j < i; j ++)
	{
	  Str << -(P.RealOffDiagonalElements[pos]) << "    ";
	  if (P.ImaginaryOffDiagonalElements[pos] > 0.0)
	    Str << -P.ImaginaryOffDiagonalElements[pos] << "i    ";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << -P.ImaginaryOffDiagonalElements[pos] << "i    ";
	    else
	      Str << "    ";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      Str << "0    ";
      pos++;
      for (int j = i + 1; j < P.NbrRow; j++)
	{
	  Str << P.RealOffDiagonalElements[pos] << "    ";
	  if (P.ImaginaryOffDiagonalElements[pos] < 0.0)
	    Str << P.ImaginaryOffDiagonalElements[pos] << "i    ";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << P.ImaginaryOffDiagonalElements[pos] << "i    ";
	    else
	      Str << "    ";
	  ++pos;
	}
      Str << endl;
    }
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// P = matrix to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexSkewSymmetricMatrix& P)
{
  Str << "{";
  for (int i = 0; i < P.NbrRow; i++)
    {
      Str << "{";
      int pos = (i - 1);
      for (int j = 0; j < i; j ++)
	{
	  if ((P.RealOffDiagonalElements[pos] != 0) || (P.ImaginaryOffDiagonalElements[pos] == 0))
	    Str << -P.RealOffDiagonalElements[pos];
	  if (P.ImaginaryOffDiagonalElements[pos] > 0.0)
	    Str << -P.ImaginaryOffDiagonalElements[pos] << "I";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << -P.ImaginaryOffDiagonalElements[pos] << "I";
	  Str << ",";
	  pos += (P.NbrRow - j - 2) + P.Increment;
	}
      Str << "0";
      if (i != (P.NbrRow - 1))
	{
	  Str << ",";	  
	  pos++;
	  for (int j = i + 1; j < (P.NbrRow - 1); j++)
	    {
	      if ((P.RealOffDiagonalElements[pos] != 0) || (P.ImaginaryOffDiagonalElements[pos] == 0))
		Str << P.RealOffDiagonalElements[pos];
	      if (P.ImaginaryOffDiagonalElements[pos] < 0.0)
		Str << P.ImaginaryOffDiagonalElements[pos] << "I";
	      else
		if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
		  Str << "+" << P.ImaginaryOffDiagonalElements[pos] << "I";
	      Str << ",";
	      ++pos;
	    }
	  Str << P.RealOffDiagonalElements[pos];
	  if (P.ImaginaryOffDiagonalElements[pos] < 0.0)
	    Str << P.ImaginaryOffDiagonalElements[pos] << "I";
	  else
	    if (P.ImaginaryOffDiagonalElements[pos] != 0.0)
	      Str << "+" << P.ImaginaryOffDiagonalElements[pos] << "I";
	  Str << "},";
	}
      else
	Str << "}";
    }
  Str << "}";
  return Str;
}

