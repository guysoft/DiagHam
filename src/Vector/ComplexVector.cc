////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                   class for n dimensional complex vector                   //
//                                                                            //
//                        last modification : 17/01/2001                      //
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


#include "Vector/ComplexVector.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"


#include <fstream>

using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;

  
// default constructor
//

ComplexVector::ComplexVector() 
{
  this->RealComponents = 0;
  this->ImaginaryComponents = 0;
  this->Flag.Initialize();
  this->Dimension = 0;
  this->TrueDimension = this->Dimension;
  this->VectorType = Vector::ComplexDatas;
}

// constructor for an empty real vector
//
// size = Vector Dimension 
// zeroFlag = true if all coordinates have to be set to zero

ComplexVector::ComplexVector(int size, bool zeroFlag) 
{
  this->Dimension = size;
  this->TrueDimension = this->Dimension;
  this->Flag.Initialize();
  this->RealComponents = new double [this->Dimension];
  this->ImaginaryComponents = new double [this->Dimension];
  this->VectorType = Vector::ComplexDatas;
  if (zeroFlag == true)
    for (int i = 0; i < this->Dimension; ++i)
      {
	this->RealComponents[i] = 0;
	this->ImaginaryComponents[i] = 0;
      }
}

// constructor from arrays of doubles
//
// real = array of doubles corresponding to real part
// imaginary = array of doubles corresponding to imaginary part
// size = Vector Dimension 

ComplexVector::ComplexVector(double* real, double* imaginary, int size) 
{
  this->Flag.Initialize();
  this->Dimension = size;
  this->VectorType = Vector::ComplexDatas;
  this->TrueDimension = this->Dimension;
  this->RealComponents = real;
  this->ImaginaryComponents = imaginary;
}

// copy constructor
//
// vector = vector to copy
// duplicateFlag = true if datas have to be duplicated

ComplexVector::ComplexVector(const ComplexVector& vector, bool duplicateFlag) 
{
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->VectorType = Vector::ComplexDatas;
  if (vector.Dimension == 0)
    {
      this->Flag.Initialize();
      this->RealComponents = 0;
      this->ImaginaryComponents = 0;
    }
  else
    if (duplicateFlag == false)
      {
	this->Flag = vector.Flag;
	this->RealComponents = vector.RealComponents;
	this->ImaginaryComponents = vector.ImaginaryComponents;
      }
    else
      {
	this->Flag.Initialize();
	this->RealComponents = new double [this->TrueDimension];
	this->ImaginaryComponents = new double [this->TrueDimension];
	for (int i = 0; i < this->Dimension; i++)
	  {
	    this->RealComponents[i] = vector.RealComponents[i];
	    this->ImaginaryComponents[i] = vector.ImaginaryComponents[i];
	  }
      }
}

// copy constructor from a real vector
//
// vector = vector to copy
// duplicateFlag = true if datas have to be duplicated

ComplexVector::ComplexVector(const RealVector& vector, bool duplicateFlag) 
{
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->VectorType = Vector::ComplexDatas;
  if (vector.Dimension == 0)
    {
      this->Flag.Initialize();
      this->RealComponents = 0;
      this->ImaginaryComponents = 0;
    }
  else
    {
      this->Flag.Initialize();
      this->RealComponents = new double [this->Dimension];
      this->ImaginaryComponents = new double [this->Dimension];
      for (int i = 0; i < this->Dimension; ++i)
	{
	  this->RealComponents[i] = vector.Components[i];
	  this->ImaginaryComponents[i] = 0.0;
	}
    }
}

// destructor
//

ComplexVector::~ComplexVector () 
{
  if ((this->RealComponents != 0) && (this->ImaginaryComponents != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealComponents;
	delete[] this->ImaginaryComponents;
      }
}

// assignement
//
// vector = vector to assign
// return value = reference on current vector

ComplexVector& ComplexVector::operator = (const ComplexVector& vector) 
{
  if ((this->RealComponents != 0) && (this->ImaginaryComponents != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealComponents;
	delete[] this->ImaginaryComponents;
      }
  this->Flag = vector.Flag;
  this->RealComponents = vector.RealComponents;
  this->ImaginaryComponents = vector.ImaginaryComponents;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  return *this;
}

// Resize vector
//
// dimension = new dimension

void ComplexVector::Resize (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      return;
    }
  double* TmpVector1 = new double [dimension];
  double* TmpVector2 = new double [dimension];
  int i = 0;
  for (; i < this->Dimension; ++i)
    {
      TmpVector1[i] = this->RealComponents[i];
      TmpVector2[i] = this->ImaginaryComponents[i];
    }
  if ((this->RealComponents != 0) && (this->ImaginaryComponents != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealComponents;
	delete[] this->ImaginaryComponents;
      }
  this->Flag.Initialize();
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->RealComponents = TmpVector1;
  this->ImaginaryComponents = TmpVector2;
}

// Resize vector and set to zero all components that have been added
//
// dimension = new dimension

void ComplexVector::ResizeAndClean (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      for (int i = this->Dimension; i < dimension; ++i)
	{
	  this->RealComponents[i] = 0.0;  
	  this->ImaginaryComponents[i] = 0.0;  
	}
      return;
    }
  double* TmpVector1 = new double [dimension];
  double* TmpVector2 = new double [dimension];
  int i = 0;
  for (; i < this->Dimension; ++i)
    {
      TmpVector1[i] = this->RealComponents[i];
      TmpVector2[i] = this->ImaginaryComponents[i];
    }
  for (; i < dimension; i++)
    {
      TmpVector1[i] = 0.0;  
      TmpVector2[i] = 0.0;  
    }
  if ((this->RealComponents != 0) && (this->ImaginaryComponents != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealComponents;
	delete[] this->ImaginaryComponents;
      }
  this->Flag.Initialize();
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->RealComponents = TmpVector1;
  this->ImaginaryComponents = TmpVector2;
}

// assignement from a real vector
//
// vector = vector to assign
// return value = reference on current vector

ComplexVector& ComplexVector::operator = (const RealVector& vector) 
{
  if ((this->RealComponents != 0) && (this->ImaginaryComponents != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealComponents;
	delete[] this->ImaginaryComponents;
      }
  this->Dimension = vector.Dimension;
  this->TrueDimension = this->Dimension;
  if (vector.Dimension == 0)
    {
      this->RealComponents = 0;
      this->ImaginaryComponents = 0;
    }
  this->RealComponents = new double [this->Dimension];
  this->ImaginaryComponents = new double [this->Dimension];
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] = vector.Components[i];
      this->ImaginaryComponents[i] = 0.0; 
    }
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

ComplexVector& ComplexVector::Copy (ComplexVector& vector, double coefficient)
{
  if (this->Dimension != vector.Dimension)
    this->Resize(vector.Dimension);
  for (int i = 0; i < this->Dimension; i++)
    {
      this->RealComponents[i] = vector.RealComponents[i] * coefficient;
      this->ImaginaryComponents[i] = vector.ImaginaryComponents[i] * coefficient;
    }
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

ComplexVector& ComplexVector::Copy (ComplexVector& vector, const Complex& coefficient)
{
  if (this->Dimension != vector.Dimension)
    this->Resize(vector.Dimension);
  for (int i = 0; i < this->Dimension; i++)
    {
      this->RealComponents[i] = vector.RealComponents[i] * coefficient.Re - vector.ImaginaryComponents[i] * coefficient.Im;
      this->ImaginaryComponents[i] = vector.RealComponents[i] * coefficient.Im + vector.ImaginaryComponents[i] * coefficient.Re;
    }
  return *this;
}

// create a new vector with same size and same type but non-initialized components
//
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* ComplexVector::EmptyClone(bool zeroFlag)
{
  return new ComplexVector(this->Dimension, zeroFlag);
}

// put all vector components to zero
//
// return value = reference on current vector

Vector& ComplexVector::ClearVector ()
{
  for (int i = 0; i < this->Dimension; i++)
    {
      this->RealComponents[i] = 0.0;  
      this->ImaginaryComponents[i] =0.0;
    }
  return *this;
}

// change sign of a vector
//
// return value = reference on current vector

ComplexVector& ComplexVector::operator - () 
{
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] *= -1;        
      this->ImaginaryComponents[i] *= -1;  
    }
  return *this;
}

// return a new vector with opposite sign form a given source vector
//
// V1 = source vector
// return value = new vector

ComplexVector operator - (const ComplexVector& V1) 
{
  if (V1.Dimension > 0)
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = -V1.RealComponents[i];
	  TmpComponents2[i] = -V1.ImaginaryComponents[i];
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// scalar product between two vectors
//
// V1 = first vector
// V2 = second vector
// return value = result of scalar product

Complex operator * (const ComplexVector& V1, const ComplexVector& V2) 
{
  int min = V1.Dimension;
  if (min > V2.Dimension)
    min = V2.Dimension;
  if (min == 0)
    return Complex();
  Complex tmp (V1.RealComponents[0] * V2.RealComponents[0] + V1.ImaginaryComponents[0] * V2.ImaginaryComponents[0], 
	       V1.RealComponents[0] * V2.ImaginaryComponents[0] - V1.ImaginaryComponents[0] * V2.RealComponents[0]);
  for (int i = 1; i < min; ++i)
    {
      tmp.Re += V1.RealComponents[i] * V2.RealComponents[i] + V1.ImaginaryComponents[i] * V2.ImaginaryComponents[i];
      tmp.Im += V1.RealComponents[i] * V2.ImaginaryComponents[i] - V1.ImaginaryComponents[i] * V2.RealComponents[i];
    }
  return tmp;
}

// scalar product between two vectors
//
// V1 = first vector
// V2 = second vector (real vector)
// return value = result of scalar product

Complex operator * (const ComplexVector& V1, const RealVector& V2) 
{
  int min = V1.Dimension;
  if (min > V2.Dimension)
    min = V2.Dimension;
  Complex tmp (V1.RealComponents[0] * V2.Components[0], - V1.ImaginaryComponents[0] * V2.Components[0]);
  for (int i = 1; i < min; ++i)
    {
      tmp.Re += V1.RealComponents[i] * V2.Components[i];
      tmp.Im -= V1.ImaginaryComponents[i] * V2.Components[i];
    }
  return tmp;
}

// scalar product between two vectors
//
// V1 = first vector (real vector)
// V2 = second vector
// return value = result of scalar product

Complex operator * (const RealVector& V1, const ComplexVector& V2) 
{
  int min = V1.Dimension;
  if (min > V2.Dimension)
    min = V2.Dimension;
  Complex tmp (V1.Components[0] * V2.RealComponents[0], V2.ImaginaryComponents[0] * V1.Components[0]);
  for (int i = 1; i < min; i++)
    {
      tmp.Re += V2.RealComponents[i] * V1.Components[i];
      tmp.Im += V2.ImaginaryComponents[i] * V1.Components[i];
    }
  return tmp;
}

// sum two vectors
//
// V1 = vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::operator += (const ComplexVector& V1) 
{
  if ((this->Dimension == 0) || (this->Dimension != V1.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] += V1.RealComponents[i];
      this->ImaginaryComponents[i] += V1.ImaginaryComponents[i];
    }
  return *this;
}

// sum two vectors
//
// V1 = real vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::operator += (const RealVector& V1) 
{
  if ((this->Dimension == 0) || (this->Dimension != V1.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    this->RealComponents[i] += V1.Components[i];
  return *this;
}

// sum two vectors
//
// vector = vector to add
// return value = reference on current vector

Vector& ComplexVector::operator += (const Vector& vector)
{
  if (vector.VectorType == Vector::RealDatas)
    return (*this += ((RealVector&) vector));
  if (vector.VectorType == Vector::ComplexDatas)
    return (*this += ((ComplexVector&) vector));
  return *this;
}

// add a linear combination to a given vector
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (double x, const ComplexVector& V) 
{
  if ((this->Dimension == 0) || (this->Dimension != V.Dimension))
    {
      return *this;
    }
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] += x * V.RealComponents[i];
      this->ImaginaryComponents[i] += x * V.ImaginaryComponents[i];
    }
  return *this;
}

// add a linear combination to a given vector, for a given range of indices
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (double x, const ComplexVector& V, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > V.Dimension))
    return *this;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->RealComponents[i] += x * V.RealComponents[i];
      this->ImaginaryComponents[i] += x * V.ImaginaryComponents[i];
    }
  return *this;
}

// add a linear combination to a given vector
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (const Complex& x, const ComplexVector& V) 
{
  if ((this->Dimension == 0) || (this->Dimension != V.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] += x.Re * V.RealComponents[i] - x.Im * V.ImaginaryComponents[i];
      this->ImaginaryComponents[i] += x.Re * V.ImaginaryComponents[i] + x.Im * V.RealComponents[i];
    }
  return *this;
}

// add a linear combination to a given vector, for a given range of indices
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (const Complex& x, const ComplexVector& V, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > V.Dimension))
    return *this;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->RealComponents[i] += x.Re * V.RealComponents[i] - x.Im * V.ImaginaryComponents[i];
      this->ImaginaryComponents[i] += x.Re * V.ImaginaryComponents[i] + x.Im * V.RealComponents[i];
    }
  return *this;
}

// add a linear combination of two vectors to a given vector
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (double x1, const ComplexVector& v1, double x2, 
						    const ComplexVector& v2)
{
  if ((this->Dimension == 0) || (this->Dimension != v1.Dimension) || (this->Dimension != v2.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] += x1 * v1.RealComponents[i] + x2 * v2.RealComponents[i];
      this->ImaginaryComponents[i] += x1 * v1.ImaginaryComponents[i] + x2 * v2.ImaginaryComponents[i];
    }
  return *this;
}

// add a linear combination of two vectors to a given vector, for a given range of indices
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (double x1, const ComplexVector& v1, double x2, 
						    const ComplexVector& v2, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > v2.Dimension) || 
      (LastComponent > v1.Dimension))
    return *this;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->RealComponents[i] += x1 * v1.RealComponents[i] + x2 * v2.RealComponents[i];
      this->ImaginaryComponents[i] += x1 * v1.ImaginaryComponents[i] + x2 * v2.ImaginaryComponents[i];
    }
  return *this;
}

// add a linear combination of two vectors to a given vector
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (const Complex& x1, const ComplexVector& v1, const Complex& x2, 
						    const ComplexVector& v2)
{
  if ((this->Dimension == 0) || (this->Dimension != v1.Dimension) || (this->Dimension != v2.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] += x1.Re * v1.RealComponents[i]  - x1.Im * v1.ImaginaryComponents[i] 
	+ x2.Re * v2.RealComponents[i] - x2.Im * v2.ImaginaryComponents[i];
      this->ImaginaryComponents[i] += x1.Im * v1.RealComponents[i]  + x1.Re * v1.ImaginaryComponents[i] 
	+ x2.Im * v2.RealComponents[i] + x2.Re * v2.ImaginaryComponents[i];
    }
  return *this;
}

// add a linear combination of two vectors to a given vector, for a given range of indices
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (const Complex& x1, const ComplexVector& v1, const Complex& x2, 
						    const ComplexVector& v2, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > v2.Dimension) || 
      (LastComponent > v1.Dimension))
    return *this;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->RealComponents[i] += x1.Re * v1.RealComponents[i]  - x1.Im * v1.ImaginaryComponents[i] 
	+ x2.Re * v2.RealComponents[i] - x2.Im * v2.ImaginaryComponents[i];
      this->ImaginaryComponents[i] += x1.Im * v1.RealComponents[i]  + x1.Re * v1.ImaginaryComponents[i] 
	+ x2.Im * v2.RealComponents[i] + x2.Re * v2.ImaginaryComponents[i];
    }
  return *this;
}

// substract two vectors
//
// V1 = first vector
// return value = reference on current vector

ComplexVector& ComplexVector::operator -= (const ComplexVector& V1) 
{
  if ((this->Dimension == 0) || (this->Dimension != V1.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] -= V1.RealComponents[i];
      this->ImaginaryComponents[i] -= V1.ImaginaryComponents[i];
    }
  return *this;
}

// substract two vectors
//
// V1 = first real vector
// return value = reference on current vector

ComplexVector& ComplexVector::operator -= (const RealVector& V1) 
{
  if ((this->Dimension == 0) || (this->Dimension != V1.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; i++)
    {
      this->RealComponents[i] -= V1.Components[i];
    }
  return *this;
}

// sum two vectors
//
// V1 = first vector
// V2 = second vector
// return value = resulting vector

ComplexVector operator + (const ComplexVector& V1, const ComplexVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = V1.RealComponents[i] + V2.RealComponents[i];
	  TmpComponents2[i] = V1.ImaginaryComponents[i] + V2.ImaginaryComponents[i];
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// sum two vectors with left one real
//
// V1 = first vector (real)
// V2 = second vector
// return value = resulting vector

ComplexVector operator + (const RealVector& V1, const ComplexVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = V1.Components[i] + V2.RealComponents[i];
	  TmpComponents2[i] = V2.ImaginaryComponents[i];
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// sum two vectors with right one real
//
// V1 = first vector
// V2 = second vector (real)
// return value = resulting vector

ComplexVector operator + (const ComplexVector& V1, const RealVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = V2.Components[i] + V1.RealComponents[i];
	  TmpComponents2[i] = V1.ImaginaryComponents[i];
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// substract two vectors
//
// V1 = first vector
// V2 = second vector
// return value = resulting vector

ComplexVector operator - (const ComplexVector& V1, const ComplexVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = V1.RealComponents[i] - V2.RealComponents[i];
	  TmpComponents2[i] = V1.ImaginaryComponents[i] - V2.ImaginaryComponents[i];
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// substract two vectors with left one real
//
// V1 = first vector (real)
// V2 = second vector
// return value = resulting vector

ComplexVector operator - (const RealVector& V1, const ComplexVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = V1.Components[i] - V2.RealComponents[i];
	  TmpComponents2[i] = -V2.ImaginaryComponents[i];
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// substract two vectors with rightt one real
//
// V1 = first vector 
// V2 = second vector (real)
// return value = resulting vector

ComplexVector operator - (const ComplexVector& V1, const RealVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = V1.RealComponents[i] - V2.Components[i];
	  TmpComponents2[i] = V1.ImaginaryComponents[i];
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// multiply a vector with a real number on the right hand side
//
// V1 = vector to multiply
// d = real to use
// return value = resulting vector

ComplexVector operator * (const ComplexVector& V1, double d) 
{
  if (V1.Dimension != 0)
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = V1.RealComponents[i]* d;
	  TmpComponents2[i] = V1.ImaginaryComponents[i] * d;
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// multiply a vector with a complex number on the right hand side
//
// V1 = vector to multiply
// d = complex to use
// return value = resulting vector

ComplexVector operator * (const ComplexVector& V1, const Complex& d) 
{
  if (V1.Dimension != 0)
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = V1.RealComponents[i]* d.Re - V1.ImaginaryComponents[i] * d.Im;
	  TmpComponents2[i] = V1.ImaginaryComponents[i] * d.Re + V1.RealComponents[i]* d.Im;
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// multiply a vector with a real number on the left hand side
//
// V1 = vector to multiply
// d = real to use
// return value = resulting vector

ComplexVector operator * (double d, const ComplexVector& V1) 
{
  if (V1.Dimension != 0)
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = V1.RealComponents[i]* d;
	  TmpComponents2[i] = V1.ImaginaryComponents[i] * d;
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// multiply a vector with a complex number on the left hand side
//
// V1 = vector to multiply
// d = complex to use
// return value = resulting vector

ComplexVector operator * (const Complex& d, const ComplexVector& V1) 
{
  if (V1.Dimension != 0)
    {
      double* TmpComponents1 = new double [V1.Dimension];
      double* TmpComponents2 = new double [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents1[i] = V1.RealComponents[i]* d.Re - V1.ImaginaryComponents[i] * d.Im;
	  TmpComponents2[i] = V1.ImaginaryComponents[i] * d.Re + V1.RealComponents[i]* d.Im;
	}
      return ComplexVector(TmpComponents1, TmpComponents2, V1.Dimension);
    }
  else
    return ComplexVector();
}

// multiply a vector with a real number on the right hand side
//
// d = real to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (double d) 
{
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] *= d;  
      this->ImaginaryComponents[i] *= d;  
    }
  return *this;
}

// divide a vector with a real number on the right hand side
//
// d = real to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator /= (double d) 
{
  d = 1.0 / d;
  for (int i = 0; i < this->Dimension; i++)
    {
      this->RealComponents[i] *= d;  
      this->ImaginaryComponents[i] *= d;  
    }
  return *this;
}

// multiply a vector with a complex number on the right hand side
//
// d = complex to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const Complex& d) 
{
  double tmp;
  for (int i = 0; i < this->Dimension; ++i)
    {
      tmp = d.Re * this->RealComponents[i] - d.Im * this->ImaginaryComponents[i];
      this->ImaginaryComponents[i] *= d.Re;
      this->ImaginaryComponents[i] += this->RealComponents[i] * d.Im;  
      this->RealComponents[i] = tmp;  
    }
  return *this;
}

// divide a vector with a complex number on the right hand side
//
// d = complex to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator /= (const Complex& d) 
{
  double fac = 1.0 / (d.Re * d.Re + d.Im * d.Im);
  double dRe = d.Re * fac;
  double dIm = - d.Im * fac;  
  double tmp;
  for (int i = 0; i < this->Dimension; i += 2)
    {
      tmp = dRe * this->RealComponents[i] - dIm * this->ImaginaryComponents[i];
      this->ImaginaryComponents[i] *= dRe;
      this->ImaginaryComponents[i] += this->RealComponents[i] * dIm;  
      this->RealComponents[i] = tmp;  
    }
  return *this;
}

// left multiply a vector with a real matrix (using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const RealMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrColumn))
    return *this;
  double* tmp1 = new double [M.NbrRow];
  double* tmp2 = new double [M.NbrRow];
  for (int i = 0; i < M.NbrRow; ++i)
    {
      tmp1[i] = 0.0;
      tmp2[i] = 0.0;
      for (int j = 0; j < M.NbrColumn; ++j)
	{
	  tmp1[i] += M.Columns[j].Components[i] * this->RealComponents[j];
	  tmp2[i] += M.Columns[j].Components[i] * this->ImaginaryComponents[j];
	}      
    }
  if ((this->RealComponents != 0) && (this->ImaginaryComponents != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealComponents;
	delete[] this->ImaginaryComponents;
      }
  this->Flag.Initialize();
  this->TrueDimension = M.NbrRow;
  this->Dimension = M.NbrRow;
  this->RealComponents = tmp1;
  this->ImaginaryComponents = tmp2;
  return *this;
}

// left multiply a vector with a complex matrix (using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const ComplexMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrColumn))
    return *this;
  double* tmp1 = new double [M.NbrRow];
  double* tmp2 = new double [M.NbrRow];
  for (int i = 0; i < M.NbrRow; ++i)
    {
      tmp1[i] = 0.0;
      tmp2[i] = 0.0;
      for (int j = 0; j < M.NbrColumn; ++j)
	{
	  tmp1[i] += (M.Columns[j].RealComponents[i] * this->RealComponents[j] - 
		      M.Columns[j].ImaginaryComponents[i] * this->ImaginaryComponents[j]);
	  tmp2[i] += (M.Columns[j].RealComponents[i] * this->ImaginaryComponents[j] + 
		      M.Columns[j].ImaginaryComponents[i] * this->RealComponents[j]);
	}      
    }
  if ((this->RealComponents != 0) && (this->ImaginaryComponents != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealComponents;
	delete[] this->ImaginaryComponents;
      }
  this->Flag.Initialize();
  this->TrueDimension = M.NbrRow;
  this->Dimension = M.NbrRow;
  this->RealComponents = tmp1;
  this->ImaginaryComponents = tmp2;
  return *this;
}

// left multiply a vector with an hermtian conjugated  complex matrix (using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator &= (const ComplexMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this;
  double* tmp1 = new double [M.NbrRow];
  double* tmp2 = new double [M.NbrRow];
  for (int i = 0; i < M.NbrColumn; i++)
    {
      tmp1[i] = 0.0;
      tmp2[i] = 0.0;
      for (int j = 0; j < M.NbrColumn; ++j)
	{
	  tmp1[i] += (M.Columns[j].RealComponents[i] * this->RealComponents[j] + 
		      M.Columns[j].ImaginaryComponents[i] * this->ImaginaryComponents[j]);
	  tmp2[i] += (M.Columns[j].RealComponents[i] * this->ImaginaryComponents[j] - 
		      M.Columns[j].ImaginaryComponents[i] * this->RealComponents[j]);
	}      
    }
  if ((this->RealComponents != 0) && (this->ImaginaryComponents != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealComponents;
	delete[] this->ImaginaryComponents;
      }
  this->Flag.Initialize();
  this->TrueDimension = M.NbrRow;
  this->Dimension = M.NbrRow;
  this->RealComponents = tmp1;
  this->ImaginaryComponents = tmp2;
  return *this;
}

// left multiply a vector with an hermitian matrix (using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const HermitianMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this;
  double* tmp1 = new double [M.NbrRow];
  double* tmp2 = new double [M.NbrRow];
  for (int i = 0; i < this->Dimension; ++i)
    {
      tmp1[i] = M.DiagonalElements[i] * this->RealComponents[i];
      tmp2[i] = M.DiagonalElements[i] * this->ImaginaryComponents[i];
      int j = 0;
      int pos = i - 1;
      for (; j < i; ++j)
	{
	  tmp1[i] += M.RealOffDiagonalElements[pos] * this->RealComponents[j] + 
	    M.ImaginaryOffDiagonalElements[pos] * this->ImaginaryComponents[j];
	  tmp2[i] += M.RealOffDiagonalElements[pos] * this->ImaginaryComponents[j] - 
	    M.ImaginaryOffDiagonalElements[pos] * this->RealComponents[j];
	  pos += this->Dimension - j - 2 + M.Increment;
	}
      ++pos;
      ++j;
      for (; j < this->Dimension; ++j)
	{
	  tmp1[i] += M.RealOffDiagonalElements[pos] * this->RealComponents[j] - 
	    M.ImaginaryOffDiagonalElements[pos] * this->ImaginaryComponents[j];
	  tmp2[i] += M.RealOffDiagonalElements[pos] * this->ImaginaryComponents[j] + 
	    M.ImaginaryOffDiagonalElements[pos] * this->RealComponents[j];
	  ++pos;
	}
    }
  if ((this->RealComponents != 0) && (this->ImaginaryComponents != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->RealComponents;
	delete[] this->ImaginaryComponents;
      }
  this->Flag.Initialize();
  this->TrueDimension = M.NbrRow;
  this->RealComponents = tmp1;
  this->ImaginaryComponents = tmp2;
  return *this;
}

// left multiply a vector with a complex tridiagonal hermitian matrix (without using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const ComplexTriDiagonalHermitianMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this;
  int ReducedDim = this->Dimension - 1;
  Complex Tmp1 (this->RealComponents[0], this->ImaginaryComponents[0]);
  Complex Tmp2;
  this->RealComponents[0] *= M.DiagonalElements[0];
  this->RealComponents[0] += this->RealComponents[1] * M.RealUpperDiagonalElements[0]
    - this->ImaginaryComponents[1] * M.ImaginaryUpperDiagonalElements[0];
  this->ImaginaryComponents[0] *= M.DiagonalElements[0];
  this->ImaginaryComponents[0] += this->ImaginaryComponents[1] * M.RealUpperDiagonalElements[0]
    + this->RealComponents[1] * M.ImaginaryUpperDiagonalElements[0];
  for (int i = 1; i < ReducedDim; i++)
    {
      Tmp2 = Complex (this->RealComponents[i], this->ImaginaryComponents[i]);
      this->RealComponents[i] *= M.DiagonalElements[i];
      this->RealComponents[i] += this->RealComponents[i + 1] * M.RealUpperDiagonalElements[i] 
	- this->ImaginaryComponents[i + 1] * M.ImaginaryUpperDiagonalElements[i] + Tmp1.Re * M.RealUpperDiagonalElements[i - 1] 
	+ Tmp1.Im * M.ImaginaryUpperDiagonalElements[i - 1];
      this->ImaginaryComponents[i] *= M.DiagonalElements[i];
      this->ImaginaryComponents[i] += this->ImaginaryComponents[i + 1] * M.RealUpperDiagonalElements[i]
	+ this->RealComponents[i + 1] * M.ImaginaryUpperDiagonalElements[i] + Tmp1.Im * M.RealUpperDiagonalElements[i - 1] 
	- Tmp1.Re * M.ImaginaryUpperDiagonalElements[i - 1];
      Tmp1 = Tmp2;
    }
  this->RealComponents[this->Dimension - 1] *= M.DiagonalElements[this->Dimension - 1];
  this->RealComponents[this->Dimension - 1] += Tmp1.Re * M.RealUpperDiagonalElements[this->Dimension - 2] 
    + Tmp1.Im * M.ImaginaryUpperDiagonalElements[this->Dimension - 2];
  this->ImaginaryComponents[this->Dimension - 1] *= M.DiagonalElements[this->Dimension - 1];  
  this->ImaginaryComponents[this->Dimension - 1] += Tmp1.Im * M.RealUpperDiagonalElements[this->Dimension - 2] 
    - Tmp1.Re * M.ImaginaryUpperDiagonalElements[this->Dimension - 2];
  return *this;
}

// left multiply a vector with a real tridiagonal symmetric matrix (without using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const RealTriDiagonalSymmetricMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this;
  int ReducedDim = this->Dimension - 1;
  Complex Tmp1 (this->RealComponents[0], this->ImaginaryComponents[1]);
  Complex Tmp2;
  this->RealComponents[0] *= M.DiagonalElements[0];
  this->RealComponents[0] += this->RealComponents[1] * M.UpperDiagonalElements[0];
  this->ImaginaryComponents[0] *= M.DiagonalElements[0];
  this->ImaginaryComponents[0] += this->ImaginaryComponents[1] * M.UpperDiagonalElements[0];
  for (int i = 1; i < ReducedDim; i++)
    {
      Tmp2 = Complex (this->RealComponents[i], this->ImaginaryComponents[i]);
      this->RealComponents[i] *= M.DiagonalElements[i];
      this->RealComponents[i] += this->RealComponents[i + 1] * M.UpperDiagonalElements[i] 
	+ Tmp1.Re * M.UpperDiagonalElements[i - 1]; 
      this->ImaginaryComponents[i] *= M.DiagonalElements[i];
      this->ImaginaryComponents[i] += this->ImaginaryComponents[i + 1] * M.UpperDiagonalElements[i]
	+ Tmp1.Im * M.UpperDiagonalElements[i - 1];
      Tmp1 = Tmp2;
    }
  this->RealComponents[this->Dimension - 1] *= M.DiagonalElements[this->Dimension - 1];
  this->RealComponents[this->Dimension - 1] += Tmp1.Re * M.UpperDiagonalElements[this->Dimension - 2];
  this->ImaginaryComponents[this->Dimension - 1] *= M.DiagonalElements[this->Dimension - 1];  
  this->ImaginaryComponents[this->Dimension - 1] += Tmp1.Im * M.UpperDiagonalElements[this->Dimension - 2];
  return *this;
}

// left multiply a vector with an hermitian matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const HermitianMatrix&  M, ComplexVector& V)
{
  if ((this->Dimension == 0) || (V.Dimension != M.NbrColumn) || (this->Dimension != M.NbrRow))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] = M.DiagonalElements[i] * V.RealComponents[i];
      this->ImaginaryComponents[i] = M.DiagonalElements[i] * V.ImaginaryComponents[i];
      int pos = i - 1;
      int j = 0;
      for (; j < i; ++j)
	{
	  this->RealComponents[i] += M.RealOffDiagonalElements[pos] * V.RealComponents[j] + 
	    M.ImaginaryOffDiagonalElements[pos ] * V.ImaginaryComponents[j];
	  this->ImaginaryComponents[i] += M.RealOffDiagonalElements[pos] * V.ImaginaryComponents[j] 
	    - M.ImaginaryOffDiagonalElements[pos] * V.RealComponents[j];
	  pos += this->Dimension - j - 2 + M.Increment;
	}
      ++pos;
      ++j;
      for (; j < this->Dimension; ++j)
	{
	  this->RealComponents[i] += M.RealOffDiagonalElements[pos] * V.RealComponents[j] - 
	    M.ImaginaryOffDiagonalElements[pos] * V.ImaginaryComponents[j];
	  this->ImaginaryComponents[i] += M.RealOffDiagonalElements[pos] * V.ImaginaryComponents[j] + 
	    M.ImaginaryOffDiagonalElements[pos] * V.RealComponents[j];
	  ++pos;
	}
    }
  return *this;
}

// do a partial left multication of a vector with an hermitian matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const HermitianMatrix&  M, ComplexVector& V, int sourceStart, 
					int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) || (this->Dimension != M.NbrRow))
    return *this;
  int Last = sourceStart + sourceNbrComponent;
  int Inc1 =  this->Dimension - sourceNbrComponent + M.Increment - 2;
  int Inc2 =  this->Dimension + M.Increment - 2;
  int j;
  double x;
  double y;
  int i = 0;
  int Pos = sourceStart - 1;
  for (; i < sourceStart; ++i)
    {
      x = 0.0;
      y = 0.0;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.RealComponents[j] - 
	    M.ImaginaryOffDiagonalElements[Pos] * V.ImaginaryComponents[j];
	  y += M.RealOffDiagonalElements[Pos] * V.ImaginaryComponents[j] 
	    + M.ImaginaryOffDiagonalElements[Pos] * V.RealComponents[j];
	  ++Pos;
	}
      Pos += Inc1 - i;
      this->RealComponents[i] = x;
      this->ImaginaryComponents[i] = y;
    }
  Inc1 = this->Dimension - Last + M.Increment;
  int Pos2 = Pos;
  int Pos3 = Pos;
  ++Pos2;
  for (; i < Last; ++i)
    {
      x = 0.0;
      y = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < i; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.RealComponents[j] + 
	    M.ImaginaryOffDiagonalElements[Pos] * V.ImaginaryComponents[j];
	  y += M.RealOffDiagonalElements[Pos] * V.ImaginaryComponents[j] 
	    - M.ImaginaryOffDiagonalElements[Pos] * V.RealComponents[j];
	  Pos += Inc2 - j;
	}
      x += M.DiagonalElements[i] * V.RealComponents[i];
      y += M.DiagonalElements[i] * V.ImaginaryComponents[i];
      ++j;
      for (; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos2] * V.RealComponents[j] - 
	    M.ImaginaryOffDiagonalElements[Pos2] * V.ImaginaryComponents[j];
	  y += M.RealOffDiagonalElements[Pos2] * V.ImaginaryComponents[j] + 
	    M.ImaginaryOffDiagonalElements[Pos2] * V.RealComponents[j];
	  ++Pos2;
	}
      Pos2 += Inc1; 
      ++Pos3;
      this->RealComponents[i] = x;
      this->ImaginaryComponents[i] = y;
    }
  for (; i < this->Dimension; ++i)
    {
      x = 0.0;
      y = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.RealComponents[j] + 
	    M.ImaginaryOffDiagonalElements[Pos] * V.ImaginaryComponents[j];
	  y += M.RealOffDiagonalElements[Pos] * V.ImaginaryComponents[j] - 
	    M.ImaginaryOffDiagonalElements[Pos] * V.RealComponents[j];
	  Pos += Inc2 - j;
	}
      ++Pos3;
      this->RealComponents[i] = x;
      this->ImaginaryComponents[i] = y;
   }
  return *this;
}

// left multiply a vector with an hermitian matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const HermitianMatrix&  M, ComplexVector& V)
{
  if ((this->Dimension == 0) || (V.Dimension != M.NbrColumn) || (this->Dimension != M.NbrRow))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] += M.DiagonalElements[i] * V.RealComponents[i];
      this->ImaginaryComponents[i] += M.DiagonalElements[i] * V.ImaginaryComponents[i];
      int pos = i - 1;
      int j = 0;
      for (; j < i; ++j)
	{
	  this->RealComponents[i] += M.RealOffDiagonalElements[pos] * V.RealComponents[j] + 
	    M.ImaginaryOffDiagonalElements[pos ] * V.ImaginaryComponents[j];
	  this->ImaginaryComponents[i] += M.RealOffDiagonalElements[pos] * V.ImaginaryComponents[j] 
	    - M.ImaginaryOffDiagonalElements[pos] * V.RealComponents[j];
	  pos += this->Dimension - j - 2 + M.Increment;
	}
      ++pos;
      ++j;
      for (; j < this->Dimension; ++j)
	{
	  this->RealComponents[i] += M.RealOffDiagonalElements[pos] * V.RealComponents[j] - 
	    M.ImaginaryOffDiagonalElements[pos] * V.ImaginaryComponents[j];
	  this->ImaginaryComponents[i] += M.RealOffDiagonalElements[pos] * V.ImaginaryComponents[j] + 
	    M.ImaginaryOffDiagonalElements[pos] * V.RealComponents[j];
	  ++pos;
	}
    }
  return *this;
}

// do a partial left multication of a vector with an hermitian matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const HermitianMatrix&  M, ComplexVector& V, int sourceStart, 
					   int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) || (this->Dimension != M.NbrRow))
    return *this;
  int Last = sourceStart + sourceNbrComponent;
  int Inc1 =  this->Dimension - sourceNbrComponent + M.Increment - 2;
  int Inc2 =  this->Dimension + M.Increment - 2;
  int j;
  double x;
  double y;
  int i = 0;
  int Pos = sourceStart - 1;
  for (; i < sourceStart; ++i)
    {
      x = 0.0;
      y = 0.0;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.RealComponents[j] - 
	    M.ImaginaryOffDiagonalElements[Pos] * V.ImaginaryComponents[j];
	  y += M.RealOffDiagonalElements[Pos] * V.ImaginaryComponents[j] 
	    + M.ImaginaryOffDiagonalElements[Pos] * V.RealComponents[j];
	  ++Pos;
	}
      Pos += Inc1 - i;
      this->RealComponents[i] += x;
      this->ImaginaryComponents[i] += y;
    }
  Inc1 = this->Dimension - Last + M.Increment;
  int Pos2 = Pos;
  int Pos3 = Pos;
  ++Pos2;
  for (; i < Last; ++i)
    {
      x = 0.0;
      y = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < i; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.RealComponents[j] + 
	    M.ImaginaryOffDiagonalElements[Pos] * V.ImaginaryComponents[j];
	  y += M.RealOffDiagonalElements[Pos] * V.ImaginaryComponents[j] 
	    - M.ImaginaryOffDiagonalElements[Pos] * V.RealComponents[j];
	  Pos += Inc2 - j;
	}
      x += M.DiagonalElements[i] * V.RealComponents[i];
      y += M.DiagonalElements[i] * V.ImaginaryComponents[i];
      ++j;
      for (; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos2] * V.RealComponents[j] - 
	    M.ImaginaryOffDiagonalElements[Pos2] * V.ImaginaryComponents[j];
	  y += M.RealOffDiagonalElements[Pos2] * V.ImaginaryComponents[j] + 
	    M.ImaginaryOffDiagonalElements[Pos2] * V.RealComponents[j];
	  ++Pos2;
	}
      Pos2 += Inc1; 
      ++Pos3;
      this->RealComponents[i] += x;
      this->ImaginaryComponents[i] += y;
    }
  for (; i < this->Dimension; ++i)
    {
      x = 0.0;
      y = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.RealComponents[j] + 
	    M.ImaginaryOffDiagonalElements[Pos] * V.ImaginaryComponents[j];
	  y += M.RealOffDiagonalElements[Pos] * V.ImaginaryComponents[j] - 
	    M.ImaginaryOffDiagonalElements[Pos] * V.RealComponents[j];
	  Pos += Inc2 - j;
	}
      ++Pos3;
      this->RealComponents[i] += x;
       this->ImaginaryComponents[i] += y;
   }
  return *this;
}

// left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const RealMatrix&  M, ComplexVector& V)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrColumn))
    return *this;
  if (V.Dimension != M.NbrRow)
    V.Resize(M.NbrRow);
  for (int i = 0; i < V.Dimension; ++i)
    {
      this->RealComponents[i] = 0.0;
      this->ImaginaryComponents[i] = 0.0;
      for (int j = 0; j < this->Dimension; ++j)
	{
	  this->RealComponents[i] += M.Columns[j].Components[i] * V.RealComponents[j];
	  this->ImaginaryComponents[i] += M.Columns[j].Components[i] * V.ImaginaryComponents[j];
	}
    }
  return *this;
}

// left multiply a vector with a complex matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const ComplexMatrix&  M, ComplexVector& V)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrColumn))
    return *this;
  if (V.Dimension != M.NbrRow)
    V.Resize(M.NbrRow);
  for (int i = 0; i < V.Dimension; ++i)
    {
      this->RealComponents[i] = 0.0;
      this->ImaginaryComponents[i] = 0.0;
      for (int j = 0; j < this->Dimension; ++j)
	{
	  this->RealComponents[i] += (M.Columns[j].RealComponents[i] * V.RealComponents[j] -  
				      M.Columns[j].ImaginaryComponents[i] * V.ImaginaryComponents[j]);
	  this->ImaginaryComponents[i] += (M.Columns[j].RealComponents[i] * V.ImaginaryComponents[j] + 
					   M.Columns[j].ImaginaryComponents[i] * V.RealComponents[j]);
	}
    }
  return *this;
}

// left multiply a vector with a matrix and use to store result in current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const Matrix&  M, ComplexVector& V)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->Multiply((HermitianMatrix&) M, V);
      break;
    case (Matrix::ComplexElements):
      return this->Multiply((ComplexMatrix&) M, V);
      break;
    case (Matrix::RealElements):
      return this->Multiply((RealMatrix&) M, V);
      break;
    default:
      return *this;
    }
  return *this;
}

// left multiply a vector with a matrix and add result to current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const Matrix&  M, ComplexVector& V)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->AddMultiply((HermitianMatrix&) M, V);
      break;
    default:
      return *this;
    }
  return *this;
}

// left multiply a vector with a matrix and use to store result in current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
					int destStart, int destStep)
{
  return *this;
}

// left multiply a vector with a matrix and add current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
					   int destStart, int destStep)
{
  return *this;
}

// left multiply a vector with a matrix and use to store result in current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceNbrComponent)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->Multiply((HermitianMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    default:
      return *this;
    }
  return *this;
}

// left multiply a vector with a matrix and add current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceNbrComponent)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->AddMultiply((HermitianMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    default:
      return *this;
    }
  return *this;
}

// left multiply a vector with a matrix and use to store result in current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
					int sourceNbrComponent, int destStart, int destStep)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->Multiply((HermitianMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
      break;
    default:
      return *this;
    }
  return *this;
}

// left multiply a vector with a matrix and add current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
					   int sourceNbrComponent, int destStart, int destStep)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->AddMultiply((HermitianMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
      break;
    default:
      return *this;
    }
  return *this;
}

// get vector norm
//
// return value = vector norm

double ComplexVector::Norm() 
{
  if (this->Dimension == 0)
    return 0.0;
  double tmp = this->RealComponents[0] * this->RealComponents[0] + this->ImaginaryComponents[0] * this->ImaginaryComponents[0];
  for (int i = 1; i  < this->Dimension; ++i)
    tmp += this->RealComponents[i] * this->RealComponents[i] + this->ImaginaryComponents[i] * this->ImaginaryComponents[i];
  return sqrt(tmp);
}
  
// get square of vector norm
//
// return value = square of vector norm

double ComplexVector::SqrNorm () 
{
  if (this->Dimension == 0)
    return 0.0;
  double tmp = this->RealComponents[0] * this->RealComponents[0] + this->ImaginaryComponents[0] * this->ImaginaryComponents[0];
  for (int i = 2; i  < this->Dimension; ++i)
    tmp += this->RealComponents[i] * this->RealComponents[i] + this->ImaginaryComponents[i] * this->ImaginaryComponents[i];
  return tmp;
}
  
// normalize vector
//
// return value = reference on current vector

ComplexVector& ComplexVector::Normalize()
{
  if (this->Dimension == 0)
    return *this;
  double tmp = this->RealComponents[0] * this->RealComponents[0] + this->ImaginaryComponents[0] * this->ImaginaryComponents[0];
  for (int i = 2; i  < this->Dimension; ++i)
    tmp += this->RealComponents[i] * this->RealComponents[i] + this->ImaginaryComponents[i] * this->ImaginaryComponents[i];
  tmp = 1.0 / sqrt(tmp);
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->RealComponents[i] *= tmp;
      this->ImaginaryComponents[i] *= tmp;
    }
  return *this;
}
  
// Extract a subvector from a given vector
//
// firstCoordinate = Coordinate where extraction has to begin
// lastCoordinate = Coordinate where extraction has to stop (extract also include this last coordinate)
// Step = distance to the next coordinate (1 means to take the following)
// return value = return corresponding subvector

ComplexVector ComplexVector::Extract(int firstCoordinate, int lastCoordinate, int step) 
{
  if (this->Dimension == 0)
    return ComplexVector();
  int TmpDimension = (lastCoordinate - firstCoordinate) / step;
  int TrueLast = (firstCoordinate + TmpDimension * step);
  TmpDimension++;
  double* TmpComponents1 = new double [TmpDimension];
  double* TmpComponents2 = new double [TmpDimension];
  int j = 0;
  for (int i = firstCoordinate; i <= TrueLast; i += step)
    {
      TmpComponents1[j] = this->RealComponents[i];
      TmpComponents2[j] = this->ImaginaryComponents[i];
      ++j;
    }
  return ComplexVector(TmpComponents1, TmpComponents2, TmpDimension);  
}
  
// Merge a subvector into a given vector
//
// V = vector to merge
// firstCoordinate = Coordinate where merge has to begin
// step = distance to the next coordinate in the destination vector (1 means to take the following)
// return value = reference to the current Vector

ComplexVector& ComplexVector::Merge(const ComplexVector& V, int firstCoordinate, int step) 
{
  if ((this->Dimension == 0) || (V.Dimension == 0))
    return *this;
  int Max = firstCoordinate + (V.Dimension * step);
  if (Max > this->Dimension)
    return *this;
  int  j = 0;
  for (int i = firstCoordinate; i < Max; i += step)
    {
      this->RealComponents[i] = V.RealComponents[j];
      this->ImaginaryComponents[i] = V.ImaginaryComponents[j];      
      ++j;
    }
  return *this;
}
  
// Merge a real subvector into a complex given vector
//
// V = real vector to merge
// firstCoordinate = Coordinate where merge has to begin
// step = distance to the next coordinate in the destination vector (1 means to take the following)
// return value = reference to the current Vector

ComplexVector& ComplexVector::Merge(const RealVector& V, int firstCoordinate, int step) 
{
  if ((this->Dimension == 0) || (V.Dimension == 0))
    return *this;
  int Max = firstCoordinate + (V.Dimension * step);
  if (Max > this->Dimension)
    return *this;
  int  j = 0;
  for (int i = firstCoordinate; i < Max; i += step)
    {
      this->RealComponents[i] = V.Components[j];
      this->ImaginaryComponents[i] = 0.0;      
      ++j;
    }
  return *this;
}
  
// Output Stream overload
//

ostream& operator << (ostream& Str, const ComplexVector& P)
{
  for (int i = 0; i < P.Dimension; ++i)
    {
      Str << P.RealComponents[i];
      if (P.ImaginaryComponents[i] < 0.0)
	Str << P.ImaginaryComponents[i] << "i    ";
      else
	if (P.ImaginaryComponents[i] != 0.0)
	  Str << "+" << P.ImaginaryComponents[i] << "i    ";
	else
	  Str << "    ";
      Str << endl;
    }
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// v = vector to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexVector& v)
{
  Str << "{";
  int i = 0;
  for (; i < (v.Dimension - 1); ++i)
    {
      Str << v.RealComponents[i];
      if (v.ImaginaryComponents[i] < 0.0)
	Str << v.ImaginaryComponents[i] << "I    ";
      else
	if (v.ImaginaryComponents[i] != 0.0)
	  Str << "+" << v.ImaginaryComponents[i] << "I,";
	else
	  Str << ",";
    }
  Str << v.RealComponents[i++];
  if (v.ImaginaryComponents[i] < 0.0)
    Str << v.ImaginaryComponents[i] << "I}";
  else
    if (v.ImaginaryComponents[i] != 0.0)
      Str << "+" << v.ImaginaryComponents[i] << "I}";
    else
      Str << "}";
  return Str;
}

// write vector in a file 
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool ComplexVector::WriteVector (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  File.write ((char*) &(this->Dimension), sizeof(int));
  for (int i = 0; i < this->Dimension; ++i)
    {
      File.write ((char*) (&(this->RealComponents[i])), sizeof(double));
      File.write ((char*) (&(this->ImaginaryComponents[i])), sizeof(double));
    }
  File.close();
  return true;
}

// write vector in a file in ascii mode
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool ComplexVector::WriteAsciiVector (char* fileName)
{
  ofstream File;
  File.precision(14);
  File.open(fileName, ios::binary | ios::out);
  int ReducedDimension = this->Dimension - 1;
  for (int i = 0; i < ReducedDimension; ++i)
    File << this->RealComponents[i] << " " << this->ImaginaryComponents[i] << " ";
  File << this->RealComponents[ReducedDimension] << " " << this->ImaginaryComponents[ReducedDimension] << endl;  
  File.close();
  return true;
}

// read vector from a file 
//
// fileName = name of the file where the vector has to be read
// return value = true if no error occurs

bool ComplexVector::ReadVector (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  int TmpDimension;
  File.read ((char*) &(TmpDimension), sizeof(int));
  this->Resize(TmpDimension);
  for (int i = 0; i < this->Dimension; ++i)
    {
      File.read ((char*) (&(this->RealComponents[i])), sizeof(double));
      File.read ((char*) (&(this->ImaginaryComponents[i])), sizeof(double));
    }
  File.close();
  return true;
}

