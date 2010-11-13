////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of rational numbers                         //
//                                                                            //
//                        last modification : 12/11/2010                      //
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
#include "MathTools/Rational.h"

#include <iostream>


  // numerator
  long Numerator;
  // denominator
  long Denominator;

// default constructor 
//

Rational::Rational()
{
  this->Numerator = 0l;
  this->Denominator = 1l;
}

// constructor from an integer
//
// x = value to assign to the rational coefficient

Rational::Rational(long x)  
{
  this->Numerator = x;
  this->Denominator = 1l;
}

// constructor from a rational number
//
// x = numerator to assign to the rational coefficient
// y = denominator to assign to the rational coefficient

Rational::Rational(long x, long y)  
{
  this->Numerator = x;
  this->Denominator = y;
}

// destructor
//

Rational::~Rational()
{
}

// assignement
//
// rational = rational coefficient to assign
// return value = reference on current rational coefficient

Rational& Rational::operator = (Rational& rational)
{
  this->Numerator = rational.Numerator;
  this->Denominator = rational.Denominator;
  return *this;
}

// assignement from integer number
//
// x = interger to assign
// return value = reference on current rational coefficient

Rational& Rational::operator = (long x)
{
  this->Numerator = x;
  this->Denominator = 1.0;
  return *this;
}

// set the coefficient to one
//
// return value = reference on current coefficient

Rational& Rational::SetToOne()
{
  this->Numerator = 1l;
  this->Denominator = 1l;
  return *this;
}

// multiply by an integer
//
// x = integer to use
// return value = reference on current coefficient

Rational& Rational::operator *= (long x)
{
  long Tmp = this->FindGCD(x, this->Denominator);
  x /= Tmp;
  this->Denominator /= Tmp;
  this->Numerator *= x;
  return *this;
}

// divide by an integer
//
// y = integer to use
// return value = reference on current coefficient

Rational& Rational::operator /= (long y)
{
  long Tmp = this->FindGCD(y, this->Denominator);
  y /= Tmp;
  this->Denominator *= y;
  this->Numerator /= Tmp;
  return *this;
}


// find greatest common divider (recurisive part of the method)
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD

long Rational::RecursiveFindGCD(long m, long n)
{
  if (m == 0)
    return n;
  else
    return this->FindGCD ((n % m), m);
}
