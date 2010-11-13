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


#ifndef RATIONAL_H
#define RATIONAL_H


#include "config.h"

#include <iostream>


class Rational
{

 private:

  // numerator
  long Numerator;
  // denominator
  long Denominator;

 public:

  // default constructor 
  //
  Rational();

  // constructor from an integer
  //
  // x = value to assign to the rational coefficient
  Rational(long x);  

  // constructor from a rational number
  //
  // x = numerator to assign to the rational coefficient
  // y = denominator to assign to the rational coefficient
  Rational(long x, long y);  

  // destructor
  //
  ~Rational();

  // assignement
  //
  // rational = rational coefficient to assign
  // return value = reference on current rational coefficient
  Rational& operator = (Rational& rational);

  // assignement from integer number
  //
  // x = interger to assign
  // return value = reference on current rational coefficient
  Rational& operator = (long x);

  // set the coefficient to one
  //
  // return value = reference on current coefficient
  Rational& SetToOne();

  // multiply by an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  Rational& operator *= (long x);

  // divide by an integer
  //
  // y = integer to use
  // return value = reference on current coefficient
  Rational& operator /= (long y);

  // return integer value associated to the coefficient numerator (0 if the coefficient can't be cast into an integer)
  //
  // return value = numerical coefficient
  long Num();

  // return integer value associated to the coefficient denominator (0 if the coefficient can't be cast into an integer)
  //
  // return value = numerical value associated to the coefficient denominator
  long Den();

  // return numerical value associated to the coefficient
  //
  // return value = numerical value associated to the coefficient
  double GetNumericalValue();

 private:

  // find greatest common divider (recursive part of the method)
  //
  // m = first integer  
  // n = second integer (must be greater than m)
  // return value = GCD
  long FindGCD(long m, long n);

  // find greatest common divider (recurisive part of the method)
  //
  // m = first integer  
  // n = second integer (must be greater than m)
  // return value = GCD
  long RecursiveFindGCD(long m, long n);

};

// return integer value associated to the coefficient numerator (0 if the coefficient can't be cast into an integer)
//
// return value = numerical coefficient

inline long Rational::Num()
{
  return this->Numerator;
}

// return integer value associated to the coefficient denominator (0 if the coefficient can't be cast into an integer)
//
// return value = numerical value associated to the coefficient denominator

inline long Rational::Den()
{
  return this->Denominator;
}

// return numerical value associated to the coefficient
//
// return value = numerical value associated to the coefficient

inline double Rational::GetNumericalValue()
{
  return (((double) this->Numerator) / ((double) this->Denominator));
}

// find greatest common divider (recursive part of the method)
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD

inline long Rational::FindGCD(long m, long n)
{
  if (m < n)
    return this->RecursiveFindGCD (m, n);
  else
    return this->RecursiveFindGCD (n, m);
  return n;
}

#endif


