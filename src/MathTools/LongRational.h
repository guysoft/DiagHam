////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of rational numbers using long long               //
//                                                                            //
//                        last modification : 17/11/2010                      //
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


#ifndef LONGRATIONAL_H
#define LONGRATIONAL_H


#include "config.h"

#include <iostream>
#include <fstream>


using std::ostream;
using std::ofstream;
using std::ifstream;


class LongRational
{

 private:

  // numerator
  LONGLONG Numerator;
  // denominator
  LONGLONG Denominator;

 public:

  // default constructor 
  //
  LongRational();

  // constructor from an integer
  //
  // x = value to assign to the rational coefficient
  LongRational(long x);  

  // constructor from a rational number
  //
  // x = numerator to assign to the rational coefficient
  // y = denominator to assign to the rational coefficient
  LongRational(long x, long y);  

  // destructor
  //
  ~LongRational();

  // assignement
  //
  // rational = rational coefficient to assign
  // return value = reference on current rational coefficient
  LongRational& operator = (const LongRational& rational);

  // assignement from integer number
  //
  // x = interger to assign
  // return value = reference on current rational coefficient
  LongRational& operator = (long x);

  // set the coefficient to one
  //
  // return value = reference on current coefficient
  LongRational& SetToOne();

  // compute the opposit number
  //
  // x = first rational
  // return value = ooposit number
  friend LongRational operator - (const LongRational& x);

  // add two rational numbers
  //
  // x = first rational
  // y = second rational
  // return value = sum
  friend LongRational operator + (const LongRational& x, const LongRational& y);

  // add a rational number and an integer
  //
  // x = rational number
  // y = integer
  // return value = sum
  friend LongRational operator + (const LongRational& x, long y);

  // add a rational number and an integer
  //
  // y = rational number
  // x = integer
  // return value = sum
  friend LongRational operator + (long y, const LongRational& x);

  // substract an rational number to another
  //
  // x = first rational
  // y = second rational
  // return value = substraction
  friend LongRational operator - (const LongRational& x, const LongRational& y);

  // substract an integer to a rational number
  //
  // x = rational number
  // y = integer
  // return value = substraction
  friend LongRational operator - (const LongRational& x, long y);

  // substract a rational number to an integer
  //
  // y = integer
  // x = rational number
  // return value = substraction
  friend LongRational operator - (long y, const LongRational& x);

  // multiply two rational numbers
  //
  // x = first rational
  // y = second rational
  // return value = product
  friend LongRational operator * (const LongRational& x, const LongRational& y);

  // multiply a rational number and an integer
  //
  // x = rational number
  // y = integer
  // return value = product
  friend LongRational operator * (const LongRational& x, long y);

  // multiply a rational number and an integer
  //
  // y = rational number
  // x = integer
  // return value = product
  friend LongRational operator * (long y, const LongRational& x);

  // divide two rational numbers
  //
  // x = first rational
  // y = second rational
  // return value = division
  friend LongRational operator / (const LongRational& x, const LongRational& y);

  // divide a rational number by an integer
  //
  // x = rational number
  // y = integer
  // return value = division
  friend LongRational operator / (const LongRational& x, long y);

  // divide an integer by a rational number
  //
  // y = rational number
  // x = integer
  // return value = division
  friend LongRational operator / (long y, const LongRational& x);

  // add a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  LongRational& operator += (const  LongRational& x);

  // add an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  LongRational& operator += (long x);

  // substract a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  LongRational& operator -= (const  LongRational& x);

  // substract an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  LongRational& operator -= (long x);

  // multiply by an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  LongRational& operator *= (long x);

  // multiply by a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  LongRational& operator *= (const LongRational& x);

  // divide by an integer
  //
  // y = integer to use
  // return value = reference on current coefficient
  LongRational& operator /= (long y);

  // divide by a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  LongRational& operator /= (const LongRational& x);

  // test is two rational numbers are identical
  // 
  // x = first rational number
  // y = first rational number
  // return value = true if the two numbers are identical
  friend bool operator == (const LongRational& x, const LongRational& y);

  // test is two rational numbers are different
  // 
  // x = first rational number
  // y = first rational number
  // return value = true if the two numbers are different
  friend bool operator != (const LongRational& x, const LongRational& y);

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

  // Output stream overload
  //
  // str = reference on output stream
  // x = rational to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, LongRational& x);

  // write a rational in binary mode
  //
  // file = reference on the output file stream
  // return value = reference on the output file stream
  ofstream& Write(ofstream& file);

  // read a rational in binary mode
  //
  // file = reference on the output file stream
  // return value = reference on the output file stream
  ifstream& Read(ifstream& file);

 private:

  // Simplify the current rational number
  //
  void Simplify();

  // find greatest common divider
  //
  // m = first integer  
  // n = second integer (must be greater than m)
  // return value = GCD
  long FindGCD(long m, long n);


};

// return integer value associated to the coefficient numerator (0 if the coefficient can't be cast into an integer)
//
// return value = numerical coefficient

inline long LongRational::Num()
{
  return this->Numerator;
}

// return integer value associated to the coefficient denominator (0 if the coefficient can't be cast into an integer)
//
// return value = numerical value associated to the coefficient denominator

inline long LongRational::Den()
{
  return this->Denominator;
}

// return numerical value associated to the coefficient
//
// return value = numerical value associated to the coefficient

inline double LongRational::GetNumericalValue()
{
  return (((double) this->Numerator) / ((double) this->Denominator));
}

// Simplify the current rational number
//

inline void LongRational::Simplify()
{
  long Tmp = this->FindGCD(this->Numerator, this->Denominator);
  this->Numerator /= Tmp;
  this->Denominator /= Tmp;
}

// find greatest common divider
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD

inline long LongRational::FindGCD(long m, long n)
{
  long Tmp;
  while (n != 0l)
    {
      Tmp = n;
      n = m % n;
      m = Tmp;
    }
  return m;
}

// test is two rational numbers are identical
// 
// x = first rational number
// y = first rational number
// return value = true if the two numbers are identical

inline bool operator == (const LongRational& x, const LongRational& y)
{
  if (((x.Numerator == y.Numerator) && (x.Denominator == y.Denominator)) || ((x.Numerator == -y.Numerator) && (x.Denominator == -y.Denominator)))
    return true;
  else
    return false;    
}


// test is two rational numbers are different
// 
// x = first rational number
// y = first rational number
// return value = true if the two numbers are different

inline bool operator != (const LongRational& x, const LongRational& y)
{
  if (((x.Numerator != y.Numerator) || (x.Denominator != y.Denominator)) && ((x.Numerator != -y.Numerator) || (x.Denominator != -y.Denominator)))
    return true;
  else
    return false;    
}

// compute the opposit number
//
// x = first rational
// return value = ooposit number

inline LongRational operator - (const LongRational& x)
{
  LongRational Tmp(x);
  Tmp.Numerator *= -1l;
  return Tmp;
}


#endif


