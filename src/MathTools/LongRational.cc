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


#include "config.h"
#include "MathTools/LongRational.h"
#include "GeneralTools/Endian.h"

#include <iostream>


using std::cout;
using std::endl;

#ifdef __GMP__


// default constructor 
//

LongRational::LongRational()
{
  mpq_init(this->Value);
}

// constructor from an integer
//
// x = value to assign to the rational coefficient

LongRational::LongRational(long x)  
{
  mpq_init(this->Value);
  mpq_set_si (this->Value, x ,1ul);
}

// constructor from a rational number
//
// x = numerator to assign to the rational coefficient
// y = denominator to assign to the rational coefficient

LongRational::LongRational(long x, long y)  
{
  mpq_init(this->Value);
  if (y >= 0l)
    mpq_set_si (this->Value, x ,(unsigned long) y);
  else
    mpq_set_si (this->Value, -x ,(unsigned long) (-y));  
  mpq_canonicalize(this->Value);
}

// copy constructor from a rational number
//
// x =  rational coefficient to copy

LongRational::LongRational(const LongRational& x)
{
  mpq_init(this->Value);
  mpq_set (this->Value, x.Value);
}

// destructor
//

LongRational::~LongRational()
{
  mpq_clear(this->Value);
}

// assignement
//
// rational = rational coefficient to assign
// return value = reference on current rational coefficient

LongRational& LongRational::operator = (const LongRational& rational)
{
  mpq_set(this->Value, rational.Value);
  return *this;
}

// assignement from integer number
//
// x = interger to assign
// return value = reference on current rational coefficient

LongRational& LongRational::operator = (long x)
{
  mpq_set_si (this->Value, x ,1ul);
  return *this;
}

// set the coefficient to one
//
// return value = reference on current coefficient

LongRational& LongRational::SetToOne()
{
  mpq_set_si (this->Value, 1l ,1ul);
  return *this;
}

// add two rational numbers
//
// x = first rational
// y = second rational
// return value = sum

LongRational operator + (const LongRational& x, const LongRational& y)
{
  LongRational Tmp;
  mpq_add(Tmp.Value, x.Value, y.Value);
  return Tmp;
}

// add a rational number and an integer
//
// x = rational number
// y = integer
// return value = sum

LongRational operator + (const LongRational& x, long y)
{
  LongRational Tmp;
  mpq_t Tmp2;
  mpq_init(Tmp2);
  mpq_set_si (Tmp2, y ,1ul);
  mpq_add(Tmp.Value, x.Value, Tmp2);  
  mpq_clear(Tmp2);  
  return Tmp;
}

// add a rational number and an integer
//
// y = rational number
// x = integer
// return value = sum

LongRational operator + (long y, const LongRational& x)
{
  LongRational Tmp;
  mpq_t Tmp2;
  mpq_init(Tmp2);
  mpq_set_si (Tmp2, y ,1ul);
  mpq_add(Tmp.Value, x.Value, Tmp2);  
  mpq_clear(Tmp2);  
  return Tmp;
}

// substract an rational number to another
//
// x = first rational
// y = second rational
// return value = substraction

LongRational operator - (const LongRational& x, const LongRational& y)
{
  LongRational Tmp;
  mpq_sub(Tmp.Value, x.Value, y.Value);
  return Tmp;
}

// substract an integer to a rational number
//
// x = rational number
// y = integer
// return value = substraction

LongRational operator - (const LongRational& x, long y)
{
  LongRational Tmp;
  mpq_t Tmp2;
  mpq_init(Tmp2);
  mpq_set_si (Tmp2, y ,1ul);
  mpq_sub(Tmp.Value, x.Value, Tmp2);  
  mpq_clear(Tmp2);  
  return Tmp;
}

// substract a rational number to an integer
//
// y = integer
// x = rational number
// return value = substraction

LongRational operator - (long y, const LongRational& x)
{
  LongRational Tmp;
  mpq_t Tmp2;
  mpq_init(Tmp2);
  mpq_set_si (Tmp2, y ,1ul);
  mpq_sub(Tmp.Value, x.Value, Tmp2);  
  mpq_clear(Tmp2);  
  mpq_neg(Tmp.Value, Tmp.Value);
  return Tmp;
}

// multiply two rational numbers
//
// x = first rational
// y = second rational
// return value = product

LongRational operator * (const LongRational& x, const LongRational& y)
{
  LongRational Tmp;
  mpq_mul(Tmp.Value, x.Value, y.Value);
  return Tmp;
}

// multiply a rational number and an integer
//
// x = rational number
// y = integer
// return value = product

LongRational operator * (const LongRational& x, long y)
{
  LongRational Tmp;
  mpq_t Tmp2;
  mpq_init(Tmp2);
  mpq_set_si (Tmp2, y ,1ul);
  mpq_mul(Tmp.Value, x.Value, Tmp2);  
  mpq_clear(Tmp2);  
  return Tmp;
}

// multiply a rational number and an integer
//
// y = rational number
// x = integer
// return value = product

LongRational operator * (long y, const LongRational& x)
{
  LongRational Tmp;
  mpq_t Tmp2;
  mpq_init(Tmp2);
  mpq_set_si (Tmp2, y ,1ul);
  mpq_mul(Tmp.Value, x.Value, Tmp2);  
  mpq_clear(Tmp2);  
  return Tmp;
}

// divide two rational numbers
//
// x = first rational
// y = second rational
// return value = division

LongRational operator / (const LongRational& x, const LongRational& y)
{
  LongRational Tmp;
  mpq_div(Tmp.Value, x.Value, y.Value);
  return Tmp;
}

// divide a rational number ny an integer
//
// x = rational number
// y = integer
// return value = division

LongRational operator / (const LongRational& x, long y)
{
  LongRational Tmp;
  mpq_t Tmp2;
  mpq_init(Tmp2);
  mpq_set_si (Tmp2, y ,1ul);
  mpq_div(Tmp.Value, x.Value, Tmp2);  
  mpq_clear(Tmp2);  
  return Tmp;
}

// divide an integer ny a rational number
//
// y = rational number
// x = integer
// return value = division

LongRational operator / (long y, const LongRational& x)
{
  LongRational Tmp;
  mpq_t Tmp2;
  mpq_init(Tmp2);
  mpq_set_si (Tmp2, y ,1ul);
  mpq_div(Tmp.Value, Tmp2, x.Value);  
  mpq_clear(Tmp2);  
  return Tmp;
}

// add a rational
//
// x = rational to use
// return value = reference on current coefficient

LongRational& LongRational::operator += (const  LongRational& x)
{
  mpq_add(this->Value, this->Value, x.Value);  
  return *this;
}

// add an integer
//
// x = integer to use
// return value = reference on current coefficient

LongRational& LongRational::operator += (long x)
{
  mpq_t Tmp;
  mpq_init(Tmp);
  mpq_set_si (Tmp, x, 1ul);
  mpq_add(this->Value, this->Value, Tmp);  
  mpq_clear(Tmp);  
  return *this;
}
 
// substract a rational
//
// x = rational to use
// return value = reference on current coefficient

LongRational& LongRational::operator -= (const  LongRational& x)
{
  mpq_sub(this->Value, this->Value, x.Value);  
  return *this;
}

// substract an integer
//
// x = integer to use
// return value = reference on current coefficient

LongRational& LongRational::operator -= (long x)
{
  mpq_t Tmp;
  mpq_init(Tmp);
  mpq_set_si (Tmp, x, 1ul);
  mpq_sub(this->Value, this->Value, Tmp);  
  mpq_clear(Tmp);  
  return *this;
}

// multiply by an integer
//
// x = integer to use
// return value = reference on current coefficient

LongRational& LongRational::operator *= (long x)
{
  mpq_t Tmp;
  mpq_init(Tmp);
  mpq_set_si (Tmp, x, 1ul);
  mpq_mul(this->Value, this->Value, Tmp);  
  mpq_clear(Tmp);  
  return *this;
}

// multiply by a rational
//
// x = rational to use
// return value = reference on current coefficient

LongRational& LongRational::operator *= (const LongRational& x)
{
  mpq_mul(this->Value, this->Value, x.Value);  
  return *this;
}

// divide by an integer
//
// y = integer to use
// return value = reference on current coefficient

LongRational& LongRational::operator /= (long y)
{
  mpq_t Tmp;
  mpq_init(Tmp);
  mpq_set_si (Tmp, y, 1ul);
  mpq_div(this->Value, this->Value, Tmp);  
  mpq_clear(Tmp);  
  return *this;
}

// divide by a rational
//
// x = rational to use
// return value = reference on current coefficient

LongRational& LongRational::operator /= (const LongRational& x)
{
  mpq_div(this->Value, this->Value, x.Value);  
  return *this;
}

// Output stream overload
//
// str = reference on output stream
// x = rational to print
// return value = reference on output stream

ostream& operator << (ostream& str, LongRational& x)
{
  str << x.Value;
  return str;
}


// write a rational in binary mode
//
// file = reference on the output file stream
// return value = reference on the output file stream

ofstream& LongRational::Write(ofstream& file)
{
  file << this->Value << endl;
  return file;
}

// read a rational in binary mode
//
// file = reference on the output file stream
// return value = reference on the output file stream

ifstream& LongRational::Read(ifstream& file)
{
  file >>  this->Value;
  return file;
}

#else

// default constructor 
//

LongRational::LongRational()
{
  this->Numerator = 0l;
  this->Denominator = 1l;
}

// constructor from an integer
//
// x = value to assign to the rational coefficient

LongRational::LongRational(long x)  
{
  this->Numerator = x;
  this->Denominator = 1l;
}

// constructor from a rational number
//
// x = numerator to assign to the rational coefficient
// y = denominator to assign to the rational coefficient

LongRational::LongRational(long x, long y)  
{
  this->Numerator = x;
  this->Denominator = y;
}

// copy constructor from a rational number
//
// x =  rational coefficient to copy

LongRational::LongRational(const LongRational& x)
{
  this->Numerator = x.Numerator;
  this->Denominator = x.Denominator;
}

// destructor
//

LongRational::~LongRational()
{
}

// assignement
//
// rational = rational coefficient to assign
// return value = reference on current rational coefficient

LongRational& LongRational::operator = (const LongRational& rational)
{
  this->Numerator = rational.Numerator;
  this->Denominator = rational.Denominator;
  return *this;
}

// assignement from integer number
//
// x = interger to assign
// return value = reference on current rational coefficient

LongRational& LongRational::operator = (long x)
{
  this->Numerator = x;
  this->Denominator = 1.0;
  return *this;
}

// set the coefficient to one
//
// return value = reference on current coefficient

LongRational& LongRational::SetToOne()
{
  this->Numerator = 1l;
  this->Denominator = 1l;
  return *this;
}

// add two rational numbers
//
// x = first rational
// y = second rational
// return value = sum

LongRational operator + (const LongRational& x, const LongRational& y)
{
  LongRational Tmp;
  Tmp.Numerator = (x.Numerator * y.Denominator) + (y.Numerator * x.Denominator);
  Tmp.Denominator = x.Denominator * y.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// add a rational number and an integer
//
// x = rational number
// y = integer
// return value = sum

LongRational operator + (const LongRational& x, long y)
{
  LongRational Tmp;
  Tmp.Numerator = x.Numerator + (y * x.Denominator);
  Tmp.Denominator = x.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// add a rational number and an integer
//
// y = rational number
// x = integer
// return value = sum

LongRational operator + (long y, const LongRational& x)
{
  LongRational Tmp;
  Tmp.Numerator = x.Numerator + (y * x.Denominator);
  Tmp.Denominator = x.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// substract an rational number to another
//
// x = first rational
// y = second rational
// return value = substraction

LongRational operator - (const LongRational& x, const LongRational& y)
{
  LongRational Tmp;
  Tmp.Numerator = (x.Numerator * y.Denominator) - (y.Numerator * x.Denominator);
  Tmp.Denominator = x.Denominator * y.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// substract an integer to a rational number
//
// x = rational number
// y = integer
// return value = substraction

LongRational operator - (const LongRational& x, long y)
{
  LongRational Tmp;
  Tmp.Numerator = x.Numerator - (y * x.Denominator);
  Tmp.Denominator = x.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// substract a rational number to an integer
//
// y = integer
// x = rational number
// return value = substraction

LongRational operator - (long y, const LongRational& x)
{
  LongRational Tmp;
  Tmp.Numerator = (y * x.Denominator) - x.Numerator;
  Tmp.Denominator = x.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// multiply two rational numbers
//
// x = first rational
// y = second rational
// return value = product

LongRational operator * (const LongRational& x, const LongRational& y)
{
  LongRational Tmp = x;
  long Tmp2 = y.Numerator;
  long Tmp3 = Tmp.FindGCD(Tmp2, Tmp.Denominator);
  long Tmp4 = y.Denominator;
  long Tmp5 = Tmp.FindGCD(Tmp4, Tmp.Numerator);
  Tmp.Denominator /= Tmp3;
  Tmp2 /= Tmp3;
  Tmp.Numerator /= Tmp5;
  Tmp4 /= Tmp5;
  Tmp.Numerator *= Tmp2;
  Tmp.Denominator *= Tmp4;
  return Tmp;
}

// multiply a rational number and an integer
//
// x = rational number
// y = integer
// return value = product

LongRational operator * (const LongRational& x, long y)
{
  LongRational Tmp = x;
  long Tmp2 = Tmp.FindGCD(y, Tmp.Denominator);
  Tmp.Denominator /= Tmp2;
  y /= Tmp2;
  Tmp.Numerator *= y;
  return Tmp;
}

// multiply a rational number and an integer
//
// y = rational number
// x = integer
// return value = product

LongRational operator * (long y, const LongRational& x)
{
  LongRational Tmp = x;
  long Tmp2 = Tmp.FindGCD(y, Tmp.Denominator);
  Tmp.Denominator /= Tmp2;
  y /= Tmp2;
  Tmp.Numerator *= y;
  return Tmp;
}

// divide two rational numbers
//
// x = first rational
// y = second rational
// return value = division

LongRational operator / (const LongRational& x, const LongRational& y)
{
  LongRational Tmp = x;
  long Tmp2 = y.Denominator;
  long Tmp3 = Tmp.FindGCD(Tmp2, Tmp.Denominator);
  long Tmp4 = y.Numerator;
  long Tmp5 = Tmp.FindGCD(Tmp4, Tmp.Numerator);
  Tmp.Denominator /= Tmp3;
  Tmp2 /= Tmp3;
  Tmp.Numerator /= Tmp5;
  Tmp4 /= Tmp5;
  Tmp.Numerator *= Tmp2;
  Tmp.Denominator *= Tmp4;
  return Tmp;
}

// divide a rational number ny an integer
//
// x = rational number
// y = integer
// return value = division

LongRational operator / (const LongRational& x, long y)
{
  LongRational Tmp (x);
  long Tmp2 = Tmp.FindGCD(y, Tmp.Numerator);
  Tmp.Numerator /= Tmp2;
  y /= Tmp2;
  Tmp.Denominator *= y;
  return Tmp;
}

// divide an integer ny a rational number
//
// y = rational number
// x = integer
// return value = division

LongRational operator / (long y, const LongRational& x)
{
  LongRational Tmp;
  Tmp.Denominator = x.Numerator;
  Tmp.Numerator = x.Denominator;
  long Tmp2 = Tmp.FindGCD(y, Tmp.Denominator);
  Tmp.Denominator /= Tmp2;
  y /= Tmp2;
  Tmp.Numerator *= y;
  return Tmp;
}

// add a rational
//
// x = rational to use
// return value = reference on current coefficient

LongRational& LongRational::operator += (const  LongRational& x)
{
  this->Numerator *= x.Denominator;
  this->Numerator += x.Numerator * this->Denominator;
  this->Denominator *= x.Denominator;
  this->Simplify();  
  return *this;
}

// add an integer
//
// x = integer to use
// return value = reference on current coefficient

LongRational& LongRational::operator += (long x)
{
  this->Numerator += x * this->Denominator;
  this->Simplify();
  return *this;
}
 
// substract a rational
//
// x = rational to use
// return value = reference on current coefficient

LongRational& LongRational::operator -= (const  LongRational& x)
{
  this->Numerator *= x.Denominator;
  this->Numerator -= x.Numerator * this->Denominator;
  this->Denominator *= x.Denominator;
  this->Simplify();  
  return *this;
}

// substract an integer
//
// x = integer to use
// return value = reference on current coefficient

LongRational& LongRational::operator -= (long x)
{
  this->Numerator -= x * this->Denominator;
  this->Simplify();
  return *this;
}

// multiply by an integer
//
// x = integer to use
// return value = reference on current coefficient

LongRational& LongRational::operator *= (long x)
{
  long Tmp = this->FindGCD(x, this->Denominator);
  x /= Tmp;
  this->Denominator /= Tmp;
  this->Numerator *= x;
  return *this;
}

// multiply by a rational
//
// x = rational to use
// return value = reference on current coefficient

LongRational& LongRational::operator *= (const LongRational& x)
{
  long Tmp = this->FindGCD(this->Numerator, x.Denominator);
  long Tmp2 = this->FindGCD(x.Numerator, this->Denominator);
  this->Numerator /= Tmp;
  this->Denominator /= Tmp2;
  this->Numerator *=  x.Numerator / Tmp2;
  this->Denominator *= x.Denominator / Tmp;
  return *this;     
}

// divide by an integer
//
// y = integer to use
// return value = reference on current coefficient

LongRational& LongRational::operator /= (long y)
{
  long Tmp = this->FindGCD(y, this->Denominator);
  y /= Tmp;
  this->Denominator *= y;
  this->Numerator /= Tmp;
  return *this;
}

// divide by a rational
//
// x = rational to use
// return value = reference on current coefficient

LongRational& LongRational::operator /= (const LongRational& x)
{
  long Tmp = this->FindGCD(this->Numerator, x.Numerator);
  long Tmp2 = this->FindGCD(x.Denominator, this->Denominator);
  this->Numerator /= Tmp;
  this->Denominator /= Tmp2;
  this->Numerator *=  x.Denominator / Tmp2;
  this->Denominator *= x.Numerator / Tmp;
  return *this;     
}

// Output stream overload
//
// str = reference on output stream
// x = rational to print
// return value = reference on output stream

ostream& operator << (ostream& str, LongRational& x)
{
   if (x.Denominator > 0l)
     {
       if (x.Denominator == 1l)
 	str << x.Numerator;
       else
 	str << x.Numerator << "/" << x.Denominator;
     }
   else
     {
       if (x.Denominator == -1l)
 	str << (-x.Numerator);
       else
 	str << (-x.Numerator) << "/" << (-x.Denominator);
     }
  return str;
}


// write a rational in binary mode
//
// file = reference on the output file stream
// return value = reference on the output file stream

ofstream& LongRational::Write(ofstream& file)
{
  WriteLittleEndian(file, this->Numerator);
  WriteLittleEndian(file, this->Denominator);
  return file;
}

// read a rational in binary mode
//
// file = reference on the output file stream
// return value = reference on the output file stream

ifstream& LongRational::Read(ifstream& file)
{
  ReadLittleEndian(file, this->Numerator);
  ReadLittleEndian(file, this->Denominator);
  return file;
}

#endif
