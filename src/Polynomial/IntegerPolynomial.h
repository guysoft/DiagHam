////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   Class of polynomial with integer coefficients            //
//                                                                            //
//                        last modification : 27/05/2005                      //
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


#ifndef INTEGERPOLYNOMIAL_H
#define INTEGERPOLYNOMIAL_H


#include "config.h"
#include "Complex.h"

#include <iostream>


using std::ostream;


class IntegerPolynomial
{

private:

  // array where polynomial coefficients are stored
  long* Coefficient;
  
  // array where roots are stored
  Complex* Root;

  // polynomial degree
  int Degree;

  // number of found roots
  int NbrRoot;

  bool RootFlag;
  bool NoRootFlag;

public:

  // default constructor
  //
  IntegerPolynomial ();
 
  // constructor
  //
  // degree = polynomial degree
  // coefficients = coefficients array ( first element is associated to the -power term)
  // flag = true if coefficients array has to be used directly and not duplicated
  IntegerPolynomial (int degree, long* coefficients, bool flag = true);

  // copy constructor
  //
  // P = polynomial to copy
  IntegerPolynomial (const IntegerPolynomial& P);
  
  // destructor
  //
  ~IntegerPolynomial();

  // assignement
  //
  // P = polynomial to copy
  // return value = reference on current polynomial 
  IntegerPolynomial& operator = (const IntegerPolynomial& P);

  // Return polynomial value at a given point
  //
  // x = point where to evaluate polynomial
  // return value = polynomial value at x
  double PolynomialEvaluate (double x);
  
  // Return polynomial value at a given point
  //
  // x = point where to evaluate polynomial
  // return value = polynomial value at x
  Complex PolynomialEvaluate (Complex x);
  
  // Evaluate polynomial derivative   
  //
  // x = position where to evaluate polynomial derivative
  // return value = polynomial derivative at x
  double DerivativeEvaluate (double x);

  // Evaluate polynomial derivative   
  //
  // x = position where to evaluate polynomial derivative
  // return value = polynomial derivative at x
  Complex DerivativeEvaluate (Complex x);

  // Evaluate polynomial n-th derivative 
  //
  // x = position where to evaluate polynomial n-th derivative
  // n = derivative order
  // return value = polynomial n-th derivative at x  
  double DerivativeEvaluate (double x, int n);

  // Evaluate polynomial n-th derivative 
  //
  // x = position where to evaluate polynomial n-th derivative
  // n = derivative order
  // return value = polynomial n-th derivative at x  
  Complex DerivativeEvaluate (Complex x, int n);

  // Return Derivative of the polynomial 
  //
  // return value = polynomial derivative
  IntegerPolynomial DerivatePolynomial ();
  
  // Evaluate polynomial n-th derivative
  //
  // n = derivative order
  // return value = polynomial n-th derivative
  IntegerPolynomial DerivatePolynomial (int n);

  // Find Roots of the polynomial
  void SolvePolynomial();
  
  // get the number of roots
  //
  // return value = number of roots
  int GetNbrRoots ();
  
  // get the degree of the polynomial
  //
  // return value = degree of the polynomial
  int GetPolynomialDegree ();
  
  // return coefficient corresponding to the nth-degree
  //
  // n = index of the corresponding coefficient
  // return value = coefficient
  long PolynomialCoefficient (int n);

  // Return the n-th root (with root(n-1) <= root(n))
  //
  // n = index of the root (0 being root with the lowest modulus)
  // return value = n-th root
  Complex PolynomialRoot (int n);

  // refine root value 
  void RefineRoot (double Epsilon, int MaxIter);

  // arithmetic operators
  friend IntegerPolynomial operator - (const IntegerPolynomial& P);
  friend IntegerPolynomial operator + (const IntegerPolynomial& P1, const IntegerPolynomial& P2);
  friend IntegerPolynomial operator - (const IntegerPolynomial& P1, const IntegerPolynomial& P2);
  friend IntegerPolynomial operator * (const IntegerPolynomial& P, const long& d);
  friend IntegerPolynomial operator * (const long& d, const IntegerPolynomial& P);
  friend IntegerPolynomial operator * (const IntegerPolynomial& P1, const IntegerPolynomial& P2);
  friend IntegerPolynomial operator / (const IntegerPolynomial& P, const double& d);
  friend IntegerPolynomial operator / (const IntegerPolynomial& P1, const IntegerPolynomial& P2); // Degree P1 > Degree P2
  friend IntegerPolynomial operator % (const IntegerPolynomial& P1, const IntegerPolynomial& P2); // Degree P1 > Degree P2
  IntegerPolynomial& operator += (const IntegerPolynomial& P);
  IntegerPolynomial& operator -= (const IntegerPolynomial& P);
  IntegerPolynomial& operator *= (const long& d);
  IntegerPolynomial& operator *= (const IntegerPolynomial& P);

  //Output Stream overload
  friend ostream& operator << (ostream& Str, const IntegerPolynomial& P);
  
private:

  // Find Root of a linear polynomial
  void SolveLinear();

  // Sort roots of a polynomial
  //
  void SortRoots ();

};

// return number of roots
//
// return value = number of roots

inline int IntegerPolynomial::GetNbrRoots ()
{
  return this->NbrRoot;
}
  
// get the degree of the polynomial
//
// return value = degree of the polynomial

inline int IntegerPolynomial::GetPolynomialDegree () 
{
  return this->Degree;
}
  
// return coefficient corresponding to the nth-degree
//
// n = index of the corresponding coefficient
// return value = coefficient

inline long IntegerPolynomial::PolynomialCoefficient (int n)
{
  return this->Coefficient[n];
}

// Return the n-th root (with root(n-1) <= root(n))
//
// n = index of the root (0 being root with the lowest modulus)
// return value = n-th root

inline Complex IntegerPolynomial::PolynomialRoot (int n) 
{
  return this->Root[n];
}

#endif

