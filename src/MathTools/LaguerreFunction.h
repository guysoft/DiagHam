////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of clebsch gordan coefficients                   //
//                                                                            //
//                        last modification : 19/06/2002                      //
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


#ifndef LAGUERREFUNCTION_H
#define LAGUERREFUNCTION_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Vector/RealVector.h"
#include <iostream>


using std::ostream;


class LaguerreFunction
{

 private:

  // degree
  int N;

  // modifier
  int Alpha;

  // array of prefactors of the polynomial terms
  double *Coefficients;

  
 public:

  // default constructor
  // creates L_0^0
  //
  LaguerreFunction();

  // constructor for a Laguerre polynomial
  //
  // n = index
  // alpha = modifier
  LaguerreFunction(int n, int alpha=0);

  // copy constructor (duplicating datas)
  //
  // coefficients = reference on function to copy
  LaguerreFunction (const LaguerreFunction& laguerre);

  // destructor
  //
  ~LaguerreFunction ();

  // get the value of the function for a given coordinate x
  //
  // x = complex coordinate
  // return = function value at z
  
  double operator()(const double &x);

  
  // get the value of the function for a given coordinate x
  //
  // x = argument
  // return = function value at x

  double GetValue(const double &x);
  
  // get the value of the function for a given coordinate z
  // values = complex vector where to store results
  // manyZ = complex vector of coordinates
  // return = function value at z
  void GetManyValues(RealVector &values, RealVector &manyX);

  // pretty-print a function value
  // str = stream to print to
  // x = point where to evaluate
  ostream& PrintValue(ostream &str, const double &x);
  
  
};

#endif
