////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                a set of functions usefull for integer algebra              //
//                                                                            //
//                        last modification : 11/09/2003                      //
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


#include "MathTools/IntegerAlgebraTools.h"


// find greatest common divider
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD

int FindGCD(int m, int n)
{
  if (m < n)
    return RecursiveFindGCD (m, n);
  else
    return RecursiveFindGCD (n, m);
  return n;
}

// find greatest common divider (recurisive part of the method)
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD

int RecursiveFindGCD(int m, int n)
{
  if (m == 0)
    return n;
  else
    return FindGCD ((n % m), m);
}

