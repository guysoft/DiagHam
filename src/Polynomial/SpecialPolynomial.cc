////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                           RayTrace version  0.10                           //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   Class of polynomial with real coefficients               //
//                                                                            //
//                        last modification : 15/01/2001                      //
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

#include "SpecialPolynomial.h"
#include "MathTools/FactorialCoefficient.h"

using std::ostream;

// return a Laguerre Polynomial of rank n (optionally associated Laguerre function)
// n = index
// alpha = modifier
Polynomial LaguerrePolynomial(int n, int alpha)
{
  double *Coefficients;
  if (n<0)
    {
      Coefficients=new double[1];
      Coefficients[0]=0.0;
      n=0;
    }
  else
    {
      Coefficients=new double[n+1];
      FactorialCoefficient MyCoefficient; 
      for (int m=0; m<=n; ++m)
	{
	  MyCoefficient.SetToOne();
	  MyCoefficient.FactorialMultiply(n+alpha);
	  MyCoefficient.FactorialDivide(m);
	  MyCoefficient.FactorialDivide(alpha+m);
	  MyCoefficient.FactorialDivide(n-m);
	  Coefficients[m]=MyCoefficient.GetNumericalValue();
	  if (m&1)
	    Coefficients[m]*=-1.0;
	}
    }
  return Polynomial (n, Coefficients, true);
}

