////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of function basis for particle on sphere             //
//                                                                            //
//                        last modification : 10/12/2002                      //
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
#include "FunctionBasis/ParticleOnCylinderFunctionBasis.h"
#include "Vector/RealVector.h"

#include <math.h>

#include <iostream>
using std::cout;
using std::endl;

// constructor
//
// maxMomentum = maximum momentum reached by a particle
// landauLevel = Landau level index
// ratio = aspect ratio of the cylinder
ParticleOnCylinderFunctionBasis::ParticleOnCylinderFunctionBasis(int maxMomentum, int landauLevel, double ratio)
{
  this->MaxMomentum = maxMomentum;
  this->LandauLevel = landauLevel;
  this->Ratio = ratio;
  this->HilbertSpaceDimension = this->MaxMomentum + 1;
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// x, y = coordinates where the function should be evaluated
// index = the function index 
// returns value of the wavefunction

Complex ParticleOnCylinderFunctionBasis::GetFunctionValue(double x, double y, double index)
{
  double L = sqrt(2.0 * M_PI * (this->MaxMomentum + 1.0) * this->Ratio);
  double kappa = 2.0 * M_PI/L;
  
  Complex Phase;
  Phase.Re = cos(kappa * index * y);
  Phase.Im = sin(kappa * index * y);

  double Normalization = 1.0/sqrt(L * sqrt(M_PI));
   
  Complex Result = Phase * exp(-0.5 * pow(x - kappa * index,2.0));

  if (this->LandauLevel == 1)
    {
      Result *= (2.0 * (x - kappa * index));    
      Normalization /= sqrt(2.0);
    }
  else if (this->LandauLevel > 1)
   {
     cout << "LL >= 2 " << endl;
     exit(1);
   }
 
  return (Result * Normalization);
}
