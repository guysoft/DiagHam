////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of Pfaffian wave function on sphere                 //
//                                                                            //
//                        last modification : 01/09/2004                      //
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
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"


// constructor
//
// nbrParticles = number of particles

PfaffianOnSphereWaveFunction::PfaffianOnSphereWaveFunction(int nbrParticles)
{
  this->NbrParticles = nbrParticles;
}

// copy constructor
//
// function = reference on the wave function to copy

PfaffianOnSphereWaveFunction::PfaffianOnSphereWaveFunction(const PfaffianOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
}

// destructor
//

PfaffianOnSphereWaveFunction::~PfaffianOnSphereWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PfaffianOnSphereWaveFunction::Clone ()
{
  return new PfaffianOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PfaffianOnSphereWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp;
  ComplexSkewSymmetricMatrix TmpPfaffian (this->NbrParticles);
  Complex WaveFunction(1.0);
  double Theta;
  double Phi;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      Theta = x[i << 1];
      Phi = x[1 + (i << 1)];
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = sin(0.5 * (x[j << 1] - Theta)) * cos(0.5 * (Phi - x[1 + (j << 1)]));
	  Tmp.Im = sin(0.5 * (Theta + x[j << 1])) * sin(0.5 * (Phi - x[1 + (j << 1)]));
	  WaveFunction *= Tmp;
	  Tmp = 1.0 / Tmp;
	  TmpPfaffian.SetMatrixElement (i , j, Tmp);
	}
    }
  WaveFunction *= TmpPfaffian.Pfaffian();
  return WaveFunction;
}

// Complex PfaffianOnSphereWaveFunction::operator ()(RealVector& x)
// {
//   Complex Tmp;
//   ComplexSkewSymmetricMatrix TmpPfaffian (this->NbrParticles);
//   Complex WaveFunction(1.0);
//   ComplexVector U(this->NbrParticles);
//   ComplexVector V(this->NbrParticles);
//   double c,s;
//   for (int i = 0; i < this->NbrParticles; ++i)
//     {
//       U[i].Re = cos(0.5 * x[i << 1]);
//       U[i].Im = U[i].Re;
//       U[i].Re *= (c=cos(0.5 * x[1 + (i << 1)]));
//       U[i].Im *= -(s=sin(0.5 * x[1 + (i << 1)]));
//       V[i].Re = sin(0.5 * x[i << 1]);
//       V[i].Im = V[i].Re;
//       V[i].Re *= c;
//       V[i].Im *= s;
//     }
//   for (int i = 0; i < this->NbrParticles; ++i)
//     {
//       for (int j = i + 1; j < this->NbrParticles; ++j)
// 	{
// 	  Tmp = U[i]*V[j]-U[j]*V[i];
// 	  WaveFunction *= Tmp*Tmp;
// 	  TmpPfaffian.SetMatrixElement (i , j, 1.0/Tmp);
// 	}
//     }
//   WaveFunction *= TmpPfaffian.Pfaffian();
//   return WaveFunction;
// }
