////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of Pfaffian wave function with two quasiholes on sphere         //
//                                                                            //
//                        last modification : 18/07/2006                      //
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
#include "Tools/QHE/QHEWaveFunction/PfaffianOnSphereTwoQuasiholeWaveFunction.h"
#include "Vector/RealVector.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"


// constructor
//
// nbrParticles = number of particles
// theta1 = position of the first quasihole (spherical coordinates, theta angle)
// phi1 = position of the first quasihole (spherical coordinates, phi angle)
// theta2 = position of the second quasihole (spherical coordinates, theta angle)
// phi2 = position of the second quasihole (spherical coordinates, phi angle)

PfaffianOnSphereTwoQuasiholeWaveFunction::PfaffianOnSphereTwoQuasiholeWaveFunction(int nbrParticles, double theta1, double phi1, double theta2, double phi2)
{
  this->NbrParticles = nbrParticles;
}

// copy constructor
//
// function = reference on the wave function to copy

PfaffianOnSphereTwoQuasiholeWaveFunction::PfaffianOnSphereTwoQuasiholeWaveFunction(const PfaffianOnSphereTwoQuasiholeWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
}

// destructor
//

PfaffianOnSphereTwoQuasiholeWaveFunction::~PfaffianOnSphereTwoQuasiholeWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PfaffianOnSphereTwoQuasiholeWaveFunction::Clone ()
{
  return new PfaffianOnSphereTwoQuasiholeWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PfaffianOnSphereTwoQuasiholeWaveFunction::operator ()(RealVector& x)
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
