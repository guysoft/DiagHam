////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of Haldane wave function on sphere                 //
//                                                                            //
//                        last modification : 11/11/2006                      //
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
#include "Tools/FQHEWaveFunction/HalperinOnSphereWaveFunction.h"
#include "Vector/RealVector.h"


// constructor
//
// nbrSpinUpParticles = number of particles with spin up
// nbrSpinDownParticles = number of particles with spin down
// m1Index = m1 index
// m2Index = m2 index
// nIndex = n index

HalperinOnSphereWaveFunction::HalperinOnSphereWaveFunction(int nbrSpinUpParticles, int nbrSpinDownParticles, int m1Index, int m2Index, int nIndex)
{
  this->M1Index = m1Index;
  this->M2Index = m2Index;
  this->NIndex = nIndex;
  this->NbrSpinUpParticles = nbrSpinUpParticles;
  this->NbrSpinDownParticles = nbrSpinDownParticles;
  this->TotalNbrParticles = this->NbrSpinUpParticles + this->NbrSpinDownParticles;
}

// copy constructor
//
// function = reference on the wave function to copy

HalperinOnSphereWaveFunction::HalperinOnSphereWaveFunction(const HalperinOnSphereWaveFunction& function)
{
  this->M1Index = function.M1Index;
  this->M2Index = function.M2Index;
  this->NIndex = function.NIndex;
  this->NbrSpinUpParticles = function.NbrSpinUpParticles;
  this->NbrSpinDownParticles = function.NbrSpinDownParticles;
  this->TotalNbrParticles = function.TotalNbrParticles;
}

// destructor
//

HalperinOnSphereWaveFunction::~HalperinOnSphereWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* HalperinOnSphereWaveFunction::Clone ()
{
  return new HalperinOnSphereWaveFunction(*this);
}

// evaluate function at a given point (the first 2*nbrSpinUpParticles coordinates correspond to the position of the spin up particles, 
//                                     the other 2*nbrSpinDownParticles are spin down particle positions)
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex HalperinOnSphereWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp;
  double Theta;
  double Phi;
  double Factor = M_PI * 0.5;
  Complex TotalWaveFunction(1.0);
  if (this->M1Index > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrSpinUpParticles; ++i)
	{
	  Theta = x[i << 1];
	  Phi = x[1 + (i << 1)];
	  for (int j = i + 1; j < this->NbrSpinUpParticles; ++j)
	    {
	      Tmp.Re = Factor * sin(0.5 * (x[j << 1] - Theta)) * cos(0.5 * (Phi - x[1 + (j << 1)]));
	      Tmp.Im = Factor * sin(0.5 * (Theta + x[j << 1])) * sin(0.5 * (Phi - x[1 + (j << 1)]));
	      WaveFunction *= Tmp;
	    }
	}
      for (int i = 0; i < this->M1Index; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->M2Index > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = this->NbrSpinUpParticles; i < this->TotalNbrParticles; ++i)
	{
	  Theta = x[i << 1];
	  Phi = x[1 + (i << 1)];
	  for (int j = i + 1; j < this->TotalNbrParticles; ++j)
	    {
	      Tmp.Re = Factor * sin(0.5 * (x[j << 1] - Theta)) * cos(0.5 * (Phi - x[1 + (j << 1)]));
	      Tmp.Im = Factor * sin(0.5 * (Theta + x[j << 1])) * sin(0.5 * (Phi - x[1 + (j << 1)]));
	      WaveFunction *= Tmp;
	    }
	}
      for (int i = 0; i < this->M1Index; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->NIndex > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrSpinUpParticles; ++i)
	{
	  Theta = x[i << 1];
	  Phi = x[1 + (i << 1)];
	  for (int j = this->NbrSpinUpParticles; j < this->TotalNbrParticles; ++j)
	    {
	      Tmp.Re = Factor * sin(0.5 * (x[j << 1] - Theta)) * cos(0.5 * (Phi - x[1 + (j << 1)]));
	      Tmp.Im = Factor * sin(0.5 * (Theta + x[j << 1])) * sin(0.5 * (Phi - x[1 + (j << 1)]));
	      WaveFunction *= Tmp;
	    }
	}
      for (int i = 0; i < this->NIndex; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  return TotalWaveFunction;
}
