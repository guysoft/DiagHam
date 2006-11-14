////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of SU(4) generalized Halperine wave function on sphere       //
//                                                                            //
//                        last modification : 12/11/2006                      //
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
#include "Tools/FQHEWaveFunction/SU4HalperinOnSphereWaveFunction.h"
#include "Vector/RealVector.h"


// constructor
//
// nbrSpinUpIsospinPlusParticles = number of particles with spin up and isopsin plus
// nbrSpinUpIsospinMinusParticles = number of particles with spin up and isopsin minus
// nbrSpinDownIsospinPlusParticles = number of particles with spin down and isopsin plus
// nbrSpinDownIsospinMinusParticles = number of particles with spin down and isopsin minus
// mupIndex = power of the laughlin-like part for spin up - isopsin plus
// mumIndex = power of the laughlin-like part for spin up - isopsin minus
// mdpIndex = power of the laughlin-like part for spin down - isopsin plus
// mdmIndex = power of the laughlin-like part for spin down - isopsin minus
// nIntraIsospinIndex = power of the intra-isospin part (i.e (z_u{pm} -z_d{pm}))
// nIntraSpinIndex = power of the intra-spin part (i.e (z_{ud}p -z_{ud}m))
// nCrossSpinIsopinIndex = power of the cross spin-isospin part (i.e (z_up -z_dm) and (z_um -z_dp))

SU4HalperinOnSphereWaveFunction::SU4HalperinOnSphereWaveFunction(int nbrSpinUpIsospinPlusParticles, int nbrSpinUpIsospinMinusParticles, 
								 int nbrSpinDownIsospinPlusParticles, int nbrSpinDownIsospinMinusParticles,
								 int mupIndex, int mumIndex, int mdpIndex, int mdmIndex
								 int nIntraIsopinIndex, int nIntraSpinIndex, int nCrossSpinIsopinIndex)
{
  this->NbrSpinUpIsospinPlusParticles = nbrSpinUpIsospinPlusParticles;
  this->NbrSpinUpIsospinMinusParticles = nbrSpinUpIsospinMinusParticles;
  this->NbrSpinDownIsospinPlusParticles = nbrSpinDownIsospinPlusParticles;
  this->NbrSpinDownIsospinMinusParticles = nbrSpinDownIsospinMinusParticles;
  this->MupIndex = mupIndex;
  this->MumIndex = mumIndex;
  this->MdpIndex = mdpIndex;
  this->MdmIndex = mdmIndex;
  this->NIntraIsopinIndex = nIntraIsopinIndex;
  this->NIntraIsopinIndex = nIntraSpinIndex;
  this->NCrossSpinIsopinIndex = nCrossSpinIsopinIndex;
  this->TotalNbrParticles = (this->NbrSpinUpIsospinPlusParticles + this->NbrSpinUpIsospinMinusParticles + 
			     this->NbrSpinDownIsospinPlusParticles + this->NbrSpinDownIsospinMinusParticles);
}

// copy constructor
//
// function = reference on the wave function to copy

SU4HalperinOnSphereWaveFunction::SU4HalperinOnSphereWaveFunction(const SU4HalperinOnSphereWaveFunction& function)
{
  this->NbrSpinUpIsospinPlusParticles = function.NbrSpinUpIsospinPlusParticles;
  this->NbrSpinUpIsospinMinusParticles = function.NbrSpinUpIsospinMinusParticles;
  this->NbrSpinDownIsospinPlusParticles = function.NbrSpinDownIsospinPlusParticles;
  this->NbrSpinDownIsospinMinusParticles = function.NbrSpinDownIsospinMinusParticles;
  this->MupIndex = function.MupIndex;
  this->MumIndex = function.MumIndex;
  this->MdpIndex = function.MdpIndex;
  this->MdmIndex = function.MdmIndex;
  this->NIntraIsopinIndex = function.NIntraIsopinIndex;
  this->NIntraIsopinIndex = function.NIntraSpinIndex;
  this->NCrossSpinIsopinIndex = function.NCrossSpinIsopinIndex;
  this->TotalNbrParticles = function.TotalNbrParticles;
}

// destructor
//

SU4HalperinOnSphereWaveFunction::~SU4HalperinOnSphereWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* SU4HalperinOnSphereWaveFunction::Clone ()
{
  return new SU4HalperinOnSphereWaveFunction(*this);
}

// evaluate function at a given point (the first 2*nbrSpinUpIsopinPlusParticles coordinates correspond to the position of the spin up - isposin plus particles, 
//                                     the following 2*nbrSpinDownIsopinPlusParticles coordinates correspond to the position of spin down - isposin plus particles,
//                                     the other coordinates obey to the same scheme for isopin minus)
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex SU4HalperinOnSphereWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp;
  double Theta;
  double Phi;
  double Factor = M_PI * 0.5;
  Complex TotalWaveFunction(1.0);
  if (this->MupIndex > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrSpinUpIsospinPlusParticles; ++i)
	{
	  Theta = x[i << 1];
	  Phi = x[1 + (i << 1)];
	  for (int j = i + 1; j < this->NbrSpinUpIsospinPlusParticles; ++j)
	    {
	      Tmp.Re = Factor * sin(0.5 * (x[j << 1] - Theta)) * cos(0.5 * (Phi - x[1 + (j << 1)]));
	      Tmp.Im = Factor * sin(0.5 * (Theta + x[j << 1])) * sin(0.5 * (Phi - x[1 + (j << 1)]));
	      WaveFunction *= Tmp;
	    }
	}
      for (int i = 0; i < this->MupIndex; ++i)
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
