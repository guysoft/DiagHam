////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of Laughlin wave function with one quasielectron on sphere      //
//                                                                            //
//                        last modification : 30/10/2008                      //
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
#include "Tools/FQHEWaveFunction/FQHESphereLaughlinOneQuasielectronWaveFunction.h"
#include "Vector/RealVector.h"
#include "GeneralTools/Endian.h"


#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;



// constructor
//
// nbrParticles = number of particles
// theta = position of the quasielectron (spherical coordinates, theta angle)
// phi = position of the quasielectron (spherical coordinates, phi angle)
// fermions = flag indicating whether to calculate bosonic or fermionic laughlin wave function

FQHESphereLaughlinOneQuasielectronWaveFunction::FQHESphereLaughlinOneQuasielectronWaveFunction(int nbrParticles, double theta, double phi, bool fermions)
{
  this->NbrParticles = nbrParticles;
  this->FermionFlag = fermions;

  this->UElectron.Re = cos(0.5*phi);
  this->UElectron.Im= -sin(0.5*phi);
  this->UElectron *= cos(0.5*theta);

  this->VElectron.Re = cos(0.5*phi);
  this->VElectron.Im = sin(0.5*phi);
  this->VElectron *= sin(0.5*theta);

  this->ConjUElectron = Conj(this->UElectron);
  this->ConjVElectron = Conj(this->VElectron);

  this->TmpJastrow = new Complex* [this->NbrParticles]; 
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpJastrow[i] = new Complex [this->NbrParticles]; 
  this->TmpWeights = new Complex [this->NbrParticles]; 
}

// copy constructor
//
// function = reference on the wave function to copy

FQHESphereLaughlinOneQuasielectronWaveFunction::FQHESphereLaughlinOneQuasielectronWaveFunction(const FQHESphereLaughlinOneQuasielectronWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;

  this->UElectron = function.UElectron;
  this->VElectron = function.VElectron;
  this->ConjUElectron = Conj(this->UElectron);
  this->ConjVElectron = Conj(this->VElectron);

  this->FermionFlag=function.FermionFlag;

  this->TmpJastrow = new Complex* [this->NbrParticles]; 
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpJastrow[i] = new Complex [this->NbrParticles]; 
  this->TmpWeights = new Complex [this->NbrParticles]; 
}

// destructor
//

FQHESphereLaughlinOneQuasielectronWaveFunction::~FQHESphereLaughlinOneQuasielectronWaveFunction()
{
  delete[] this->TmpWeights;
  for (int i = 0; i < this->NbrParticles; ++i)
    delete[] this->TmpJastrow[i]; 
  delete[] this->TmpJastrow;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHESphereLaughlinOneQuasielectronWaveFunction::Clone ()
{
  return new FQHESphereLaughlinOneQuasielectronWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex FQHESphereLaughlinOneQuasielectronWaveFunction::operator ()(RealVector& x)
{
  Complex WaveFunction(1.0);
  return WaveFunction;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)

Complex FQHESphereLaughlinOneQuasielectronWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  Complex Tmp;
  Complex TmpU1;
  Complex TmpV1;
  Complex WaveFunction(1.0);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      TmpU1 = uv[i << 1];
      TmpV1 = uv[1 + (i << 1)];
      this->TmpWeights[i] = ((this->ConjUElectron * TmpU1) + (this->ConjVElectron * TmpV1));
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp = (TmpU1 * uv[1 + (j << 1)]) - (uv[j << 1] * TmpV1);
	  this->TmpJastrow[i][j] = Tmp;
	  this->TmpJastrow[j][i] = -Tmp;
	  WaveFunction *= Tmp;
	}
    }

//   if (this->FermionFlag == true)
//     WaveFunction *= WaveFunction;
  
  if (this->FermionFlag == false)
    WaveFunction *= WaveFunction;
  else
    WaveFunction *= WaveFunction * WaveFunction;
  
  Complex WaveFunction2(0.0);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      for (int j = 0; j < i; ++j)
	WaveFunction2 += this->TmpWeights[j] / this->TmpJastrow[i][j];
      for (int j = i + 1; j < this->NbrParticles; ++j)
	WaveFunction2 += this->TmpWeights[j] / this->TmpJastrow[i][j];
   }


  //  cout << WaveFunction << " " << WaveFunction2 << endl;
  WaveFunction *= WaveFunction2;
  return WaveFunction;
}

