////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                             coulombian interaction                         //
//                                                                            //
//                        last modification : 18/07/2002                      //
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


#include "GeneralTools/StringTools.h"
#include "Hamiltonian/ParticleOnTorusPerturbedCoulombHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Polynomial/SpecialPolynomial.h"

#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define M1_12 0.08333333333333333


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// landauLevel = landauLevel to be simulated
// coulombPrefactor = prefactor multiplying Coulomb terms
// nbrPseudopotentials = number of perturbation in terms of additional pseudopotentials
// pseudopotentials = additional pseudopotential coefficients
// perturbationStrength = prefactor multiplying perturbation terms
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
ParticleOnTorusPerturbedCoulombHamiltonian::ParticleOnTorusPerturbedCoulombHamiltonian(ParticleOnTorus* particles, int nbrParticles, int maxMomentum, double ratio,
										       int landauLevel, double coulombPrefactor, int nbrPseudopotentials, double* pseudopotentials,
										       double perturbationStrength, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->LandauLevel = landauLevel;
  if (this->LandauLevel>=0)
    {
      this->FormFactor=LaguerrePolynomial(this->LandauLevel);
    }
  else
    {
      this->FormFactor=0.5*(LaguerrePolynomial(abs(this->LandauLevel))+LaguerrePolynomial(abs(this->LandauLevel)-1));
    }
  this->CoulombPrefactor=coulombPrefactor;
  cout << "FormFactor=" << this->FormFactor << endl;
  
  this->NbrPseudopotentials = nbrPseudopotentials;
  this->Pseudopotentials = new double[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->Pseudopotentials[i] = perturbationStrength*pseudopotentials[i];
  this->LaguerrePolynomials =new Polynomial[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  PrintMemorySize(cout,TmpMemory)<<endl;
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnTorusPerturbedCoulombHamiltonian::~ParticleOnTorusPerturbedCoulombHamiltonian() 
{
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusPerturbedCoulombHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Precision;
//  cout << "new coef====================================" << m1 << " "  << m2 << " "  << m3 << " "  << m4 << endl;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  // Coulomb terms
	  Coefficient = this->GetVofQ(PIOnM*Q2);
	  Precision = Coefficient;	  
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnM*Q2); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM*Q2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->NbrLzValue;
    }
  N2 = (double) (m1 - m4 - this->NbrLzValue);
  Coefficient = 1.0;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient =  this->GetVofQ(PIOnM*Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnM*Q2); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM*Q2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->NbrLzValue;
    }
  return (Sum / (2.0 * this->NbrLzValue));
}


// get Fourier transform of interaction
//
// Q2_half = one half of q^2 value
// return value = Fourier tranform

double ParticleOnTorusPerturbedCoulombHamiltonian::GetVofQ(double Q2_half)
{
  // Coulomb term
  double CoulombTerm = 1.0;
  if ((Q2_half != 0.0) && (this->CoulombPrefactor != 0.0))
    {
      CoulombTerm = this->FormFactor(Q2_half);
      CoulombTerm *= this->CoulombPrefactor*CoulombTerm;
      CoulombTerm /= sqrt(2.0*Q2_half);
    }
  else CoulombTerm=0.0;
  // additional pseudopotential terms
  double TmpInteraction = 0.0;
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    if (this->Pseudopotentials[i] != 0.0)
      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0*Q2_half);
  
  return exp(-2.0*Q2_half) * (TmpInteraction+CoulombTerm);
}
