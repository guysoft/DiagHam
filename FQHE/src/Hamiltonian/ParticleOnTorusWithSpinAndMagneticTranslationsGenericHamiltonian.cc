////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//      two body pseudopotential interaction and magnetic translations        //
//                                                                            //
//                        last modification : 07/03/2014                      //
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


#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Architecture/AbstractArchitecture.h"

#include "Polynomial/SpecialPolynomial.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define M1_12 0.08333333333333333


// constructor from pseudopotentials
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// nbrPseudopotentialsUpUp = number of pseudopotentials for up-up interaction
// pseudopotentialsUpUp = pseudopotential coefficients for up-up interaction
// nbrPseudopotentialsDownDown = number of pseudopotentials for down-down interaction
// pseudopotentialsDownDown = pseudopotential coefficients for down-down interaction
// nbrPseudopotentialsUpDown = number of pseudopotentials for up-down interaction
// pseudopotentialsUpDown = pseudopotential coefficients for up-down interaction
// spinFluxUp = additional inserted flux for spin up
// spinFluxDown = additional inserted flux for spin down
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian::ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, 
																   int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
																   int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
																   int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
																   double spinFluxUp, double spinFluxDown, 
																   AbstractArchitecture* architecture, long memory, char* precalculationFileName, double* oneBodyPotentielUpUp, double* oneBodyPotentielDownDown, double* oneBodyPotentielUpDown)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;  
  this->InvRatio = 1.0 / ratio;
  double WignerEnergy = 0.0;
  this->SpinFluxUp = spinFluxUp;
  this->SpinFluxDown = spinFluxDown;
  this->Architecture = architecture;
  this->NbrPseudopotentialsUpUp = nbrPseudopotentialsUpUp;
  this->PseudopotentialsUpUp = new double[this->NbrPseudopotentialsUpUp];
  for (int i = 0; i < this->NbrPseudopotentialsUpUp; ++i)
    this->PseudopotentialsUpUp[i] = pseudopotentialsUpUp[i];
  this->NbrPseudopotentialsDownDown = nbrPseudopotentialsDownDown;
  this->PseudopotentialsDownDown = new double[this->NbrPseudopotentialsDownDown];
  for (int i = 0; i < this->NbrPseudopotentialsDownDown; ++i)
    this->PseudopotentialsDownDown[i] = pseudopotentialsDownDown[i];
  this->NbrPseudopotentialsUpDown = nbrPseudopotentialsUpDown;
  this->PseudopotentialsUpDown = new double[this->NbrPseudopotentialsUpDown];
  for (int i = 0; i < this->NbrPseudopotentialsUpDown; ++i)
    this->PseudopotentialsUpDown[i] = pseudopotentialsUpDown[i];

  this->MaxNbrPseudopotentials = this->NbrPseudopotentialsUpUp;
  if (this->NbrPseudopotentialsDownDown > this->MaxNbrPseudopotentials)
    this->MaxNbrPseudopotentials = this->NbrPseudopotentialsDownDown;
  if (this->NbrPseudopotentialsUpDown > this->MaxNbrPseudopotentials)
    this->MaxNbrPseudopotentials = this->NbrPseudopotentialsUpDown;
  this->LaguerrePolynomials =new Polynomial[this->MaxNbrPseudopotentials];
  for (int i = 0; i < this->MaxNbrPseudopotentials; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  if(oneBodyPotentielUpUp != 0)
    {
      this->OneBodyInteractionFactorsUpUp = new double[this->NbrLzValue];
      for(int i = 0; i < this->NbrLzValue; i++)
	this->OneBodyInteractionFactorsUpUp[i] = oneBodyPotentielUpUp[i];
    }
  this->OneBodyInteractionFactorsDownDown = 0;
  if(oneBodyPotentielDownDown != 0)
    {
      this->OneBodyInteractionFactorsDownDown = new double[this->NbrLzValue];
      for(int i = 0; i < this->NbrLzValue; i++)
	this->OneBodyInteractionFactorsDownDown[i] = oneBodyPotentielDownDown[i];
    } 
  this->OneBodyInteractionFactorsUpDown = 0;
  if(oneBodyPotentielUpDown != 0)
    {
      this->OneBodyInteractionFactorsUpDown = new double[this->NbrLzValue];
      for(int i = 0; i < this->NbrLzValue; i++)
	{
	  this->OneBodyInteractionFactorsUpDown[i] = oneBodyPotentielUpDown[i];
	  cout << this->OneBodyInteractionFactorsUpDown[i]<<" ";
	}
      cout <<endl;
    } 

  this->EvaluateInteractionFactors();
  this->EnergyShift = ((double) this->NbrParticles) * WignerEnergy;
  this->CosinusTable = new double [this->MaxMomentum];
  this->SinusTable = new double [this->MaxMomentum];
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      this->CosinusTable[i] = cos(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->MaxMomentum));
      this->SinusTable[i] = sin(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->MaxMomentum));
    }
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  int TmpMemory = this->FastMultiplicationMemory(memory);
	  if (TmpMemory < 1024)
	    cout  << "fast = " <<  TmpMemory << "b ";
	  else
	    if (TmpMemory < (1 << 20))
	      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	    else
	      if (TmpMemory < (1 << 30))
		cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	      else
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
	  cout << endl;
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

ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian::~ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian() 
{
  delete[] M12IntraValue;
  for (int i=0; i<NbrM12IntraIndices; ++i)
    delete[] M34IntraValues[i];
  delete[] M34IntraValues;
  delete[] NbrM34IntraValues;

  delete[] M12InterValue;
  for (int i=0; i<NbrM12InterIndices; ++i)
    delete[] M34InterValues[i];
  delete[] M34InterValues;
  delete[] NbrM34InterValues;

  delete[] this->CosinusTable;
  delete[] this->SinusTable;

  if (this->PseudopotentialsUpUp != 0)
    delete[] this->PseudopotentialsUpUp;
  if (this->PseudopotentialsDownDown != 0)
    delete[] this->PseudopotentialsDownDown;
  if (this->PseudopotentialsUpDown != 0)
    delete[] this->PseudopotentialsUpDown;
  if (this->LaguerrePolynomials != 0)
    delete[] this->LaguerrePolynomials;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] InteractionFactorsUpUp;
  delete[] InteractionFactorsDownDown;
  delete[] InteractionFactorsUpDown;

  delete[] OneBodyInteractionFactorsUpUp;
  delete[] OneBodyInteractionFactorsDownDown;  

  if (this->FastMultiplicationFlag == true)
    {
      for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
    }
  this->FastMultiplicationFlag = false;
  this->Particles = (ParticleOnTorusWithSpinAndMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian::ShiftHamiltonian (double shift)
{
}
  
// evaluate all interaction factors
//   

void ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian::EvaluateInteractionFactors()
{
  unsigned L16Mask = (1u<<16)-1;
  int Pos = 0;
  int M12Index = 0;
  int m4;
  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
  double MaxCoefficient = 0.0;  
  if (this->Particles->GetParticleStatistic() == ParticleOnTorusWithSpinAndMagneticTranslations::FermionicStatistic)
    {
      this->NbrM12IntraIndices = ((this->NbrLzValue-1) * (this->NbrLzValue - 2)) / 2;
      this->M12IntraValue = new unsigned [this->NbrM12IntraIndices];
      this->NbrM34IntraValues = new int [this->NbrM12IntraIndices];
      this->M34IntraValues = new unsigned* [this->NbrM12IntraIndices];
      for (int i = 0; i < this->NbrM12IntraIndices; ++i)
	this->M34IntraValues[i] = new unsigned[this->NbrLzValue];
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  {
	    for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	      {
		m4 = m1 + m2 - m3;
		if (m4 < 0)
		  m4 += this->MaxMomentum;
		else
		  if (m4 >= this->MaxMomentum)
		    m4 -= this->MaxMomentum;
		if (m3 > m4)
		  {
// 		    TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, 0.0)
// 					   + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, 0.0)
// 					   - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, 0.0)
// 					   - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, 0.0));
		    if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		      MaxCoefficient = fabs(TmpCoefficient[Pos]);
		    ++Pos;
		  }
	      }
	  }
      cout << "Max Nbr InteractionUpUp = " << Pos << endl;            
      MaxCoefficient *= MACHINE_PRECISION;
      InteractionFactorsUpUp = new double[Pos];
      InteractionFactorsDownDown = new double[Pos];
      M12Index = 0;
      Pos = 0;
      int TmpNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  {
	    this->M12IntraValue[M12Index] = (m1&L16Mask)|((m2&L16Mask)<<16);
	    this->NbrM34IntraValues[M12Index]=0;
	    for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	      {		
		m4 = m1 + m2 - m3;
		if (m4 < 0)
		  m4 += this->MaxMomentum;
		else
		  if (m4 >= this->MaxMomentum)
		    m4 -= this->MaxMomentum;
		if (m3 > m4)
		{
		  if  (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      this->InteractionFactorsUpUp[TmpNbrInteractionFactors] = TmpCoefficient[Pos];
		      this->InteractionFactorsDownDown[TmpNbrInteractionFactors] = TmpCoefficient[Pos];
		      this->M34IntraValues[M12Index][this->NbrM34IntraValues[M12Index]]
			= (m3&L16Mask)|((m4&L16Mask)<<16);
		      ++TmpNbrInteractionFactors;
		      ++this->NbrM34IntraValues[M12Index];
		    }
		  ++Pos;
		}
	    }
	    ++M12Index;
	  }
      cout << "Actual Nbr InteractionUpUp = " << TmpNbrInteractionFactors << endl;
      // matrix elements for different spin
      this->NbrM12InterIndices = (this->NbrLzValue-1) * (this->NbrLzValue-1);
      this->M12InterValue = new unsigned [this->NbrM12InterIndices];
      this->NbrM34InterValues = new int [this->NbrM12InterIndices];
      this->M34InterValues = new unsigned*[this->NbrM12InterIndices];      
      for (int i = 0; i < this->NbrM12InterIndices; ++i)
	this->M34InterValues[i] = new unsigned[this->NbrLzValue];
      Pos = 0;
      M12Index = 0;
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < this->MaxMomentum; ++m2)
	  {
	    for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	      {
		m4 = m1 + m2 - m3;
		if (m4 < 0)
		  m4 += this->MaxMomentum;
		else
		  if (m4 >= this->MaxMomentum)
		    m4 -= this->MaxMomentum;
		if ((m1 != m2) || (m3 != m4))
		  {
// 		    TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->LayerSeparation)
// 					   + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->LayerSeparation));
		  }
		else
		  {
// 		    TmpCoefficient[Pos] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->LayerSeparation);
		  }
		if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		  MaxCoefficient = fabs(TmpCoefficient[Pos]);
		++Pos;
	      }
	  }
      cout << "Max Nbr InteractionUpDown = " << Pos << endl;            
      MaxCoefficient *= MACHINE_PRECISION;
      InteractionFactorsUpDown = new double[Pos];
      M12Index = 0;
      Pos = 0;
      TmpNbrInteractionFactors = 0;
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < this->MaxMomentum; ++m2)
	  {
	    this->M12InterValue[M12Index] = (m1&L16Mask)|((m2&L16Mask)<<16);
	    this->NbrM34InterValues[M12Index]=0;
	    for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	      {		
		m4 = m1 + m2 - m3;
		if (m4 < 0)
		  m4 += this->MaxMomentum;
		else
		  if (m4 >= this->MaxMomentum)
		    m4 -= this->MaxMomentum;
		if  (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		  {
		    // swap 3,4, introduce additional minus sign
		    this->InteractionFactorsUpDown[TmpNbrInteractionFactors] = -1.0*TmpCoefficient[Pos];
		    this->M34InterValues[M12Index][this->NbrM34InterValues[M12Index]]
		      = (m4&L16Mask)|((m3&L16Mask)<<16); // rather than (m3&L16Mask)|((m4&L16Mask)<<16);
		    ++TmpNbrInteractionFactors;
		    ++this->NbrM34InterValues[M12Index];
		  }
		++Pos;
	      }
	    ++M12Index;
	  }
      cout << "Actual Nbr InteractionUpDown = " << TmpNbrInteractionFactors << endl;
    }
  else
    {
    }
  cout << "====================================" << endl;
  delete[] TmpCoefficient;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = pseudopotential coefficients
// spinFluxM1 = additional inserted flux for m1
// spinFluxM2 = additional inserted flux for m2
// spinFluxM3 = additional inserted flux for m3
// spinFluxM4 = additional inserted flux for m4
// return value = numerical coefficient

double ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
													double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double Factor =  - (((double) (m1-m3)) + spinFluxM1 - spinFluxM3) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = ((double) (m1 - m4)) + spinFluxM1 - spinFluxM4;
  double N1;
  double Q2;
  double Precision;
  double TmpInteraction;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0* PIOnM * Q2);
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
          if (fabs(Coefficient) != 0.0)
 	    Precision = Coefficient;
          else
            Precision = 1.0;
	}
       else
 	{
	  Precision = 1.0;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
	  Precision = 2.0 * exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->NbrLzValue;
    }
  N2 = (double) (m1 - m4 - this->NbrLzValue) + spinFluxM1 - spinFluxM4;
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  TmpInteraction = 0.0;
	  for (int i=0; i< nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
          if (fabs(Coefficient) != 0.0)
	    Precision = Coefficient;
          else
            Precision = 1.0;
	}
       else
 	{
	  Precision = 1.0;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
	  Precision = 2.0 *  exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->NbrLzValue;
    }
  //Normalize per flux (gives correct energy scale for 2-particle problem)
  return (Sum / ((double) this->NbrLzValue));
}
