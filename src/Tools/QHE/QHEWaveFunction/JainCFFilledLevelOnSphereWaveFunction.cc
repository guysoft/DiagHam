////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 16/09/2004                      //
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
#include "Tools/QHE/QHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
#include "Matrix/ComplexMatrix.h"
#include "MathTools/FactorialCoefficient.h"
#include "Vector/RealVector.h"


// constructor
//
// nbrParticles = number of particles
// nbrLandauLevel = number of Landau levels filled with composite fermions
// jastrowPower = power to which the Jastrow factor has to be raised

JainCFFilledLevelOnSphereWaveFunction::JainCFFilledLevelOnSphereWaveFunction(int nbrParticles, int nbrLandauLevels, int jastrowPower)
{
  this->NbrParticles = nbrParticles;
  this->NbrLandauLevels = nbrLandauLevels;
  this->JastrowPower = jastrowPower;
  this->Flag.Initialize();
  this->EvaluateNormalizationPrefactors();  
  this->EvaluateSumPrefactors();
}

// copy constructor
//
// function = reference on the wave function to copy

JainCFFilledLevelOnSphereWaveFunction::JainCFFilledLevelOnSphereWaveFunction(const JainCFFilledLevelOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrLandauLevels = function.NbrLandauLevels;
  this->JastrowPower = function.JastrowPower;
  this->Flag = function.Flag;
  this->NormalizationPrefactors = function.NormalizationPrefactors;
  this->SumPrefactors = function.SumPrefactors;
}

// destructor
//

JainCFFilledLevelOnSphereWaveFunction::~JainCFFilledLevelOnSphereWaveFunction()
{
  if (this->Flag.Shared() == false)
    {
      delete[] this->NormalizationPrefactors;
      for (int i = 0; i < this->NbrLandauLevels; ++i)
	{
	  for (int j = 0; j < this->NbrLandauLevels; ++j)
	    delete[] this->SumPrefactors[i][j];
	  delete[] this->SumPrefactors[i];
	}
      delete this->SumPrefactors;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* JainCFFilledLevelOnSphereWaveFunction::Clone ()
{
  return new JainCFFilledLevelOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex JainCFFilledLevelOnSphereWaveFunction::operator ()(RealVector& x)
{
  Complex* SpinorUCoordinates = new Complex[NbrParticles];
  Complex* SpinorVCoordinates = new Complex[NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
      SpinorUCoordinates[i].Im = SpinorUCoordinates[i].Re;
      SpinorUCoordinates[i].Re *= cos(0.5 * x[1 + (i << 1)]);
      SpinorUCoordinates[i].Im *= sin(0.5 * x[1 + (i << 1)]);
      SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
      SpinorVCoordinates[i].Im = SpinorVCoordinates[i].Re;
      SpinorVCoordinates[i].Re *= cos(0.5 * x[1 + (i << 1)]);
      SpinorVCoordinates[i].Im *= -sin(0.5 * x[1 + (i << 1)]);
    }

  Complex JastrowFactor(1.0);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  JastrowFactor *= ((SpinorUCoordinates[i] * SpinorVCoordinates[j]) - (SpinorUCoordinates[j] * SpinorVCoordinates[i]));
	}
    }
  Complex Tmp = JastrowFactor;
  for (int i = 1; i < this->JastrowPower; ++i)
    {
      JastrowFactor *= Tmp;
    }
  
  ComplexMatrix Slater (this->NbrParticles, this->NbrParticles);

  int MaxMomentum = (this->NbrParticles / this->NbrLandauLevels) - this->NbrLandauLevels;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      int Index = 0;
      int MaxMomentum2 = MaxMomentum;
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  for (int k = 0; k <= MaxMomentum2; ++k)
	    {
	      Slater.SetMatrixElement(Index, i, EvaluateCFMonopoleHarmonic(SpinorUCoordinates, SpinorVCoordinates, i, k, j, MaxMomentum2));
	      ++Index;
	    }
	  MaxMomentum2 += 2;
	}
    }  

  delete[] SpinorUCoordinates;
  delete[] SpinorVCoordinates;
  return JastrowFactor * Slater.Determinant();
}



// evaluate normalization factors of projected monopole harmonics
//

void JainCFFilledLevelOnSphereWaveFunction::EvaluateNormalizationPrefactors()
{
  this->NormalizationPrefactors = new double* [this->NbrLandauLevels];
  int TwiceQ = (this->NbrParticles / this->NbrLandauLevels) - this->NbrLandauLevels;
  int MaxMomentum = TwiceQ;
  FactorialCoefficient Coef;
  double Factor = ((double) (TwiceQ + 1)) / (4.0  * M_PI);
  this->NormalizationPrefactors[0] = new double[MaxMomentum + 1];
  for (int j = 0; j <= MaxMomentum; ++j)
    {
      Coef.SetToOne();
      Coef.PartialFactorialDivide(TwiceQ - j + 1, TwiceQ);
      Coef.FactorialMultiply(j);
      this->NormalizationPrefactors[0][j] = sqrt(Factor * Coef.GetNumericalValue());
      if ((j & 1) != 0)
	{
	  this->NormalizationPrefactors[0][j] *= -1.0;
	}
    }
  MaxMomentum += 2;
  for (int i = 1; i < this->NbrLandauLevels; ++i)  
    {
      double Factor = ((double) (TwiceQ + (2 * this->NbrLandauLevels) - 1)) / (4.0  * M_PI);
      this->NormalizationPrefactors[i] = new double[MaxMomentum + 1];
      for (int j = 0; j <= MaxMomentum; ++j)
	{	  
	  Coef.SetToOne();
	  Coef.FactorialMultiply(TwiceQ + 2 * (this->NbrLandauLevels - 1) - j);
	  Coef.FactorialDivide(this->NbrLandauLevels - 1);
	  Coef.PartialFactorialDivide(j + 1, TwiceQ + this->NbrLandauLevels - 1);
	  this->NormalizationPrefactors[i][j] = sqrt(Factor * Coef.GetNumericalValue());
	  if (((j + this->NbrLandauLevels - 1) & 1) != 0)
	    {
	      this->NormalizationPrefactors[i][j] *= -1.0;
	    }
	}
      MaxMomentum += 2;
    }
}

// evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
//

void JainCFFilledLevelOnSphereWaveFunction::EvaluateSumPrefactors()
{
  this->SumPrefactors = new double** [this->NbrLandauLevels];
  int TwiceQ = (this->NbrParticles / this->NbrLandauLevels) - this->NbrLandauLevels;
  int MaxMomentum = TwiceQ + 2 * (this->NbrLandauLevels - 1);
  FactorialCoefficient Coef;
  double Factor = 1.0;
  int ReducedNbrLandauLevels = this->NbrLandauLevels - 1;
  for (int i = 0; i < this->NbrLandauLevels; ++i)  
    {
      this->SumPrefactors[i] = new double* [i + 1];
      Factor = 1.0;
      for (int k = 0; k <= i; ++i)  
	{
	  this->SumPrefactors[i][k] = new double [MaxMomentum + 1];
	  for (int j = 0; j < (this->NbrLandauLevels - k); ++j)
	    this->SumPrefactors[i][k][j] = 0.0;
	  for (int j = this->NbrLandauLevels - k; j <= (MaxMomentum - k); ++j)
	    {
	      Coef.SetToOne();
	      Coef.PartialFactorialDivide(TwiceQ + 2 * (1 + this->JastrowPower * (this->NbrParticles - 1)), 
					   TwiceQ + 2 * (this->JastrowPower * (this->NbrParticles - 1)) + this->NbrLandauLevels);
	      Coef.PartialFactorialMultiply(k + 1, ReducedNbrLandauLevels);
	      Coef.FactorialDivide(k);
	      Coef.PartialFactorialMultiply(j + k - ReducedNbrLandauLevels, TwiceQ + ReducedNbrLandauLevels);	  
	      Coef.FactorialDivide(TwiceQ + (2 * ReducedNbrLandauLevels) -j - k);	  
	      this->SumPrefactors[i][k][j] = Factor * Coef.GetNumericalValue();
	    }
	  for (int j = MaxMomentum - k + 1; j <= MaxMomentum; ++j)
	    this->SumPrefactors[i][k][j] = 0.0;
	  Factor *= -1.0;
	}
    }
}


// evaluate composite fermion monopole spherical harmonic 
//
// spinorUCoordinates = spinor u coordinates where the function has to be evalauted
// spinorVCoordinates = spinor v coordinates where the function has to be evalauted
// coordinate = index of the main coordinate (aka coordinate before project onto the lowest Landau level)
// momentum = monopole spherical harmonic Lz momentum (plus S shift)
// landauLevel = index of the pseudo Landau level
// maximumMomentum = maxixum momentum that can be reached in the current pseudo Landau level
// return value = value of the monopole spherical harmonic at the givne point

Complex JainCFFilledLevelOnSphereWaveFunction::EvaluateCFMonopoleHarmonic (Complex* spinorUCoordinates, Complex* spinorVCoordinates,
									   int coordinate, int momentum, int landauLevel, int maximumMomentum)
{
  Complex Z (this->NormalizationPrefactors[landauLevel][momentum]);
  Complex Tmp = spinorUCoordinates[coordinate];
  if (momentum >= 0)
    {
      for (int i = 0; i < momentum; ++i)
	Z *= Tmp;
    }
  else
    {
      for (int i = momentum; i < 0; ++i)
	Z /= Tmp;
    }
  Tmp = spinorVCoordinates[coordinate];
  if (momentum <= (maximumMomentum - landauLevel))
    {
      for (int i = momentum; i < maximumMomentum; ++i)
	Z *= Tmp;  
    }
  else
    {
      for (int i = momentum; i < maximumMomentum; ++i)
	Z /= Tmp;
    }
  Tmp = 1.0;
  for (int i = 0; i <= landauLevel; ++i)
    {
      Tmp += this->SumPrefactors[landauLevel][i][maximumMomentum];// * this->EvaluateCFMonopoleHarmonicDerivative();
    }
  return (Z * Tmp);
}
