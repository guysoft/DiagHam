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

#include <iostream>

using std::cout;
using std::endl;


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
  this->JastrowPowerPowers = new double [this->NbrLandauLevels];
  this->JastrowPowerPowers[0] = 1.0;
  for (int i = 1; i < this->NbrLandauLevels; ++i)
    this->JastrowPowerPowers[i] = this->JastrowPowerPowers[i - 1] * ((double) this->JastrowPower);
  this->Flag.Initialize();
  this->EvaluateNormalizationPrefactors();  
  this->EvaluateSumPrefactors();
  this->JastrowFactorElements = new Complex*[this->NbrParticles];
  for (int i = 1; i < this->NbrParticles; ++i)
    this->JastrowFactorElements[i] = new Complex[i];
  this->DerivativeFactors = new Complex** [this->NbrParticles];
  this->DerivativeFactors2 = new Complex** [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->DerivativeFactors[i] = new Complex* [this->NbrLandauLevels];
      this->DerivativeFactors2[i] = new Complex* [this->NbrLandauLevels];
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  this->DerivativeFactors[i][j] = new Complex [this->NbrLandauLevels];
	  this->DerivativeFactors2[i][j] = new Complex [this->NbrLandauLevels];
	}
    }
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
  this->JastrowFactorElements = new Complex*[this->NbrParticles - 1];
  for (int i = 1; i < this->NbrParticles; ++i)
    this->JastrowFactorElements[i] = new Complex[i];
  this->DerivativeFactors = new Complex** [this->NbrParticles];
  this->DerivativeFactors2 = new Complex** [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->DerivativeFactors[i] = new Complex* [this->NbrLandauLevels];
      this->DerivativeFactors2[i] = new Complex* [this->NbrLandauLevels];
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  this->DerivativeFactors[i][j] = new Complex [this->NbrLandauLevels];
	  this->DerivativeFactors2[i][j] = new Complex [this->NbrLandauLevels];
	}
    }
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
      delete[] this->SumPrefactors;
      delete[] this->JastrowPowerPowers;
    }
  for (int i = 0; i < this->NbrParticles; ++i)
    delete[] this->JastrowFactorElements[i];
  delete[] this->JastrowFactorElements;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  delete[] this->DerivativeFactors[i][j];
	  delete[] this->DerivativeFactors2[i][j];
	}
      delete[] this->DerivativeFactors[i];
      delete[] this->DerivativeFactors2[i];
    }
  delete[] this->DerivativeFactors;
  delete[] this->DerivativeFactors2;
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
  Complex Tmp;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      for (int j = 0; j < i; ++j)
	{
	  Tmp = ((SpinorUCoordinates[i] * SpinorVCoordinates[j]) - (SpinorUCoordinates[j] * SpinorVCoordinates[i]));
	  this->JastrowFactorElements[i][j] = 1.0 / Tmp;
	  cout << this->JastrowFactorElements[i][j] << " ";
	  JastrowFactor *= Tmp;
	}
      cout << endl;
    }
  Tmp = JastrowFactor;
  for (int i = 1; i < this->JastrowPower; ++i)
    {
      JastrowFactor *= Tmp;
    }


  Complex Tmp2;
  for (int i = 0; i < this->NbrParticles; ++i)
    {     
      Complex** TmpDerivativeFactors = this->DerivativeFactors[i];
      Complex** TmpDerivativeFactors2 = this->DerivativeFactors2[i];
      for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
	for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
	  TmpDerivativeFactors[k1][k2] = 0.0;

      int Index = 0;
      for (int j = 1; j < this->NbrParticles; ++j)
	{
	  if (Index == i)
	    ++Index;
	  for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
	    {
	      Tmp = 1.0;
	      if (Index > i)
		Tmp2 = SpinorVCoordinates[Index] * SpinorUCoordinates[i] * this->JastrowFactorElements[Index][i];
	      else
		Tmp2 = -SpinorVCoordinates[Index] * SpinorUCoordinates[i] * this->JastrowFactorElements[i][Index];
	      for (int l = 0; l < k1; ++l)
		Tmp *= Tmp2;
	      for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
		{
		  TmpDerivativeFactors2[k1][k2] = Tmp;
		}
	    }
	  for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
	    {
	      Tmp = 1.0;
	      if (Index > i)
		Tmp2 = - SpinorUCoordinates[Index] * SpinorVCoordinates[i] * this->JastrowFactorElements[Index][i];
	      else
		Tmp2 = SpinorUCoordinates[Index] * SpinorVCoordinates[i] * this->JastrowFactorElements[i][Index];
	      for (int l = 0; l < k1; ++l)
		Tmp *= Tmp2;
	      for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
		{
		  TmpDerivativeFactors2[k2][k1] *= Tmp;
		}
	    }
	  for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
	    for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
	      TmpDerivativeFactors[k1][k2] += TmpDerivativeFactors2[k1][k2];
	  ++Index;
	}
      for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
	{
	  for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
	    cout << TmpDerivativeFactors[k1][k2] << " ";
	  cout << endl;
	}
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

  cout << Slater << endl;

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
	  Coef.FactorialMultiply(j);
	  Coef.FactorialDivide(TwiceQ + this->NbrLandauLevels - 1);
	  this->NormalizationPrefactors[i][j] = sqrt(Factor * Coef.GetNumericalValue());
	  if (((j + this->NbrLandauLevels - 1) & 1) != 0)
	    {
	      this->NormalizationPrefactors[i][j] *= -1.0;
	    }
	  cout <<  this->NormalizationPrefactors[i][j] << " ";
	}
      cout << endl;
      MaxMomentum += 2;
    }
  cout << endl;
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
      for (int k = 0; k <= i; ++k)  
	{
	  this->SumPrefactors[i][k] = new double [MaxMomentum + 1];
	  for (int j = 0; j < (i - k); ++j)
	    this->SumPrefactors[i][k][j] = 0.0;
	  for (int j = i - k; j <= (MaxMomentum - k); ++j)
	    {
	      Coef.SetToOne();
	      Coef.PartialFactorialDivide(TwiceQ + 2 * (1 + this->JastrowPower * (this->NbrParticles - 1)), 
					   TwiceQ + 2 * (this->JastrowPower * (this->NbrParticles - 1)) + i + 1);
	      Coef.PartialFactorialMultiply(k + 1, i);
	      Coef.FactorialDivide(i - k);
	      Coef.PartialFactorialMultiply(j + k - i, TwiceQ + i);	  
	      Coef.FactorialDivide(TwiceQ + (2 * i) -j - k);	  
	      this->SumPrefactors[i][k][j] = Factor * Coef.GetNumericalValue();
	      cout <<  i << " " << j << " " << k << " " << this->SumPrefactors[i][k][j] << endl;
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
// return value = value of the monopole spherical harmonic at the given point

Complex JainCFFilledLevelOnSphereWaveFunction::EvaluateCFMonopoleHarmonic (Complex* spinorUCoordinates, Complex* spinorVCoordinates,
									   int coordinate, int momentum, int landauLevel, int maximumMomentum)
{
  Complex Z (this->NormalizationPrefactors[landauLevel][momentum]);
  Complex Tmp = spinorUCoordinates[coordinate];
  cout << coordinate << " " << momentum << " " << landauLevel << endl;
  momentum -= landauLevel;
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
  cout << coordinate << " " << momentum << " " << landauLevel << endl;
  momentum = maximumMomentum - momentum - 2 * landauLevel;
  Tmp = spinorVCoordinates[coordinate];
  if (momentum > 0)
    {
      for (int i = 0; i < momentum; ++i)
	Z *= Tmp;  
    }
  else
    {
      for (int i = momentum; i < 0; ++i)
	Z /= Tmp;
    }
  cout << coordinate << " " << momentum << " " << landauLevel << " Z=" << Z << endl;
  Tmp = 0.0;
  for (int i = 0; i <= landauLevel; ++i)
    {
      Tmp += this->SumPrefactors[landauLevel][i][maximumMomentum] * this->EvaluateCFMonopoleHarmonicDerivative(coordinate, i, landauLevel - i);
      cout << this->SumPrefactors[landauLevel][i][maximumMomentum] << " " << this->EvaluateCFMonopoleHarmonicDerivative(coordinate, i, landauLevel - i) << endl;
    }
  return (Z * Tmp);
}


// evaluate derivative part of the composite fermion monopole spherical harmonic 
//
// index = particle index
// alpha = number of (d\du) derivates
// beta = number of (d\dv) derivates
// return value = derivative contribution to the monopole spherical harmonic

Complex JainCFFilledLevelOnSphereWaveFunction::EvaluateCFMonopoleHarmonicDerivative(int index, int alpha, int beta)
{
  Complex Tmp;
  switch (alpha)
    {
    case 0:
      {
	switch (beta)
	  {
	  case 0:
	    Tmp = 1.0;
	    break;
	  case 1:
	    Tmp = this->JastrowPowerPowers[1] * this->DerivativeFactors[index][0][1];
	    break;
	  case 2:
	    Tmp = (this->JastrowPowerPowers[2] * this->DerivativeFactors[index][0][1] * this->DerivativeFactors[index][0][1] 
		   - this->JastrowPowerPowers[1] * this->DerivativeFactors[index][0][2]);
	    break;
	  case 3:
	    Tmp = (this->JastrowPowerPowers[3] * this->DerivativeFactors[index][0][1] * this->DerivativeFactors[index][0][1] * this->DerivativeFactors[index][0][1]
		   - 3 * this->JastrowPowerPowers[2] * this->DerivativeFactors[index][0][2] * this->DerivativeFactors[index][0][1]
		   + 2 * this->JastrowPowerPowers[1] * this->DerivativeFactors[index][0][3]);
	    break;
	  }
      }
      break;
    case 1:
      {
	switch (beta)
	  {
	  case 0:
	    Tmp = this->JastrowPowerPowers[1] * this->DerivativeFactors[index][1][0];
	    break;
	  case 1:
	    Tmp = (this->JastrowPowerPowers[2] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][0][1] 
		   - this->JastrowPowerPowers[1] * this->DerivativeFactors[index][1][1]);
	    break;
	  case 2:
	    Tmp = (this->JastrowPowerPowers[3] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][0][1] * this->DerivativeFactors[index][0][1]
		   - this->JastrowPowerPowers[2] * (2 * this->DerivativeFactors[index][1][1] * this->DerivativeFactors[index][0][1]
						    + this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][0][2])
		   - 2 * this->JastrowPowerPowers[1] * this->DerivativeFactors[index][1][2]);
	    break;
	  }
      }
      break;
    case 2:
      {
	switch (beta)
	  {
	  case 0:
	    Tmp = (this->JastrowPowerPowers[2] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][1][0] 
		   - this->JastrowPowerPowers[1] * this->DerivativeFactors[index][2][0]);
	    break;
	  case 1:
	    Tmp = (this->JastrowPowerPowers[3] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][0][1]
		   - this->JastrowPowerPowers[2] * (2 * this->DerivativeFactors[index][1][1] * this->DerivativeFactors[index][1][0]
						    + this->DerivativeFactors[index][0][1] * this->DerivativeFactors[index][2][0])
		   - 2 * this->JastrowPowerPowers[1] * this->DerivativeFactors[index][2][1]);
	    break;
	  }
      }
      break;
    case 3:
      {
	switch (beta)
	  {
	  case 0:
	    Tmp = (this->JastrowPowerPowers[3] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][1][0]
		   - 3 * this->JastrowPowerPowers[2] * this->DerivativeFactors[index][2][0] * this->DerivativeFactors[index][1][0]
		   + 2 * this->JastrowPowerPowers[1] * this->DerivativeFactors[index][3][0]);
	    break;
	  }
      }
      break;	    
    }
  return Tmp;
}
