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
  this->EvaluateNormalizationPrefactors();  
}

// copy constructor
//
// function = reference on the wave function to copy

JainCFFilledLevelOnSphereWaveFunction::JainCFFilledLevelOnSphereWaveFunction(const JainCFFilledLevelOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrLandauLevels = function.NbrLandauLevels;
  this->JastrowPower = function.JastrowPower;
  this->NormalizationPrefactors = function.NormalizationPrefactors;
}

// destructor
//

JainCFFilledLevelOnSphereWaveFunction::~JainCFFilledLevelOnSphereWaveFunction()
{
  delete[] this->NormalizationPrefactors;
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
  
  switch (this->NbrLandauLevels)
    {
    case 1:
      {
	Complex TmpU;
	Complex TmpV;
	double* TmpPrefactors = this->NormalizationPrefactors[0];
	for (int i = 0; i < this->NbrParticles; ++i)
	  {
	    TmpU = SpinorUCoordinates[i];
	    TmpV = SpinorVCoordinates[i];
	    for (int j = i; j < this->NbrParticles; ++j)
	      {
		Tmp = 1.0;
		for (int k = 0; k < j; ++k)
		  Tmp *= TmpU;
		for (int k = this->NbrParticles - j - 1; k > 0; --k)
		  Tmp *= TmpV;
		Slater.SetMatrixElement(i, j, Tmp) * TmpPrefactors[j];
	      }
	    for (int j = 0; j < i; ++j)
	      {
		Tmp = 1.0;
		for (int k = 0; k < j; ++k)
		  Tmp *= TmpU;
		for (int k = this->NbrParticles - j - 1; k > 0; --k)
		  Tmp *= TmpV;
		Slater.SetMatrixElement(i, j, Tmp) * TmpPrefactors[j];
	      }
	  }
      }
      break;
    case 2:
      {
	Complex TmpU;
	Complex TmpV;
	int MaxLL0 = (this->NbrParticles >> 1) - 1;
	int MaxLL1 =  MaxLL0 + 2;
	int Max;	
	for (int i = 0; i < this->NbrParticles; ++i)
	  {
	    TmpU = SpinorUCoordinates[i];
	    TmpV = SpinorVCoordinates[i];
	    int Max = i + MaxLL0;
	    if (Max >= this->NbrParticles)
	      Max = this->NbrParticles
	    for (int j = i; j < Max; ++j)
	      {
		Tmp = 1.0;
		for (int k = 0; k < j; ++k)
		  Tmp *= TmpU;
		for (int k = this->NbrParticles - j - 1; k > 0; --k)
		  Tmp *= TmpV;
		Slater.SetMatrixElement(i, j, Tmp);
	      }
	    for (int j = 0; j < i; ++j)
	      {
		Tmp = 1.0;
		for (int k = 0; k < j; ++k)
		  Tmp *= TmpU;
		for (int k = this->NbrParticles - j - 1; k > 0; --k)
		  Tmp *= TmpV;
		Slater.SetMatrixElement(i, j, Tmp);
	      }
	  }
      }
      break;
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
  int NbrMomentum = TwiceQ;
  FactorialCoefficient Coef;
  double Factor = ((double) (TwiceQ + 1)) / (4.0  * M_PI);
  this->NormalizationPrefactors[0] = new double[MaxMomentum + 1];
  for (int j = 0; j <= MaxMomentum; ++j)
    {
      Coef.SetToOne();
      Coef.PartialFactorialMultiply(Q - j + 1, Q);
      Coef.FactorialDivide(j);
      this->NormalizationPrefactors[0][j] = sqrt(Factor * Coef.GetNumericalValue());
    }
  MaxMomentum += 2;
  for (int i = 1; i < this->NbrLandauLevels; ++i)  
    {
      double Factor = ((double) (TwiceQ + (2 * this->NbrLandauLevels) - 1)) / (4.0  * M_PI);
      this->NormalizationPrefactors[i] = new double[MaxMomentum + 1];
      for (int j = 0; j <= MaxMomentum; ++j)
	{
	  Coef.SetToOne();
	  Coef.FactorialMultiply(Q + 2 * (this->NbrLandauLevels - 1) - j);
	  Coef.FactorialDivide(this->NbrLandauLevels - 1);
	  Coef.PartialFactorialDivide(j + 1, Q + this->NbrLandauLevels - 1);
	  this->NormalizationPrefactors[i][j] = sqrt(Factor * Coef.GetNumericalValue());
	}
      MaxMomentum += 2;
    }
}

// evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
//

void JainCFFilledLevelOnSphereWaveFunction::EvaluateSumPrefactors()
{
  this->SumPrefactors = new double* [this->NbrLandauLevels];
  int TwiceQ = (this->NbrParticles / this->NbrLandauLevels) - this->NbrLandauLevels;
  int MaxMomentum = TwiceQ + 2 * (this->NbrLandauLevels - 1);
  FactorialCoefficient Coef;
  double Factor = 1.0;
  int ReducedNbrLandauLevels = this->NbrLandauLevels - 1;
  for (int i = 0; i < this->NbrLandauLevels; ++i)  
    {
      this->SumPrefactors[i] = new double [MaxMomentum + 1];
      for (int j = 0; j <= MaxMomentum; ++j)
	{
	  Coef.SetToOne();
	  Coef.PartialFactorialMultiply(j + 1, ReducedNbrLandauLevels);
	  Coef.FactorialDivide(j);
	  Coef.PartialFactorialMultiply(j + i, Q + ReducedNbrLandauLevels);	  
	  this->SumPrefactors[i][j] = sqrt(Factor * Coef.GetNumericalValue());
	}
      Factor *= -1.0;
    }

}
