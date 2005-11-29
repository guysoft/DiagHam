////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of function basis for particle on sphere             //
//                            in a given Landau level                         //
//                                                                            //
//                        last modification : 29/11/2005                      //
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
#include "FunctionBasis/QHEFunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <math.h>


// constructor
//
// lzMax = twice the maximum Lz value reached by a particle in the lowest Landau level (i.e. number of flux quanta)
// landauLevel = index of the Landau level to consider (0 for the lowest Landau level)

ParticleOnSphereGenericLLFunctionBasis::ParticleOnSphereGenericLLFunctionBasis(int lzMax, int landauLevel)
{
  this->LzMax = lzMax;
  this->LandauLevel = landauLevel;
  this->HilbertSpaceDimension = this->LzMax + (2 * this->LandauLevel) + 1;
  this->EvaluateSumPrefactors();
  this->EvaluateNormalizationPrefactors();
}

// destructor
//

ParticleOnSphereGenericLLFunctionBasis::~ParticleOnSphereGenericLLFunctionBasis ()
{
  delete[] this->NormalizationPrefactors;
  int MaxMomentum = this->LzMax + (2 * this->LandauLevel);      
  for (int i = 0; i <= MaxMomentum; ++i)
    delete[] this->SumPrefactors[i];
  delete[] this->SumPrefactors;
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void ParticleOnSphereGenericLLFunctionBasis::GetFunctionValue(RealVector& value, Complex& result, int index)
{
  Complex TmpU (cos(0.5 * value[1]), sin(0.5 * value[1]));
  Complex TmpV = Conj(TmpU);
  TmpU *= cos (0.5 * value[0]);
  TmpV *= sin (0.5 * value[0]);
  result.Re = 0.0;
  result.Im = 0.0;
  double* TmpCoefficients = this->SumPrefactors[index];
  for (int i = 0; i <= this->LandauLevel; ++i)
    {
    }  
  result *= this->NormalizationPrefactors[index];
}



// evaluate normalization factors of monopole harmonics
//

void ParticleOnSphereGenericLLFunctionBasis::EvaluateNormalizationPrefactors()
{
  int MaxMomentum = this->LzMax + (2 * this->LandauLevel);
  FactorialCoefficient Coef;
  this->NormalizationPrefactors = new double[MaxMomentum + 1];
  if (this->LandauLevel == 0)
    {
      double Factor = ((double) (MaxMomentum + 1)) / (4.0  * M_PI);
      for (int j = 0; j <= MaxMomentum; ++j)
	{
	  Coef.SetToOne();
	  Coef.PartialFactorialDivide(MaxMomentum - j + 1, MaxMomentum);
	  Coef.FactorialMultiply(j);
	  this->NormalizationPrefactors[j] = sqrt(Factor * Coef.GetNumericalValue());
	  if ((j & 1) != 0)
	    {
	      this->NormalizationPrefactors[j] *= -1.0;
	    }
	}
    }
  else
    {
      double Factor = ((double) (MaxMomentum + 1)) / (4.0  * M_PI);
      double Sign = 1.0;
      if ((this->LzMax & 1) != 0)
	Sign = -1.0;
      for (int j = 0; j <= MaxMomentum; ++j)
	{	  
	  Coef.SetToOne();
	  Coef.FactorialMultiply(MaxMomentum - j);
	  Coef.FactorialDivide(this->LandauLevel);
	  Coef.FactorialMultiply(j);
	  Coef.FactorialDivide(MaxMomentum - this->LandauLevel);
	  this->NormalizationPrefactors[j] = sqrt(Factor * Coef.GetNumericalValue());
	  this->NormalizationPrefactors[j] *= Sign;
	  Sign *= -1.0;
	}
    }
}

// evaluate constant factors that appears in the sum of projected monopole harmonic
//

void ParticleOnSphereGenericLLFunctionBasis::EvaluateSumPrefactors()
{
  if (this->LandauLevel == 0)
    {
      double TmpFactor = ((double) (this->LzMax + 1)) / (4.0 * M_PI);
      double TmpBinomial = 1.0;
      this->SumPrefactors = new double* [this->LzMax + 1];
      this->SumPrefactors[0] = new double [1];
      this->SumPrefactors[0][0] = sqrt (TmpBinomial * TmpFactor);
      for (int i = 1; i < this->HilbertSpaceDimension; ++i)
	{
	  this->SumPrefactors[i] = new double [1];
	  TmpBinomial *= this->LzMax - ((double) i) + 1.0;
	  TmpBinomial /= ((double) i);
	  this->SumPrefactors[i][0] = sqrt (TmpBinomial * TmpFactor);
	}
    }
  else
    {
      this->SumPrefactors = new double* [this->LandauLevel + 1];
      FactorialCoefficient Coef;
      int MaxMomentum = this->LzMax + (2 * this->LandauLevel);      
      double Factor = 1.0;
      if ((this->LandauLevel & 1) != 0)
	Factor = -1.0;
      for (int j = 0; j <= MaxMomentum; ++j)
	{
	  this->SumPrefactors[j] = new double [this->LandauLevel + 1];
	  double Factor = 1.0;
	  for (int k = 0; k <= this->LandauLevel; ++k)  
	    {
	      Coef.SetToOne();
	      Coef.PartialFactorialMultiply(k + 1, this->LandauLevel);
	      Coef.FactorialDivide(this->LandauLevel - k);
	      Coef.PartialFactorialMultiply(j + k - this->LandauLevel + 1, this->LzMax + this->LandauLevel);	  
	      Coef.FactorialDivide(this->LzMax + (2 * this->LandauLevel) -j - k);	  
	      this->SumPrefactors[k][j] = Factor * Coef.GetNumericalValue();
	      Factor *= -1.0;	      
	    }
	}
    }
}


