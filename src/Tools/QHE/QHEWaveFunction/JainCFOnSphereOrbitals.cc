////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2006 Nicolas Regnault and Gunnar Moeller       //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 07/05/2006                      //
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
#include "Tools/QHE/QHEWaveFunction/JainCFOnSphereOrbitals.h"
#include "Matrix/ComplexMatrix.h"
#include "MathTools/FactorialCoefficient.h"
#include "Vector/RealVector.h"

#include "DerivativeProductFactor.h"
#include "DerivativeProduct.h"
#include "SumDerivativeProduct.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//

JainCFOnSphereOrbitals::JainCFOnSphereOrbitals()
{
  this->NbrParticles = 0;
}

// constructor
//
// nbrParticles = number of particles
// nbrLandauLevel = number of Landau levels filled with composite fermions
// jastrowPower = power to which the Jastrow factor has to be raised

JainCFOnSphereOrbitals::JainCFOnSphereOrbitals(int nbrParticles, int nbrLandauLevels,
					       int nbrEffectiveFlux, int jastrowPower)
{
  this->NbrParticles = nbrParticles;
  this->NbrLandauLevels = nbrLandauLevels;
  this->TwiceS = abs(nbrEffectiveFlux);
  if (nbrEffectiveFlux<0) this->ReverseFluxFlag=true;
  else this->ReverseFluxFlag=false;
  this->NbrOrbitals = NbrLandauLevels*(NbrLandauLevels+TwiceS);
  this->Orbitals=new ComplexMatrix(NbrParticles,NbrOrbitals);
  this->ActualJastrowPower = jastrowPower;
  if ((this->ActualJastrowPower & 1) == 1)
    {
      this->JastrowPower = (jastrowPower + 1) >> 1;
    }
  else
    {
      this->JastrowPower = jastrowPower >> 1;
    }
  this->Flag.Initialize();
  this->EvaluateNormalizationPrefactors();
  this->EvaluateSumPrefactors();
  this->JastrowFactorElements = new Complex*[this->NbrParticles];
  for (int i = 1; i < this->NbrParticles; ++i)
    this->JastrowFactorElements[i] = new Complex[i];
  this->DerivativeFactors2 = new Complex** [this->NbrParticles];
  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
  this->SpinorUCoordinatePower = new Complex*[NbrParticles];
  this->SpinorVCoordinatePower = new Complex*[NbrParticles];
  if (ReverseFluxFlag)
    {
      this->MaxSpinorPower = (this->NbrLandauLevels - 1);
      this->MaxDerivativeNum = this->NbrLandauLevels + this->TwiceS;
      this->LLLDerivativeNum = this->TwiceS;
    }
  else
    {
      this->MaxSpinorPower = this->TwiceS + (this->NbrLandauLevels - 1)*2;
      this->MaxDerivativeNum = this->NbrLandauLevels;
      this->LLLDerivativeNum = 0;
    }
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinatePower[i] = new Complex[this->MaxSpinorPower + 1];
      this->SpinorVCoordinatePower[i] = new Complex[this->MaxSpinorPower + 1];      
      this->DerivativeFactors2[i] = new Complex* [this->MaxDerivativeNum];
      for (int j = 0; j < this->MaxDerivativeNum; ++j)
	{
	  this->DerivativeFactors2[i][j] = new Complex [this->MaxDerivativeNum-j];
	}
    }
  FactorialSignFactors = new double[MaxDerivativeNum];
  FactorialSignFactors[0]=1.0;
  FactorialSignFactors[1]=this->JastrowPower;
  for (int i=2; i<MaxDerivativeNum;++i)
    FactorialSignFactors[i] = (1-i)*FactorialSignFactors[i-1];
  MaxDerivativePower = new int*[MaxDerivativeNum];
  for (int j = 0; j < this->MaxDerivativeNum; ++j)
    MaxDerivativePower[j] = new int[MaxDerivativeNum-j];
  this->EvaluateDerivativeStructure();
}

// copy constructor
//
// function = reference on the wave function to copy

JainCFOnSphereOrbitals::JainCFOnSphereOrbitals(const JainCFOnSphereOrbitals& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrLandauLevels = function.NbrLandauLevels;
  this->NbrOrbitals =  function.NbrOrbitals;
  this->TwiceS = function.TwiceS;
  this->JastrowPower = function.JastrowPower;
  this->Flag = function.Flag;
  this->ReverseFluxFlag=function.ReverseFluxFlag;
  this->MaxSpinorPower = function.MaxSpinorPower;
  this->MaxDerivativeNum = function.MaxDerivativeNum;
  this->NormalizationPrefactors = function.NormalizationPrefactors;
  this->SumPrefactors = function.SumPrefactors;
  this->FactorialSignFactors = function.FactorialSignFactors;
  this->Orbitals=new ComplexMatrix(NbrParticles,NbrOrbitals);
  this->JastrowFactorElements = new Complex*[this->NbrParticles];
  for (int i = 1; i < this->NbrParticles; ++i)
    this->JastrowFactorElements[i] = new Complex[i];
  this->DerivativeFactors2 = new Complex** [this->NbrParticles];
  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
  this->SpinorUCoordinatePower = new Complex*[NbrParticles];
  this->SpinorVCoordinatePower = new Complex*[NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinatePower[i] = new Complex[MaxSpinorPower + 1];
      this->SpinorVCoordinatePower[i] = new Complex[MaxSpinorPower + 1];
      this->DerivativeFactors2[i] = new Complex* [this->MaxDerivativeNum];
      for (int j = 0; j < this->MaxDerivativeNum; ++j)
	this->DerivativeFactors2[i][j] = new Complex [this->MaxDerivativeNum-j];
    }
  for (int j = 0; j < this->MaxDerivativeNum; ++j)
    MaxDerivativePower[j] = new int[MaxDerivativeNum];
  this->EvaluateDerivativeStructure();
}

// destructor
//

JainCFOnSphereOrbitals::~JainCFOnSphereOrbitals()
{
  if (this->Flag.Shared() == false)
    {      
      for (int i = 0; i < this->NbrLandauLevels; ++i)
	{
	  delete[] this->NormalizationPrefactors[i];
	  for (int j = 0; j <= i; ++j)
	    delete[] this->SumPrefactors[i][j];
	  delete[] this->SumPrefactors[i];
	}
      delete[] this->NormalizationPrefactors;
      delete[] this->SumPrefactors;
      delete[] this->FactorialSignFactors;
    }
  for (int i = 1; i < this->NbrParticles; ++i)
    delete[] this->JastrowFactorElements[i];
  delete[] this->JastrowFactorElements;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      for (int j = 0; j < this->MaxDerivativeNum; ++j)
	delete[] this->DerivativeFactors2[i][j];
      delete[] this->DerivativeFactors2[i];
      delete[] this->SpinorUCoordinatePower[i];
      delete[] this->SpinorVCoordinatePower[i];
    }
  delete[] this->SpinorUCoordinates;
  delete[] this->SpinorVCoordinates;
  delete[] this->SpinorUCoordinatePower;
  delete[] this->SpinorVCoordinatePower;
  delete[] this->DerivativeFactors2;
  
  for (int k1=0; k1<MaxDerivativeNum;++k1)
    {
      for (int k2=0; k2<MaxDerivativeNum-k1;++k2)
	{
	  for (int pwr=0; pwr<MaxDerivativePower[k1][k2];++pwr)
	    delete[] DerivativeFactors[k1][k2][pwr];
	  delete[] DerivativeFactors[k1][k2];
	}
      delete[] DerivativeFactors[k1];
    }
  delete[] DerivativeFactors;
  delete Orbitals;
  for (int j = 0; j < this->MaxDerivativeNum; ++j)
    delete[] MaxDerivativePower[j];
  delete[] MaxDerivativePower;
  for (int i = 0; i < this->NbrLandauLevels; ++i)
    delete[] this->DerivativeStructure[i];
  delete[] this->DerivativeStructure;
}


// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = Slater-Determinant with CF Orbitals
ComplexMatrix& JainCFOnSphereOrbitals::operator ()(RealVector& x)
{
  int Index = 0;
  int MaxMomentum = this->TwiceS;
  this->EvaluateTables(x);
  for (int j = 0; j < this->NbrLandauLevels; ++j)
    {
      for (int k = 0; k <= MaxMomentum; ++k)
	{
	  this->EvaluateOrbitals(Index, k, j, MaxMomentum);
	  ++Index;
	}
      MaxMomentum += 2;
    }
  return *Orbitals;
}



// evaluate precalculation tables used during wave function evaluation (called at each evaluation)
//
// x = point where the function has to be evaluated
// derivativeFlag = indicate if precalculation tables invloved in derivative evaluation have to be calculated
// return value = value of the Jastrow factor

Complex JainCFOnSphereOrbitals::EvaluateTables(RealVector& x, bool derivativeFlag)
{
  double s,c;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
      this->SpinorUCoordinates[i].Im = this->SpinorUCoordinates[i].Re;
      this->SpinorUCoordinates[i].Re *= (c=cos(0.5 * x[1 + (i << 1)]));
      this->SpinorUCoordinates[i].Im *= (s=sin(0.5 * x[1 + (i << 1)]));
      this->SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
      this->SpinorVCoordinates[i].Im = this->SpinorVCoordinates[i].Re;
      this->SpinorVCoordinates[i].Re *= c;
      this->SpinorVCoordinates[i].Im *= -s;
    }

  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinatePower[i][0] = 1.0;
      this->SpinorVCoordinatePower[i][0] = 1.0;
      Complex TmpU = this->SpinorUCoordinates[i];
      Complex TmpV = this->SpinorVCoordinates[i];
      for (int j = 1; j <= this->MaxSpinorPower; ++j)
	{
	  this->SpinorUCoordinatePower[i][j] = this->SpinorUCoordinatePower[i][j - 1] * TmpU;
	  this->SpinorVCoordinatePower[i][j] = this->SpinorVCoordinatePower[i][j - 1] * TmpV;
	}
    }

  Complex JastrowFactor(1.0);
  Complex Tmp;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      for (int j = 0; j < i; ++j)
	{
	  Tmp = ((this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]));
	  this->JastrowFactorElements[i][j] = 1.0 / Tmp;
	  JastrowFactor *= Tmp;
	}
    }
  Tmp = JastrowFactor;
  for (int i = 1; i < this->ActualJastrowPower; ++i)
    {
      JastrowFactor *= Tmp;
    }

  if (derivativeFlag == true)
    {
      Complex Tmp2;
      for (int i = 0; i < this->NbrParticles; ++i)
	{     	  
	  Complex** TmpDerivativeFactors2 = this->DerivativeFactors2[i]; 
	  for (int k1 = 0; k1 < this->MaxDerivativeNum; ++k1)
	    for (int k2 = 0; k2 < this->MaxDerivativeNum-k1; ++k2)
	      DerivativeFactors[k1][k2][0][i] = 0.0;
	  
	  int Index = 0;
	  for (int j = 1; j < this->NbrParticles; ++j)
	    {
	      if (Index == i)
		++Index;
	      if (Index > i)
		Tmp2 = this->SpinorVCoordinates[Index] * this->JastrowFactorElements[Index][i];
	      else
		Tmp2 = -this->SpinorVCoordinates[Index] * this->JastrowFactorElements[i][Index];
	      Tmp = 1.0;
	      for (int k1 = 0; k1 < this->MaxDerivativeNum; ++k1)
		{
		  for (int k2 = 0; k2 < this->MaxDerivativeNum-k1; ++k2) // check!
		    {
		      TmpDerivativeFactors2[k1][k2] = Tmp;
		    }
		  Tmp *= Tmp2;
		}
	      if (Index > i)
		Tmp2 = - this->SpinorUCoordinates[Index] * this->JastrowFactorElements[Index][i];
	      else
		Tmp2 = this->SpinorUCoordinates[Index] * this->JastrowFactorElements[i][Index];
	      Tmp = 1.0;
	      for (int k1 = 0; k1 < this->MaxDerivativeNum; ++k1)
		{
		  for (int k2 = 0; k2 < this->MaxDerivativeNum-k1; ++k2) // check!
		    {
		      TmpDerivativeFactors2[k2][k1] *= Tmp;
		    }
		  Tmp *= Tmp2;
		}
	      for (int k1 = 0; k1 < this->MaxDerivativeNum; ++k1)
		for (int k2 = 0; k2 < this->MaxDerivativeNum-k1; ++k2) // check! 
		  DerivativeFactors[k1][k2][0][i] += TmpDerivativeFactors2[k1][k2];	      
	      ++Index;
	    }
	  for (int k1 = 0; k1 < this->MaxDerivativeNum; ++k1)
	    for (int k2 = 0; k2 < this->MaxDerivativeNum-k1; ++k2)
	      {
		DerivativeFactors[k1][k2][0][i] *= FactorialSignFactors[k1+k2]; // set proper prefactors...
		for (int pwr = 1; pwr < this->MaxDerivativePower[k1][k2]; ++pwr)
		  DerivativeFactors[k1][k2][pwr][i]=DerivativeFactors[k1][k2][0][i]*DerivativeFactors[k1][k2][pwr-1][i];
	      }
	}
    }
  return JastrowFactor;
}


// evaluate normalization factors of projected monopole harmonics
//

void JainCFOnSphereOrbitals::EvaluateNormalizationPrefactors()
{
  this->NormalizationPrefactors = new double* [this->NbrLandauLevels];
  int MaxMomentum = this->TwiceS;
  FactorialCoefficient Coef;
  double Factor = ((double) (this->TwiceS + 1)) / (4.0  * M_PI);
  this->NormalizationPrefactors[0] = new double[MaxMomentum + 1];
  if (ReverseFluxFlag)
    {  // normalization factors if q<0
      for (int j = 0; j <= MaxMomentum; ++j)
	{
	  Coef.SetToOne();
	  Coef.PartialFactorialDivide(this->TwiceS - j + 1, this->TwiceS);
	  Coef.FactorialMultiply(j);
	  this->NormalizationPrefactors[0][j] = sqrt(Factor * Coef.GetNumericalValue());
	}
      MaxMomentum += 2;
      for (int i = 1; i < this->NbrLandauLevels; ++i)  
	{
	  double Factor = ((double) (this->TwiceS + (2 * i) + 1)) / (4.0  * M_PI);
	  this->NormalizationPrefactors[i] = new double[MaxMomentum + 1];
	  double Sign = 1.0;
	  if (i & 1) Sign = -1.0; // sign (-1) ^ i
	  for (int j = 0; j <= MaxMomentum; ++j)
	    {	  
	      Coef.SetToOne();
	      Coef.FactorialMultiply(this->TwiceS + 2 * i - j);
	      Coef.FactorialDivide(i);
	      Coef.FactorialMultiply(j);
	      Coef.FactorialDivide(this->TwiceS + i);
	      this->NormalizationPrefactors[i][j] = Sign * sqrt(Factor * Coef.GetNumericalValue());
	    }
	  MaxMomentum += 2;
	}
    } // end q<0
  else
    { // normalization factors if q>=0
      for (int j = 0; j <= MaxMomentum; ++j)
	{
	  Coef.SetToOne();
	  Coef.PartialFactorialDivide(this->TwiceS - j + 1, this->TwiceS);
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
	  double Factor = ((double) (this->TwiceS + (2 * i) + 1)) / (4.0  * M_PI);
	  this->NormalizationPrefactors[i] = new double[MaxMomentum + 1];
	  double Sign = 1.0;
	  if ((this->TwiceS & 1) != 0)
	    Sign = -1.0;
	  for (int j = 0; j <= MaxMomentum; ++j)
	    {	  
	      Coef.SetToOne();
	      Coef.FactorialMultiply(this->TwiceS + 2 * i - j);
	      Coef.FactorialDivide(i);
	      Coef.FactorialMultiply(j);
	      Coef.FactorialDivide(this->TwiceS + i);
	      this->NormalizationPrefactors[i][j] = sqrt(Factor * Coef.GetNumericalValue());
	      this->NormalizationPrefactors[i][j] *= Sign;
	      Sign *= -1.0;
	    }
	  MaxMomentum += 2;
	}
    } // end q>0
}

// evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
//

void JainCFOnSphereOrbitals::EvaluateSumPrefactors()
{
  this->SumPrefactors = new double** [this->NbrLandauLevels];
  int MaxMomentum = this->TwiceS;
  FactorialCoefficient Coef;
  double Factor = 1.0;
  int TwiceBigQ;
  int TwiceQPrime=0;
  
  if (ReverseFluxFlag)
    {
      TwiceBigQ = 2*this->JastrowPower*( this->NbrParticles - 1 ) - TwiceS;
      TwiceQPrime = 2*this->JastrowPower*( this->NbrParticles - 1 );
    }
  else
    TwiceBigQ = 2*this->JastrowPower*( this->NbrParticles - 1 ) + TwiceS;
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
	      if (ReverseFluxFlag)
		{ // q<0
		  if (TwiceQPrime + i + 1 <= TwiceBigQ +1)
		    // = (TwiceBigQ + 1)! / (TwiceQPrime + i + 1)!
		    Coef.PartialFactorialMultiply( TwiceQPrime + i + 2 , TwiceBigQ +1 );  
		  else
		    // = (TwiceBigQ + 1)! / (TwiceQPrime + i + 1)!
		    Coef.PartialFactorialDivide( TwiceBigQ + 2, TwiceQPrime + i + 1);  
		}
	      else  // q>=0
		// this is (2Q+1)!/(2Q+n+1)!
		Coef.PartialFactorialDivide(this->TwiceS + 2 * (1 + this->JastrowPower * (this->NbrParticles - 1)), 
					    this->TwiceS + 2 * (this->JastrowPower * (this->NbrParticles - 1)) + i + 1);
	      Coef.PartialFactorialMultiply(k + 1, i);
	      Coef.FactorialDivide(i - k);
	      Coef.PartialFactorialMultiply(j + k - i + 1, this->TwiceS + i);	  
	      Coef.FactorialDivide(this->TwiceS + (2 * i) -j - k);	  
	      this->SumPrefactors[i][k][j] = Factor * Coef.GetNumericalValue();
	    }
	  for (int j = MaxMomentum - k + 1; j <= MaxMomentum; ++j)
	    this->SumPrefactors[i][k][j] = 0.0;
	  Factor *= -1.0;
	}
      MaxMomentum += 2;
    }
}

// recursive functions calculating the DerivativeStructure: Derivatives in V
//
SumDerivativeProduct JainCFOnSphereOrbitals::VRecursion(int VDerivatives)
{
  if (VDerivatives>1)
    {
      SumDerivativeProduct tmp=VRecursion(VDerivatives-1);
      SumDerivativeProduct rst=tmp;
      rst*=DerivativeProductFactor(this,0,1);
      rst+=tmp.Derivative(0,1);
      return (rst);
    }
  else if (VDerivatives==1)
    return SumDerivativeProduct(DerivativeProductFactor(this,0,1)); // = G
  else if (VDerivatives==0) return SumDerivativeProduct(DerivativeProductFactor(this,0,0)); // = 1
  else return SumDerivativeProduct(this); // =0
}

// recursive functions calculating the DerivativeStructure: Derivatives in U
//
SumDerivativeProduct JainCFOnSphereOrbitals::URecursion(int UDerivatives, int VDerivatives)
{
  if (UDerivatives>0)
    {
      SumDerivativeProduct tmp=URecursion(UDerivatives-1, VDerivatives);
      SumDerivativeProduct rst=tmp;
      rst*=DerivativeProductFactor(this,1,0);
      SumDerivativeProduct tmp2 = tmp.Derivative(1,0);
      rst+=tmp2;
      return (rst);
    }
  else if (UDerivatives==0) return VRecursion(VDerivatives);
  else return SumDerivativeProduct(this); // =0
}


// evaluate Structure of Derivatives
//
// DerivativeStructure[n][dU] initialized with J^-p (d/du)^dU (d/dv)^(LLLDerivativeNum+n-dU) J^p
//
// reserves proper memory for DerivativeFactor after counting the required powers...

void JainCFOnSphereOrbitals::EvaluateDerivativeStructure()
{
  int MaxDerivative=LLLDerivativeNum;
  
  DerivativeStructure = new SumDerivativeProduct*[this->NbrLandauLevels];
  for (int i=0; i<this->NbrLandauLevels; ++i)
    {
      DerivativeStructure[i] = new SumDerivativeProduct[MaxDerivative+1];
      cout << "For Landau Level " << i << endl;
      for (int j=0; j <= MaxDerivative; ++j)
	{
	  DerivativeStructure[i][j] = URecursion(j, MaxDerivative-j);
	  cout <<  "J^-p (d/du)^"<<j<<" (d/dv)^"<<MaxDerivative-j<<" J^p =" << DerivativeStructure[i][j] << endl;	  
	}
      ++MaxDerivative;
    }
  // find the highest Powers of Derivatives occurring in the sums:
  for (int i=0; i<MaxDerivativeNum;++i)
    for (int j=0; j<MaxDerivativeNum-i;++j)
      MaxDerivativePower[i][j]=0;
  for (int j=0; j < MaxDerivative; ++j)
    DerivativeStructure[this->NbrLandauLevels-1][j].TestHighestPowers();
  MaxDerivativePower[0][0]=1; // Add this by convention
  // Reserve DerivativeFactors as needed!
  DerivativeFactors = new Complex***[MaxDerivativeNum];
  for (int k1=0; k1<MaxDerivativeNum;++k1)
    {
      DerivativeFactors[k1] = new Complex**[MaxDerivativeNum-k1];
      for (int k2=0; k2<MaxDerivativeNum-k1;++k2)
	{
	  DerivativeFactors[k1][k2] = new Complex*[MaxDerivativePower[k1][k2]];
	  for (int pwr=0; pwr<MaxDerivativePower[k1][k2];++pwr)
	    DerivativeFactors[k1][k2][pwr] = new Complex[NbrParticles];
	}
    }
}



// evaluate composite fermion monopole spherical harmonic 
// Index = number of the orbital to be evaluated
// momentum = monopole spherical harmonic Lz momentum (plus S shift)
// landauLevel = index of the pseudo Landau level
// maximumMomentum = maxixum momentum that can be reached in the current pseudo Landau level

void JainCFOnSphereOrbitals::EvaluateOrbitals (int Index, int momentum, int landauLevel, int maximumMomentum)
{
  Complex *Tmp, *Tmp2;
  int s = landauLevel - momentum;
  if (s < 0)
    s = 0;
  int Max = maximumMomentum - momentum;
  int UPower, VPower, UDerivatives;
  if (this->ReverseFluxFlag)
    {
      UDerivatives=momentum + s - landauLevel;
      UPower = s;
      VPower = landauLevel - s;
    }
  else
    {
      UDerivatives = s;
      UPower = momentum + s - landauLevel;
      VPower = Max - s;
    }
  if (Max > landauLevel)
    Max = landauLevel;
  // get first term in the sum:
  Tmp = this->DerivativeStructure[landauLevel][UDerivatives].getValues();
  for (int i=0; i<NbrParticles; ++i) Tmp[i] *=
    this->SpinorUCoordinatePower[i][UPower] * this->SpinorVCoordinatePower[i][VPower] *
    this->SumPrefactors[landauLevel][s][momentum];
  ++s; ++UDerivatives; ++UPower; --VPower;
  for (; s <= Max; ++s)
    {
      Tmp2 = this->DerivativeStructure[landauLevel][UDerivatives].getValues();
      for (int i=0; i<NbrParticles; ++i)
	Tmp[i] += this->SpinorUCoordinatePower[i][UPower] * this->SpinorVCoordinatePower[i][VPower] *
	this->SumPrefactors[landauLevel][s][momentum]*Tmp2[i];
      ++UDerivatives; ++UPower; --VPower;
    }
  for (int i=0; i<NbrParticles; ++i)
    Orbitals->SetMatrixElement(i,Index,Tmp[i]*this->NormalizationPrefactors[landauLevel][momentum]);
}


