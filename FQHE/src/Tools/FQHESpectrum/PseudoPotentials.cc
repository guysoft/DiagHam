#include "PseudoPotentials.h"

#include "Vector/RealVector.h"

#include "MathTools/Complex.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include <iostream>

using std::cout;
using std::endl;

// evalute pseudopotentials for coulomb interaction in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that conatins the pseudopotentials

double* EvaluatePseudopotentials(int nbrFlux, int landauLevel, double layerSeparation, bool quiet)
{
  cout.precision(14);
  int MaxMomentum = nbrFlux + (landauLevel << 1);
  double* Pseudopotentials = new double [MaxMomentum + 1];
  ClebschGordanCoefficients MainCoefficients(MaxMomentum, MaxMomentum);
  ClebschGordanCoefficients* Coefficients = new ClebschGordanCoefficients[MaxMomentum + 1];
  // new formfactors for finite thickness/separation:
  double *FormFactors = new double[MaxMomentum + 1];
  double dd = layerSeparation*layerSeparation;
  double base = ( sqrt(2.0*nbrFlux + dd) - layerSeparation ) / ( sqrt(2.0*nbrFlux + dd) + layerSeparation );
  FormFactors[0]=sqrt(base);
  for (int l = 1; l <= MaxMomentum; ++l)
    FormFactors[l] = FormFactors[l-1]*base;
  for (int l = 0; l <= MaxMomentum; ++l)
    Coefficients[l] = ClebschGordanCoefficients(MaxMomentum, l << 1);  
  for (int l = 0; l <= MaxMomentum; ++l)
    {
      double TmpPseudopotentials = 0.0;
      for (int m1 = -MaxMomentum; m1 <= MaxMomentum; m1 +=2)
	for (int m2 = -MaxMomentum; m2 <= MaxMomentum; m2 +=2)
	  {
	    double TmpCoef = 0.0;
	    double TmpCoef2;
	    int Min = abs(m1 - m2) >> 1;
	    double Sign = 1.0;
	    for (int j = Min; j <= MaxMomentum; ++j)
	      {
		TmpCoef2 = Coefficients[j].GetCoefficient(m1, m2 - m1, MaxMomentum) * Coefficients[j].GetCoefficient(nbrFlux, 0, MaxMomentum);
		TmpCoef += FormFactors[j] * Sign * TmpCoef2 * TmpCoef2;
		Sign *= -1.0;
	      }
	    TmpPseudopotentials += (MainCoefficients.GetCoefficient(m1, -m1, l << 1) * 
				    MainCoefficients.GetCoefficient(m2, -m2, l << 1) * TmpCoef);
	  }
      Pseudopotentials[MaxMomentum - l] = TmpPseudopotentials / sqrt (0.5 * nbrFlux);
      if (quiet == false)
	cout << "V[" << (MaxMomentum - l) << "] = " << Pseudopotentials[MaxMomentum - l] << endl;
    }
  delete[] Coefficients;
  delete[] FormFactors;
  return Pseudopotentials;
}

// evalute one body potentials for two impurities located at the poles in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// northPolePotential = potential of the impurity located at the north pole
// southPolePotential = potential of the impurity located at the south pole
// return value = array that conatins the pseudopotentials

double* EvaluateOneBodyPotentials(int nbrFlux, int landauLevel, double northPolePotential, double southPolePotential)
{
  int MaxMomentum = nbrFlux + (landauLevel << 1);
  double* OneBodyPotentials = new double [MaxMomentum + 1];
  for (int i = 0; i <= MaxMomentum; ++i)
   OneBodyPotentials[i] = 0.0;
  ParticleOnSphereGenericLLFunctionBasis Basis (nbrFlux, landauLevel);
  RealVector Value(2, true);
  Complex TmpValue;
  Value[0] = M_PI;
  Basis.GetFunctionValue(Value, TmpValue, landauLevel);
  OneBodyPotentials[landauLevel] = southPolePotential * SqrNorm(TmpValue);
  Value[0] = 0.0;
  Basis.GetFunctionValue(Value, TmpValue, MaxMomentum - landauLevel);
  OneBodyPotentials[MaxMomentum - landauLevel] = northPolePotential * SqrNorm(TmpValue);  
  return OneBodyPotentials;
}


// evaluate pseudopotentials coefficients of the monomials r^n in units of 1/R
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle of the corresponding LLL problem)
// exponentN = exponent of the monomial
// onlyOdd = boolean indidicating whether it's sufficient to reproduce only the odd pseudopotentials V_(2m+1)
// return value = array that conatins the coefficients V_m(r^n)
// where m runs over 0,...,nbrFlux, or if option onlyOdd given, from 0 to nbrFlux/2 with entries V_(2m+1)(r^n)
//
double* GetMonomialPseudopotentials(int nbrFlux, int exponentN, bool onlyOdd, bool verbose)
{
  FactorialCoefficient TmpCoeff;
  int sizeRst = (nbrFlux+1);
  if (onlyOdd) // fit only the odd pseudopotentials?
    {
      sizeRst /= 2;
      if (nbrFlux&1==0) --sizeRst;
    }    
  int spacing = (onlyOdd ?2:1);
  int M= (onlyOdd?1:0);
  double *rst = new double[sizeRst];  
  if ( exponentN%2 == 0)  // even exponents:
    {      
      TmpCoeff.Power2Multiply(exponentN);
      TmpCoeff.FactorialDivide(nbrFlux+exponentN/2 +1);
      TmpCoeff.FactorialMultiply(nbrFlux+ 1);
      TmpCoeff.FactorialDivide(nbrFlux+exponentN/2 +1);
      TmpCoeff.FactorialMultiply(nbrFlux+ 1);      
      double Prefactor =  TmpCoeff.GetNumericalValue();
      if (verbose) cout << "Pseudopotential Coefficients of <r^" << exponentN<<">"<<endl;
      for (int i=0; i<sizeRst; ++i, M+=spacing)
	{
	  int J=nbrFlux-M;	  
	  TmpCoeff.SetToOne();
	  TmpCoeff.FactorialMultiply(nbrFlux+exponentN/2-J);
	  TmpCoeff.FactorialDivide(nbrFlux-J);
	  TmpCoeff.FactorialMultiply(nbrFlux+exponentN/2+J+1);
	  TmpCoeff.FactorialDivide(nbrFlux+J+1);	  
	  rst[i]= Prefactor * TmpCoeff.GetNumericalValue() / sqrt (0.5 * nbrFlux);
	  if (verbose) cout << "V_"<<M<<"="<<rst[i]<<endl;
	}      
    }
  else // odd exponents:
    {
      int m=(exponentN+1)/2; // N = 2m-1
      if (verbose) cout << "Pseudopotential Coefficients of <r^" << exponentN<<">"<<endl;
      for (int i=0; i<sizeRst; ++i, M+=spacing)
	{
	  int J=nbrFlux-M;	  
	  TmpCoeff.SetToOne();
	  TmpCoeff.Power2Multiply(2*m+1);
	  
	  TmpCoeff.FactorialMultiply(2*nbrFlux+2*m+2*J+2);
	  TmpCoeff.FactorialDivide(nbrFlux+J+1);
	  TmpCoeff.FactorialDivide(nbrFlux+m+J+1);
	  
	  TmpCoeff.FactorialMultiply(2*nbrFlux+2*m-2*J);
	  TmpCoeff.FactorialDivide(nbrFlux-J);
	  TmpCoeff.FactorialDivide(nbrFlux+m-J);
	  
	  TmpCoeff.FactorialMultiply(nbrFlux+1);
	  TmpCoeff.FactorialMultiply(nbrFlux+m+1);
	  TmpCoeff.FactorialDivide(2*nbrFlux+2*m+2);

	  TmpCoeff.FactorialMultiply(nbrFlux+1);
	  TmpCoeff.FactorialMultiply(nbrFlux+m+1);
	  TmpCoeff.FactorialDivide(2*nbrFlux+2*m+2);
	  
	  rst[i]= TmpCoeff.GetNumericalValue() / sqrt (0.5 * nbrFlux);
	  if (verbose) cout << "V_"<<M<<"="<<rst[i]<<endl;
	}
    }
  return rst;
}

