#include "config.h"
#include "Complex.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// evaluate the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
//
// coefficients = bi-dimensional array where the numerical coefficients have to be stored
// nbrStates = number of states in the Lowest Landau Level
void EvaluateNumericalCoefficients (double** coefficients, int nbrStates);

// evaluate the derivative of the Phi^4 contribution
//
// waveFunctionCoefficients = array containing the wave function coefficients
// numericalCoefficients = bi-dimensional array containing the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
// nbrStates = number of states in the Lowest Landau Level 
// return value = Phi^4 contribution
double EvaluatePhi4 (Complex* waveFunctionCoefficients, double** numericalCoefficients, int nbrStates);

// evaluate the derivative of the Phi^4 contribution with respect to one of the wave function coefficients
//
// waveFunctionCoefficients = array containing the wave function coefficients
// index = index of the contribution with respect to one of the wave function coefficients for which derivative has to be evaluated
// numericalCoefficients = bi-dimensional array containing the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
// nbrStates = number of states in the Lowest Landau Level 
// return value = derivative of the Phi^4 contribution
Complex EvaluatePhi4Derivative (Complex* waveFunctionCoefficients, int index, double** numericalCoefficients, int nbrStates);



int main(int argc, char** argv)
{
  cout.precision(14);
  int NbrStates = ;
  int MaxNbrIteration = ;
  double Precision = 1e-8;
  double** Coefficients = new double* [NbrStates];
  for (int i = 0; i < NbrStates; ++i)
    Coefficients[i] = new double [i + 1];
  EvaluateNumericalCoefficients(Coefficients, NbrStates);
  int NbrIteration = 0;
}


// evaluate the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
//
// coefficients = bi-dimensional array where the numerical coefficients have to be stored
// nbrStates = number of states in the Lowest Landau Level

void EvaluateNumericalCoefficients (double** coefficients, int nbrStates)
{
  int TwiceS = nbrStates - 1;
  double Factor = 2.0 / M_PI;
  for (int i = 0; i < nbrStates; ++i)
    for (int j = 0; j < i; ++j)
      {
	Coef.SetToOne();
	Coef.PartialFactorialMultiply(nbrStates - i , 2 * TwiceS - (i + j));
	Coef.PartialFactorialMultiply(nbrStates - j , 2 * TwiceS - (i + j));
	Coef.FactorialDivide(TwiceS + i);
	Coef.FactorialDivide(TwiceS + j);
	Coef.PartialFactorialDivide(nbrStates + 1, TwiceS + 1);	
	Coef.FactorialMultiply(nbrStates);
	coefficients[i][j] = sqrt(Coef.GetNumericalValue() * Factor);
      }
  for (int i = 0; i < nbrStates; ++i)
    {
      Coef.SetToOne();
      Coef.PartialFactorialMultiply(nbrStates - i , 2 * (TwiceS - i));
      Coef.FactorialDivide(TwiceS + i);
      coefficients[i][i] = 0.5 * Coef.GetNumericalValue();
      Coef.SetToOne();
      Coef.PartialFactorialDivide(nbrStates + 1, TwiceS + 1);	
      Coef.FactorialMultiply(nbrStates);
      coefficients[i][i] = sqrt(Coef.GetNumericalValue() * Factor);
    }
}

// evaluate the derivative of the Phi^4 contribution
//
// waveFunctionCoefficients = array containing the wave function coefficients
// numericalCoefficients = bi-dimensional array containing the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
// nbrStates = number of states in the Lowest Landau Level 
// return value = Phi^4 contribution

double EvaluatePhi4 (Complex* waveFunctionCoefficients, double** numericalCoefficients, int nbrStates)
{
  double Phi4 = 0.0;
  int i3;
  int i4;
  Complex Tmp;
  Complex Tmp2;
  for (int i1 = 0; i1 < nbrStates; ++i1)
    for (int i2 = 0; i2 <= i1; ++i2)
      {
	Tmp = (waveFunctionCoefficients[i1] * waveFunctionCoefficients[i2]) * numericalCoefficients[i1][i2];
	i3 = i1 + 2;
	if ((i3 & 1) == 0)
	  i3 >>= 1;
	else
	  {
	    i3 >>= 1;
	    i3 += 1;
	  }
	for (; i3 < nbrStates; ++i3)
	  {
	    i4 =  i1 + i2 - i3;
	    Tmp2 = Tmp;
	    Tmp2.ConjugateProduct(waveFunctionCoefficients[i3]);
	    Phi4 += (Tmp2.Re * waveFunctionCoefficients[i4].Re + Tmp2.Im * waveFunctionCoefficients[i4].Im) * numericalCoefficients[i3][i4];
	  }
      }
  return Phi4;    
}

// evaluate the derivative of the Phi^4 contribution with respect to one of the wave function coefficients
//
// waveFunctionCoefficients = array containing the wave function coefficients
// index = index of the contribution with respect to one of the wave function coefficients for which derivative has to be evaluated
// numericalCoefficients = bi-dimensional array containing the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
// nbrStates = number of states in the Lowest Landau Level 
// return value = derivative of the Phi^4 contribution

Complex EvaluatePhi4Derivative (Complex* waveFunctionCoefficients, int index, double** numericalCoefficients, int nbrStates)
{
  Complex DPhi4;
  int i3;
  int i4;
  Complex Tmp;
  Complex Tmp2;

  for (int i1 = 0; i1 < nbrStates; ++i1)
    {
      Tmp2 = (4.0 * numericalCoefficients[i1][index]) * waveFunctionCoefficients[i1];
      i3 = i1 + index;
      if ((i3 & 1) == 0)
	i3 >>= 1;
      else
	{
	  i3 >>= 1;
	  i3 += 1;
	}
      for (; i3 < nbrStates; ++i3)
	{
	  i4 =  i1 + index - i3;
	  Tmp = (waveFunctionCoefficients[i4] * waveFunctionCoefficients[i3]);
	  DPhi4 += (Tmp.Re * Tmp2) * numericalCoefficients[i3][i4];
	}      
    }
  return DPhi4;    
}

