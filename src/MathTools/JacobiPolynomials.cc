#include "JacobiPolynomials.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::floor;


// maximum number of iterations before convergence has to be reached
#define ITERMAX 500


// default constructor
//
JacobiPolynomials::JacobiPolynomials()
{
  this->ParameterA = 0.0;
  this->ParameterB = 0.0;
  this->MaxDegreeN = 1;
  this->FunctionValues = new double[MaxDegreeN+1];
  this->FunctionValues[0]=1.0;
  this->CNPlusOne = NULL;
  this->CNConst = NULL;
  this->CNLinear = NULL;
  this->CNMinusOne = NULL;
  this->A1Zero = 0.5*(this->ParameterA-this->ParameterB);
  this->A1One = 0.5*(this->ParameterA+this->ParameterB)+1.0;
  this->LastN = -1;

}

// constructor for a general theta function \theta[^a_b](z|tau)
// maxDegreeN = maximum degree of the function
// param_a = parameter a
// param_b = parameter b
// tau = modulus
// precision = required precision
JacobiPolynomials::JacobiPolynomials(int maxDegreeN, double a, double b)
{
  this->ParameterA=a;
  this->ParameterB=b;
  this->MaxDegreeN=maxDegreeN;
  this->FunctionValues = new double[MaxDegreeN+1];
  
  this->CNPlusOne = new double[MaxDegreeN-1];
  this->CNConst = new double[MaxDegreeN-1];
  this->CNLinear = new double[MaxDegreeN-1];
  this->CNMinusOne = new double[MaxDegreeN-1];

  this->FunctionValues[0]=1.0;
  this->A1Zero = 0.5*(this->ParameterA-this->ParameterB);
  this->A1One = 0.5*(this->ParameterA+this->ParameterB)+1.0;

  this->InitializeRecursion();

  this->LastN = -1;  
}

// copy constructor (with duplicating of datas)
//
// 
JacobiPolynomials::JacobiPolynomials (const JacobiPolynomials& p)
{
  this->ParameterA=p.ParameterA;
  this->ParameterB=p.ParameterB;
  this->MaxDegreeN=p.MaxDegreeN;
  this->FunctionValues = new double[MaxDegreeN+1];
  for (int i=0; i<=MaxDegreeN; ++i)
    this->FunctionValues[i] = p.FunctionValues[i];
  if (MaxDegreeN>1)
    {
      this->CNPlusOne = new double[MaxDegreeN-1];
      this->CNConst = new double[MaxDegreeN-1];
      this->CNLinear = new double[MaxDegreeN-1];
      this->CNMinusOne = new double[MaxDegreeN-1];
      for (int i=0; i<MaxDegreeN-1; ++i)
	{
	  this->CNPlusOne[i] = p.CNPlusOne[i];
	  this->CNConst[i] = p.CNConst[i];
	  this->CNLinear[i] = p.CNLinear[i];
	  this->CNMinusOne[i] = p.CNMinusOne[i];
	}
    }

  this->A1Zero = 0.5*(this->ParameterA-this->ParameterB);
  this->A1One = 0.5*(this->ParameterA+this->ParameterB)+1.0;
  this->LastArgument = p.LastArgument;
  this->LastN = p.LastN;

}

// destructor
//
JacobiPolynomials::~JacobiPolynomials ()
{
  if (MaxDegreeN>1)
    {
      delete [] this->CNPlusOne;
      delete [] this->CNConst;
      delete [] this->CNLinear;
      delete [] this->CNMinusOne;
    }
  delete [] this->FunctionValues;
}


// get the value of the function for a given coordinate z
//
// n = degree (must be <=MaxDegreeN)
// x = argument
// return = function value of P_n(x)
double JacobiPolynomials::GetValue(int n, double x)
{
  if ((x!=LastArgument)||(n>LastN))
    {
      LastArgument=x;
      LastN=n;
      this->RunRecursion(n,x);
    }
  return FunctionValues[n];
}

// get the value of the function for a given coordinate z
//
// x = argument
// return = function value of P_n(x)
double* JacobiPolynomials::GetValues(double x)
{
  if ((x!=LastArgument)||(LastN<MaxDegreeN))
    {
      LastArgument=x;
      LastN=MaxDegreeN;
      this->RunRecursion(MaxDegreeN,x);
    }
  return this->FunctionValues;
}



  // pretty-print a function value
  // str = stream to print to
  // z = point where to evaluate
ostream& JacobiPolynomials::PrintValues(ostream &str, double x)
{
  this->GetValues(x);
  for (int i=0; i<=MaxDegreeN; ++i)
    str << "P_"<<i<<"^("<<this->ParameterA<<", "<<this->ParameterB<<") ("<<x<<")="<<this->FunctionValues[i]<<endl;
  return str;
}


// evaluate recursively up to degree n
void JacobiPolynomials::RunRecursion(int &limitN, double x)
{
  this->FunctionValues[0]=1.0;
  this->FunctionValues[1]=A1Zero+A1One*x;
  if (limitN>this->MaxDegreeN)
    {
      cout << "Attention: JacobiPolynomials object was initialized to calculate polynomials up to degree "<<this->MaxDegreeN<<endl;      
      limitN=this->MaxDegreeN;
    }
  for (int n=2; n<=limitN; ++n)
    this->FunctionValues[n] = ((CNConst[n-2] + CNLinear[n-2]*x)*this->FunctionValues[n-1] - CNMinusOne[n-2]*this->FunctionValues[n-2])/CNPlusOne[n-2];
}


// initialize recursion coefficients
void JacobiPolynomials::InitializeRecursion()
{
  FactorialCoefficient Pochammer;
  double SumAB = this->ParameterA+this->ParameterB;
  for (int n=1; n<this->MaxDegreeN; ++n)
    {
      CNPlusOne[n-1] = 2.0*(n+1.0)*(n+SumAB+1.0)*(2.0*n+SumAB);
      CNConst[n-1] = (2.0*n+SumAB+1.0)*(this->ParameterA*this->ParameterA-this->ParameterB*this->ParameterB);
      CNLinear[n-1] = (2.0*n+SumAB)*(2.0*n+SumAB+1.0)*(2.0*n+SumAB+2.0);
      CNMinusOne[n-1] = 2.0*(n+this->ParameterA)*(n+this->ParameterB)*(2.0*n+SumAB+2.0);
    }
}
