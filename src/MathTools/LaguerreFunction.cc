#include "LaguerreFunction.h"
#include "MathTools/FactorialCoefficient.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::floor;



// default constructor
//
LaguerreFunction::LaguerreFunction()
{
  this->N = 0;
  this->Alpha = 0;
  this->Coefficients = new double;
  this->Coefficients[0] = 1.0;
}


// constructor for a Laguerre polynomial
//
// n = index
// alpha = modifier
LaguerreFunction::LaguerreFunction(int n, int alpha)
{
  this->N=n;
  this->Alpha=alpha;
  if (N>0)
    this->Coefficients=new double[N+1];
  else
    this->Coefficients=new double;
  FactorialCoefficient MyCoefficient; 
  for (int m=0; m<=N; ++m)
    {
      MyCoefficient.SetToOne();
      MyCoefficient.FactorialMultiply(N+Alpha);
      MyCoefficient.FactorialDivide(m);
      MyCoefficient.FactorialDivide(Alpha+m);
      MyCoefficient.FactorialDivide(N-m);
      this->Coefficients[m]=MyCoefficient.GetNumericalValue();
      if (m&1)
	this->Coefficients[m]*=-1.0;
    }
}


// copy constructor (duplicating datas)
//
// coefficients = reference on function to copy
LaguerreFunction::LaguerreFunction (const LaguerreFunction& laguerre)
{
  this->N=laguerre.N;
  this->Alpha=laguerre.Alpha;
  if (N>0)
    this->Coefficients=new double[N+1];
  else
    this->Coefficients=new double;
  for (int i=0; i<=N; ++i)
    this->Coefficients[i]=laguerre.Coefficients[i];
}

// destructor
//
LaguerreFunction::~LaguerreFunction ()
{
  if (this->N>0)
    delete [] this->Coefficients;
  else
    delete this->Coefficients;
}


// get the value of the function for a given coordinate x
//
// x = function argument
// return = function value at x

double LaguerreFunction::operator()(const double &x)
{
  double Res = this->Coefficients[this->N];
  for (int i = this->N - 1 ; i >= 0; i--)
    {
      Res *= x;
      Res += this->Coefficients[i];
    }
  return Res;
}


// get the value of the function for a given coordinate z
//
// x = function argument
// return = function value at x

double LaguerreFunction::GetValue(const double &x)
{
  double Res = this->Coefficients[this->N];
  for (int i = this->N - 1 ; i >= 0; i--)
    {
      Res *= x;
      Res += this->Coefficients[i];
    }
  return Res;
}

// get the value of the function for a given coordinate z
// values = vector where to store results
// manyX = vector of coordinates
void LaguerreFunction::GetManyValues(RealVector &values, RealVector &manyX)
{
  int Dim = manyX.GetVectorDimension();
  double *TmpValues = &(values[0]);
  double *TmpX = &(manyX[0]);
  for (int n=0; n<Dim; ++n)
    TmpValues[n] = this->Coefficients[this->N];
  for (int i = this->N - 1 ; i >= 0; i--)
    {
      for (int n=0; n<Dim; ++n)
	{
	  TmpValues[n] *= TmpX[n];
	  TmpValues[n] += this->Coefficients[i];
	}
    }
}

  // pretty-print a function value
  // str = stream to print to
  // z = point where to evaluate
ostream& LaguerreFunction::PrintValue(ostream &str, const double &x)
{
  str << "L_"<<this->N;
  if (this->Alpha!=0)
    str << "^"<<this->Alpha;
  str<<"("<<x<<")="<<this->GetValue(x);
  return str;
}
