#include "Options/Options.h"
#include "MathTools/BinomialCoefficients.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#ifdef __GSL__
#include <gsl/gsl_sf_hyperg.h>
#endif

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// compute the one-body matrix element for a step confining potential
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutPosition = x position of the step
// return value = one-body matrix element
double FQHECylinderComputeStepPotentialCoefficients (double orbitalIndex, double perimeter, double cutPosition);

// compute the one-body matrix element for a polynomial right confining potential
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// alpha = exponent of the polynomial
// cutPosition = x position below which the potential should be set to zero 
// binomials = reference to the binomial coefficients
// return value = one-body matrix element
double FQHECylinderComputePolynomialRightPotentialCoefficients (double orbitalIndex, double perimeter, int alpha, double cutPosition, BinomialCoefficients& binomials);

// compute the one-body matrix element for a polynomial left confining potential
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// alpha = exponent of the polynomial
// cutPosition = x position above which the potential should be set to zero 
// binomials = reference to the binomial coefficients
// return value = one-body matrix element
double FQHECylinderComputePolynomialLeftPotentialCoefficients (double orbitalIndex, double perimeter, int alpha, double cutPosition, BinomialCoefficients& binomials);

// compute the one-body matrix element for a polynomial right confining potential defined in momentum space
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// alpha = exponent of the polynomial
// cutPosition = x position below which the potential should be set to zero 
// return value = one-body matrix element
double FQHECylinderComputePolynomialRightMomentumPotentialCoefficients (double orbitalIndex, double perimeter, int alpha, double cutPosition);

// compute the one-body matrix element for a polynomial left confining potential
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// alpha = exponent of the polynomial
// cutPosition = x position above which the potential should be set to zero 
// return value = one-body matrix element
double FQHECylinderComputePolynomialLeftMomentumPotentialCoefficients (double orbitalIndex, double perimeter, int alpha, double cutPosition);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHECylinderConfiningPotentialCoefficients" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (if zero, assume an infinite cylinder)", 0);
  (*SystemGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the cylinder", 1);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "error", "error below which a matrix element is considered equal to zero", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "confining-rightpower", "integer exponent of the right confining potential", 2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "confining-rightstrength", "strength of the right confining potential", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "confining-rightoffset", "position (with respect to half of the cylinder) below which the right confining potential is equal to zero", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "confining-leftpower", "integer exponent of the left confining potential", 2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "confining-leftstrength", "strength of the left confining potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "confining-leftoffset", "position (with respect to half of the cylinder) above which the left confining potential is equal to zero", 0.0);
  (*SystemGroup) += new  BooleanOption ('\n', "confining-momentum", "if set, the confining potential is defined in the momentum space");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "optional output file name instead of the default one");
  (*OutputGroup) += new BooleanOption ('\n', "spinful", "the output file will be used for a spinful system, assuming the same confining potential for all species");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderConfiningPotentialCoefficients -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");
  double Error = Manager.GetDouble("error");
  if ((NbrFluxQuanta == 0) && (Error == 0.0))
    {
      Error = MACHINE_PRECISION;
    }

  double Perimeter = Manager.GetDouble("cylinder-perimeter");
  if (Perimeter == 0.0)
    {
      if (NbrFluxQuanta == 0)
	{
	  cout << "error, nbr-flux has to be provided if cylinder-perimeter is zero" << endl;
	  return 0;
	}
      else
	{
	  Perimeter = sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) * Manager.GetDouble("aspect-ratio"));
	}
    }

  int RightAlpha = Manager.GetInteger("confining-rightpower");
  double RightV0 = Manager.GetDouble("confining-rightstrength");
  double RightShift = Manager.GetDouble("confining-rightoffset");
  int LeftAlpha = Manager.GetInteger("confining-leftpower");
  double LeftV0 = Manager.GetDouble("confining-leftstrength");
  double LeftShift = Manager.GetDouble("confining-leftoffset");

  int MaxBinomial = RightAlpha;
  if (LeftAlpha > MaxBinomial)
    MaxBinomial = LeftAlpha;
  BinomialCoefficients TmpBinomialCoefficients(MaxBinomial);

  double CutPosition = 0.0;
  char* OutputFile = 0;
  if (Manager.GetString("output-file") == 0)
    {
      OutputFile = new char[512];
      if (Manager.GetBoolean("confining-momentum") == true)
	{
	  sprintf (OutputFile, "confining_momentum_cylinder_perimeter_%.6f_2s_%d_alphar_%d_x0r_%.6f_v0r_%.6f_alphal_%d_x0l_%.6f_v0l_%.6f.dat", 
		   Perimeter, NbrFluxQuanta, RightAlpha, RightShift, RightV0, LeftAlpha, LeftShift, LeftV0);
	}
      else
	{
	  sprintf (OutputFile, "confining_cylinder_perimeter_%.6f_2s_%d_alphar_%d_x0r_%.6f_v0r_%.6f_alphal_%d_x0l_%.6f_v0l_%.6f.dat", 
		   Perimeter, NbrFluxQuanta, RightAlpha, RightShift, RightV0, LeftAlpha, LeftShift, LeftV0);
	}
    }
  else
    {
      OutputFile = new char[strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFile, Manager.GetString("output-file"));
    }
  ofstream File;
  File.open(OutputFile, ios::binary | ios::out);
  File.precision(14);
  if (Manager.GetBoolean("confining-momentum") == true)
    {
      File << "# confining potential in momentum space defined by :" << endl;
    }
  else
    {
      File << "# confining potential defined by :" << endl;
    }
  File << "# right alpha = " << RightAlpha << ", right V0 = " << RightV0 << ", right shift = " << RightShift << endl;
  File << "# left alpha = " << LeftAlpha << ", left V0 = " << LeftV0 << ", left shift = " << LeftShift << endl;
  if (NbrFluxQuanta == 0)
    {
      File << "# on an infinite cylinder with perimeter L=" << Perimeter;
    }
  else
    {
      File << "# on a cylinder with perimeter L=" << Perimeter << " and N_phi=" << NbrFluxQuanta;
    }
  File << endl;
  int NbrCoefficients = 0;
  double* Coefficients = 0;
  if (NbrFluxQuanta == 0)
    {
      int MaxNbrCoefficients = 1000; 
      double* TmpCoefficients = new double [MaxNbrCoefficients];
      for (int i = 0; i < MaxNbrCoefficients; ++i)
	{
	  TmpCoefficients[i] = FQHECylinderComputeStepPotentialCoefficients(i, Perimeter, CutPosition);
	  if (TmpCoefficients[i] < Error)
	    {
	      i = MaxNbrCoefficients;
	    }
	  else
	    {
	      ++NbrCoefficients;
	    }	    
	}
      Coefficients = new double [2 * NbrCoefficients - 1];
      Coefficients[NbrCoefficients - 1] =  TmpCoefficients[0];
      for (int i = 1; i < NbrCoefficients; ++i)
	{
	  Coefficients[NbrCoefficients - 1 + i] =  TmpCoefficients[i];
	  Coefficients[NbrCoefficients - 1 - i] =  1.0 - TmpCoefficients[i];
	}
      delete[] TmpCoefficients;
    }
  else
    {
      Coefficients = new double [NbrFluxQuanta + 1];
      for (NbrCoefficients = 0; NbrCoefficients <= NbrFluxQuanta; ++NbrCoefficients)
	{
	  Coefficients[NbrCoefficients] = 0.0;
	}
      if (RightV0 != 0.0)
	{
	  for (NbrCoefficients = 0; NbrCoefficients <= NbrFluxQuanta; ++NbrCoefficients)
	    {
	      double TmpCoefficient;
	      if (Manager.GetBoolean("confining-momentum") == true)
		{
		  TmpCoefficient = RightV0 * FQHECylinderComputePolynomialRightMomentumPotentialCoefficients(((double) NbrCoefficients) - 0.5 * ((double) NbrFluxQuanta), 
													     Perimeter, RightAlpha, RightShift);
		}
	      else
		{
		  if (RightAlpha == 0)
		    {
		      TmpCoefficient = RightV0 * (1.0 - FQHECylinderComputeStepPotentialCoefficients(((double) NbrCoefficients) - 0.5 * ((double) NbrFluxQuanta), Perimeter, RightShift));
		    }
		  else
		    {
		      TmpCoefficient = RightV0 * FQHECylinderComputePolynomialRightPotentialCoefficients(((double) NbrCoefficients) - 0.5 * ((double) NbrFluxQuanta), 
													 Perimeter, RightAlpha, RightShift, TmpBinomialCoefficients);
		    }
		}
	      if (TmpCoefficient > Error)
		{
		  Coefficients[NbrCoefficients] += TmpCoefficient;
		}
	    }
	}
      if (LeftV0 != 0.0)
	{
	  for (NbrCoefficients = 0; NbrCoefficients <= NbrFluxQuanta; ++NbrCoefficients)
	    {
	      double TmpCoefficient;
	      if (Manager.GetBoolean("confining-momentum") == true)
		{
		  TmpCoefficient = LeftV0 * FQHECylinderComputePolynomialLeftMomentumPotentialCoefficients(((double) NbrCoefficients) - 0.5 * ((double) NbrFluxQuanta), 
													     Perimeter, LeftAlpha, LeftShift);
		}
	      else
		{
		  // 	      if (LeftAlpha == 0)
		  // 		{
		  // 		  TmpCoefficient = LeftV0 * FQHECylinderComputeStepPotentialCoefficients(((double) NbrCoefficients) - 0.5 * ((double) NbrFluxQuanta), Perimeter, LeftShift);
		  // 		}
		  // 	      else
		  {
		    TmpCoefficient = LeftV0 * FQHECylinderComputePolynomialLeftPotentialCoefficients(((double) NbrCoefficients) - 0.5 * ((double) NbrFluxQuanta), 
												     Perimeter, LeftAlpha, LeftShift, TmpBinomialCoefficients);
		  }
		}
	      // 	      if (TmpCoefficient > Error)
	      {
		  Coefficients[NbrCoefficients] += TmpCoefficient;
		}
	    }
	}
    }
  if (Manager.GetBoolean("spinful") == false)
    {
      File << "OneBodyPotentials =";
      for (int i = 0; i < NbrCoefficients; ++i)
	File << " " << Coefficients[i];
      File << endl;
    }
  else
    {
      File << "OneBodyPotentialUpUp     =";
      for (int i = 0; i < NbrCoefficients; ++i)
	File << " " << Coefficients[i];
      File << endl;
      File << "OneBodyPotentialDownDown =";
      for (int i = 0; i < NbrCoefficients; ++i)
	File << " " << Coefficients[i];
      File << endl;
    }
  File.close();
  delete[] Coefficients;
  return 0;
}


// compute the one-body matrix element for a step confining potential
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutPosition = x position of the step
// return value = one-body matrix element

double FQHECylinderComputeStepPotentialCoefficients (double orbitalIndex, double perimeter, double cutPosition)
{  
  return (0.5 * (1.0 + erf(cutPosition - (orbitalIndex * 2.0 * M_PI / perimeter))));
}

// compute the one-body matrix element for a polynomial right confining potential
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// alpha = exponent of the polynomial
// cutPosition = x position below which the potential should be set to zero 
// binomials = reference to the binomial coefficients
// return value = one-body matrix element

double FQHECylinderComputePolynomialRightPotentialCoefficients (double orbitalIndex, double perimeter, int alpha, double cutPosition, BinomialCoefficients& binomials)
{  
  cutPosition -= (orbitalIndex * 2.0 * M_PI / perimeter);
  double Tmp = 0.0;
  if (cutPosition >= 0.0)
    {
      for (int i = 0; i <= alpha; ++i)
	{
#ifdef __GSL__
	  Tmp += binomials(alpha, i) * pow(cutPosition, (double) (alpha - i)) * gsl_sf_hyperg_U(-0.5 * (((double) i) - 1.0), -0.5 * (((double) i) - 1.0), cutPosition * cutPosition);
#else
	  cout << "error, GSL library is needed" << endl;
#endif
	}
      Tmp *= (0.5 * exp(- cutPosition * cutPosition) / sqrt(M_PI));
    }
  else
    {
      for (int i = 0; i <= alpha; i += 2)
	{
#ifdef __GSL__
	  Tmp += (binomials(alpha, i) * pow(-cutPosition, (double) (alpha - i)) * 
		  ((2.0 * gsl_sf_hyperg_U(-0.5 * (((double) i) - 1.0), -0.5 * (((double) i) - 1.0), 0.0))
		   - (exp(-cutPosition * cutPosition) * 
		      gsl_sf_hyperg_U(-0.5 * (((double) i) - 1.0), -0.5 * (((double) i) - 1.0), cutPosition * cutPosition))));
#else
	  cout << "error, GSL library is needed" << endl;
#endif
	}
      for (int i = 1; i <= alpha; i += 2)
	{
#ifdef __GSL__
	  Tmp += (binomials(alpha, i) * pow(-cutPosition, (double) (alpha - i)) * 
		  (exp(-cutPosition * cutPosition) * gsl_sf_hyperg_U(-0.5 * (((double) i) - 1.0), -0.5 * (((double) i) - 1.0), cutPosition * cutPosition)));
#else
	  cout << "error, GSL library is needed" << endl;
#endif
	}
      Tmp *= (0.5  / sqrt(M_PI));
    }
  return Tmp;
}

// compute the one-body matrix element for a polynomial left confining potential
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// alpha = exponent of the polynomial
// cutPosition = x position above which the potential should be set to zero 
// binomials = reference to the binomial coefficients
// return value = one-body matrix element

double FQHECylinderComputePolynomialLeftPotentialCoefficients (double orbitalIndex, double perimeter, int alpha, double cutPosition, BinomialCoefficients& binomials)
{  
  cutPosition -= (orbitalIndex * 2.0 * M_PI / perimeter);
  double Tmp = 0.0;
  if (cutPosition >= 0.0)
    {
      for (int i = 0; i <= alpha; ++i)
	{
#ifdef __GSL__
	  Tmp += binomials(alpha, i) * pow(cutPosition, (double) (alpha - i)) * gsl_sf_hyperg_U(-0.5 * (((double) i) - 1.0), -0.5 * (((double) i) - 1.0), cutPosition * cutPosition);
#else
	  cout << "error, GSL library is needed" << endl;
#endif
	}
      Tmp *= (0.5 * exp(- cutPosition * cutPosition) / sqrt(M_PI));
    }
  else
    {
      for (int i = 0; i <= alpha; i += 2)
	{
#ifdef __GSL__
	  Tmp += (binomials(alpha, i) * pow(-cutPosition, (double) (alpha - i)) * 
		  ((2.0 * gsl_sf_hyperg_U(-0.5 * (((double) i) - 1.0), -0.5 * (((double) i) - 1.0), 0.0))
		   - (exp(-cutPosition * cutPosition) * 
		      gsl_sf_hyperg_U(-0.5 * (((double) i) - 1.0), -0.5 * (((double) i) - 1.0), cutPosition * cutPosition))));
#else
	  cout << "error, GSL library is needed" << endl;
#endif
	}
      for (int i = 1; i <= alpha; i += 2)
	{
#ifdef __GSL__
	  Tmp += (binomials(alpha, i) * pow(-cutPosition, (double) (alpha - i)) * 
		  (exp(-cutPosition * cutPosition) * gsl_sf_hyperg_U(-0.5 * (((double) i) - 1.0), -0.5 * (((double) i) - 1.0), cutPosition * cutPosition)));
#else
	  cout << "error, GSL library is needed" << endl;
#endif
	}
      Tmp *= (0.5  / sqrt(M_PI));
    }
  return Tmp;
}

// compute the one-body matrix element for a polynomial right confining potential defined in momentum space
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// alpha = exponent of the polynomial
// cutPosition = x position below which the potential should be set to zero 
// return value = one-body matrix element

double FQHECylinderComputePolynomialRightMomentumPotentialCoefficients (double orbitalIndex, double perimeter, int alpha, double cutPosition)
{
  orbitalIndex -= cutPosition * perimeter / (2.0 * M_PI);
  if (orbitalIndex >= 0.0)
    {
      return pow(orbitalIndex, alpha);
    }
  else
    {
      return 0.0;
    }
}

// compute the one-body matrix element for a polynomial left confining potential
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// alpha = exponent of the polynomial
// cutPosition = x position above which the potential should be set to zero 
// return value = one-body matrix element

double FQHECylinderComputePolynomialLeftMomentumPotentialCoefficients (double orbitalIndex, double perimeter, int alpha, double cutPosition)
{
  orbitalIndex -= cutPosition * perimeter / (2.0 * M_PI);
  orbitalIndex *= -1.0;
  if (orbitalIndex >= 0.0)
    {
      return pow(orbitalIndex, alpha);
    }
  else
    {
      return 0.0;
    }
}
