#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// compute the weight of a given orbital for a sharp real space cut
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutPosition = x position of the cut
// return value = square of the orbital weight
double FQHECylinderComputeSharpRealSpaceCutCoefficient (double orbitalIndex, double perimeter, double cutPosition);

// compute the weight of a given orbital for a sharp real space cut with patch having a finite length along the cylinder axis
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutInitialPosition = leftmost x position of the cut
// cutFinalPosition = leftmost x position of the cut
// return value = square of the orbital weight
double FQHECylinderComputeSharpRealSpaceCutCoefficient (double orbitalIndex, double perimeter, double cutInitialPosition, double cutFinalPosition);

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHECylinderRealSpacePartitionCoefficients" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (if zero, assume an infinite cylinder)", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "error", "error below which a coefficient is consider as 0 (or 1)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the cylinder", 1);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
  (*SystemGroup) += new BooleanOption ('\n', "finite-patch", "use a patch with a finite length and optionally a finite width");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "patch-center", "position of the patch center along the cylinder axis", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "patch-length", "length of the patch along the cylinder axis", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "patch-width", "width of the patch along the cylinder perimeter (0 if it should extend along the whole cylinder perimeter)", 0.0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "optional output file name (default is realspace_cylinder_l_*_perimeter_*_2s_*.dat)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderRealSpacePartitionCoefficients -h" << endl;
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
  char* OutputFile = 0;

  if ((Manager.GetBoolean("finite-patch") == false) || (Manager.GetDouble("patch-width") == 0.0))
    {
      double CutPosition = 0.0;
      double CutLength = 0.0;
      if (Manager.GetString("output-file") == 0)
	{
	  OutputFile = new char[512];
	  sprintf (OutputFile, "realspace_cylinder_l_%.6f_perimeter_%.6f_2s_%d.dat", CutPosition, Perimeter, NbrFluxQuanta);
	}
      else
	{
	  OutputFile = new char[strlen(Manager.GetString("output-file")) + 1];
	  strcpy (OutputFile, Manager.GetString("output-file"));
	}

      ofstream File;
      File.open(OutputFile, ios::binary | ios::out);
      File.precision(14);
      if (Manager.GetBoolean("finite-patch") == false)
	{
	  File << "# real space coefficients for a sharp cut at x=" << CutPosition << " on ";
	}
      else
	{
	  CutPosition = Manager.GetDouble("patch-center");
	  CutLength = Manager.GetDouble("patch-length");
	  File << "# real space coefficients for a sharp cut centered at x=" << CutPosition << " and of length " << CutLength << " on ";
	}
      if (NbrFluxQuanta == 0)
	{
	  File << "an infinite cylinder with perimeter L=" << Perimeter;
	}
      else
	{
	  File << "a cylinder with perimeter L=" << Perimeter << " and N_phi=" << NbrFluxQuanta;
	}
      File << endl << "OrbitalSquareWeights =";
      int NbrCoefficients = 0;
      double* Coefficients = 0;
      if (NbrFluxQuanta == 0)
	{
	  int MaxNbrCoefficients = 1000; 
	  double* TmpCoefficients = new double [MaxNbrCoefficients];
	  if (Manager.GetBoolean("finite-patch") == false)
	    {
	      for (int i = 0; i < MaxNbrCoefficients; ++i)
		{
		  TmpCoefficients[i] = FQHECylinderComputeSharpRealSpaceCutCoefficient(i, Perimeter, CutPosition);
		  if (TmpCoefficients[i] < Error)
		    {
		      i = MaxNbrCoefficients;
		    }
		  else
		    {
		      ++NbrCoefficients;
		    }	    
		}
	    }
	  else
	    {
	      for (int i = 0; i < MaxNbrCoefficients; ++i)
		{
		  TmpCoefficients[i] = FQHECylinderComputeSharpRealSpaceCutCoefficient(i, Perimeter, CutPosition - (0.5 * CutLength), 
										       CutPosition + (0.5 * CutLength));
		  if (TmpCoefficients[i] < Error)
		    {
		      i = MaxNbrCoefficients;
		    }
		  else
		    {
		      ++NbrCoefficients;
		    }
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
	  if (Manager.GetBoolean("finite-patch") == false)
	    {
	      for (NbrCoefficients = 0; NbrCoefficients <= NbrFluxQuanta; ++NbrCoefficients)
		{
		  Coefficients[NbrCoefficients] = FQHECylinderComputeSharpRealSpaceCutCoefficient(((double) NbrCoefficients) - 0.5 * ((double) NbrFluxQuanta), Perimeter, CutPosition);
		}
	    }
	  else
	    {
	      for (NbrCoefficients = 0; NbrCoefficients <= NbrFluxQuanta; ++NbrCoefficients)
		{
		  Coefficients[NbrCoefficients] = FQHECylinderComputeSharpRealSpaceCutCoefficient(((double) NbrCoefficients) - 0.5 * ((double) NbrFluxQuanta), Perimeter, 
												  CutPosition - (0.5 * CutLength), CutPosition + (0.5 * CutLength));
		}
	    }
	}
      for (int i = 0; i < NbrCoefficients; ++i)
	File << " " << Coefficients[i];
      File << endl;
      File.close();
      delete[] Coefficients;
    }
  return 0;
}


// compute the weight of a given orbital for a sharp real space cut
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutPosition = x position of the cut
// return value = square of the orbital weight

double FQHECylinderComputeSharpRealSpaceCutCoefficient (double orbitalIndex, double perimeter, double cutPosition)
{  
  return (0.5 * (1.0 + erf(cutPosition - (orbitalIndex * 2.0 * M_PI / perimeter))));
}

// compute the weight of a given orbital for a sharp real space cut with patch having a finite length along the cylinder axis
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutInitialPosition = leftmost x position of the cut
// cutFinalPosition = leftmost x position of the cut
// return value = square of the orbital weight

double FQHECylinderComputeSharpRealSpaceCutCoefficient (double orbitalIndex, double perimeter, double cutInitialPosition, double cutFinalPosition)
{  
  return (0.5 * (erf((orbitalIndex * 2.0 * M_PI / perimeter) - cutInitialPosition) + erf(cutFinalPosition - (orbitalIndex * 2.0 * M_PI / perimeter))));
}
