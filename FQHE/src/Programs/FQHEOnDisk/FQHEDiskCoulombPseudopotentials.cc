#include "Options/Options.h"

#include "Tools/FQHESpectrum/PseudoPotentials.h"
#include "Tools/FQHESpectrum/AbstractZDensityProfile.h"

#include "Vector/RealVector.h"
 
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "GeneralTools/FilenameTools.h"


#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

using std::ofstream;
using std::ios;
using std::cout;
using std::endl;

// compute the Haldane pseudopotentials on the disk for the Coulomb interaction in the lowest Landau level
//
// lzMax = maximum angular momentum to evaluate
// return value = array that contains lzMax + 1 pseudopotentials
double* FQHEDiskCoulombPseudopotentialLLL (int lzMax);

// compute the Haldane pseudopotentials on the disk for the Coulomb interaction in the second Landau level
//
// lzMax = maximum angular momentum to evaluate
// return value = array that contains lzMax + 1 pseudopotentials
double* FQHEDiskCoulombPseudopotential2LL (int lzMax);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEDiskCoulombPseudopotentials" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;


  (*SystemGroup) += new SingleIntegerOption  ('l', "landau-level", "index of the Landau level (0 for the lowest Landau level)", 0, true, 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "lz-max", "maximum Lz angular momentum that has to be evaluated", 8);
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name (default is pseudopotential_coulomb_l_x_2s_y[_v_xxx_w_yyy].dat)");
  (*SystemGroup) += new BooleanOption ('\n', "std-output", "use standard output instead of an output file");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskCoulombPseudopotentials -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int LandauLevel = Manager.GetInteger("landau-level");
  int MaxMomentum = Manager.GetInteger("lz-max");
      
  char* OutputFile;

  if (Manager.GetString("output") == 0l)
    {
      OutputFile = Manager.GetFormattedString("pseudopotential_disk_coulomb_l_%landau-level%_2s_%lz-max%.dat");
    }
  else
    {
      OutputFile = new char [strlen(Manager.GetString("output")) + 1];
      strcpy (OutputFile, Manager.GetString("output"));
    }
    
  double* Pseudopotentials = 0;
  if (LandauLevel == 0)
    {
      Pseudopotentials = FQHEDiskCoulombPseudopotentialLLL(MaxMomentum);
    }
  else
    {
      if (LandauLevel == 1)
	{
	  Pseudopotentials = FQHEDiskCoulombPseudopotential2LL(MaxMomentum);
	}
      else
	{
	  cout << "Landau levels higher than the second one are not supported" << endl;
	  return -1;
	}
    }
  if (Manager.GetBoolean("std-output") == false)
    {
      ofstream File;
      File.open(OutputFile, ios::binary | ios::out);
      File.precision(14);
      File << "# pseudopotentials on the disk for coulomb interaction " << endl
	   << "# in the Landau level N=" << LandauLevel << " and up to angular momentum Lz=" << MaxMomentum << endl;
      File << "#" << endl
	   << "# Pseudopotentials = V_0 V_1 ..." << endl << endl
	   << "Pseudopotentials =";
      for (int i = 0; i <= MaxMomentum; ++i)
	File << " " << Pseudopotentials[i];
      File.close();
    }
  else
    {
      cout.precision(14);
      cout << "# pseudopotentials on the disk for coulomb interaction "
	   << "# in the Landau level N=" << LandauLevel << " and up to angular momentum Lz=" << MaxMomentum << endl;
      cout << "#" << endl
	   << "# Pseudopotentials = V_0 V_1 ..." << endl << endl
	   << "Pseudopotentials =";
      for (int i = 0; i <= MaxMomentum; ++i)
	cout << " " << Pseudopotentials[i];
      cout << endl;
    }

  return 0;
}

// compute the Haldane pseudopotentials on the disk for the Coulomb interaction in the lowest Landau level
//
// lzMax = maximum angular momentum to evaluate
// return value = array that contains lzMax + 1 pseudopotentials

double* FQHEDiskCoulombPseudopotentialLLL (int lzMax)
{
  double* TmpPseudopotentials = new double [lzMax];
  double Gamma1 = sqrt(M_PI);
  double Gamma2 = 1.0;
  TmpPseudopotentials[0] = 0.5 * Gamma1 / Gamma2;
  for (int i = 1; i <= lzMax; ++i)
    {
      Gamma1 *= ((double) i) - 0.5;
      Gamma2 *= (double) i; 
      TmpPseudopotentials[i] = 0.5 * Gamma1 / Gamma2;      
    }
  return TmpPseudopotentials;
}

// compute the Haldane pseudopotentials on the disk for the Coulomb interaction in the second Landau level
//
// lzMax = maximum angular momentum to evaluate
// return value = array that contains lzMax + 1 pseudopotentials

double* FQHEDiskCoulombPseudopotential2LL (int lzMax)
{
  double* TmpPseudopotentials = new double [lzMax];
  double Gamma1 = sqrt(M_PI);
  double Gamma2 = 2.0;
  TmpPseudopotentials[0] = 11.0 * sqrt(M_PI) / 32.0;
  TmpPseudopotentials[1] = 15.0 * sqrt(M_PI) / 64.0;
  TmpPseudopotentials[2] = 65.0 * sqrt(M_PI) / 256.0;
  for (int i = 3; i <= lzMax; ++i)
    {
      Gamma1 *= ((double) i) - 2.5;
      Gamma2 *= (double) i; 
      TmpPseudopotentials[i] = (8.0 * ((double) i) - 11.0) * (8.0 * ((double) i) - 3.0)  * Gamma1 / (128.0 * Gamma2);      
    }
  return TmpPseudopotentials;
}
