#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnDisk.h"
#include "HilbertSpace/FermionOnDisk.h"
#include "HilbertSpace/FermionOnDiskUnlimited.h"
#include "HilbertSpace/ParticleOnDisk.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Operator/ParticleOnDiskDensityOperator.h"
#include "FunctionBasis/ParticleOnDiskFunctionBasis.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEDensity" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PlotOptionGroup = new OptionGroup ("plot options");  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += PlotOptionGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (overriding the one found in the vector file name if greater than 0)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "momentum", "single particle momentum (overriding the one found in the vector file name if greater than 0)", 0, true, 0);
  (*SystemGroup) += new BooleanOption  ('\n', "bosons", "use boson statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "fermions", "use fermion statistics");
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the vector file describing the state whose density has to be plotted");

  (*PlotOptionGroup) += new SingleStringOption ('\n', "output", "output file name", "density.dat");
  (*PlotOptionGroup) += new SingleIntegerOption ('\n', "nbr-samples", "number of samples in radial direction", 1000, true, 10);
  (*PlotOptionGroup) += new BooleanOption  ('\n', "profile", "only draw density profile");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEDensity -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "QHEBosonsCorrelation requires a state" << endl;
      return -1;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["momentum"])->GetInteger();
  bool BosonFlag = ((BooleanOption*) Manager["bosons"])->GetBoolean();
  bool FermionFlag = ((BooleanOption*) Manager["fermions"])->GetBoolean();
  bool ProfileFlag = ((BooleanOption*) Manager["profile"])->GetBoolean();
  char* OutputName = ((SingleStringOption*) Manager["output"])->GetString();
  int NbrSamples = ((SingleIntegerOption*) Manager["nbr-samples"])->GetInteger();
  
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["state"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
      return -1;      
    }

  char* StrNbrParticles;
  if (NbrParticles == 0)
    {
      StrNbrParticles = strstr(((SingleStringOption*) Manager["state"])->GetString(), "_n_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 3;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      NbrParticles = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess number of particles from file name " << ((SingleStringOption*) Manager["state"])->GetString() << endl
	       << "use --nbr-particles option" << endl;
	  return -1;            
	}
    }
  if (Lz == 0)
    {
      StrNbrParticles = strstr(((SingleStringOption*) Manager["state"])->GetString(), "_lz_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 4;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '.') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      Lz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess momentum from file name " << ((SingleStringOption*) Manager["state"])->GetString() << endl
	       << "use --momentum option" << endl;
	  return -1;            
	}
    }
  if ((FermionFlag == false) && (BosonFlag == false))
    {
       if (strncmp(((SingleStringOption*) Manager["state"])->GetString(), "bosons", 6) == 0)
	 BosonFlag = true;
       else
	 if (strncmp(((SingleStringOption*) Manager["state"])->GetString(), "bosons", 6) == 0)
	   FermionFlag = true;
	 else
	   {
	     cout << "can't guess statistics from file name " << ((SingleStringOption*) Manager["state"])->GetString() << endl
		  << "use --bosons or --fermions option" << endl;
	     return -1;            
	   }
    }
  ParticleOnDiskFunctionBasis Basis(Lz);
  ParticleOnDisk* Space = 0;
  if (BosonFlag == true)
    {
      Space = new BosonOnDisk(NbrParticles, Lz);
    }
  else
    {
#ifdef __64_BITS__
      if ((Lz - (((NbrParticles - 1) * (NbrParticles - 2)) / 2)) < 63)      
#else
      if ((Lz - (((NbrParticles - 1) * (NbrParticles - 2)) / 2)) < 31)
#endif
	Space = new FermionOnDisk (NbrParticles, Lz);
      else
	Space = new FermionOnDiskUnlimited (NbrParticles, Lz);      
    }
  Complex* PrecalculatedValues = new Complex [Lz + 1];
	  
  for (int i = 0; i <= Lz; ++i)
    {
      ParticleOnDiskDensityOperator Operator (Space, i);
      PrecalculatedValues[i] = Operator.MatrixElement(State, State);
    }
  
  if (ProfileFlag == true)
    {
      Complex Tmp (0.0);
      double RInc = 1.2 * ((double) Lz) / ((double) NbrSamples);
      ofstream File;
      File.precision(14);
      File.open(OutputName, ios::binary | ios::out);
      NbrSamples *= 2;      
      RealVector Value(2, true);
      Value[0] = - 1.2 * ((double) Lz);
      Complex TmpValue;
      for (int i = 0; i <= NbrSamples; ++i)
	{
	  Tmp = 0.0;
	  for (int j = 0; j <= Lz; ++j)
	    {
	      Basis.GetFunctionValue(Value, TmpValue, j);
	      Tmp += PrecalculatedValues[j] * SqrNorm(TmpValue);
	    }
	  File << Value[0] << " " << (Tmp.Re * exp (-0.5 * (Value[0] * Value[0]))) << endl;
	  Value[0] += RInc;
	}
      File.close();
    }

  delete[] PrecalculatedValues;
  delete Space;
}
