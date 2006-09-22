#include "HilbertSpace/BosonOnSphere.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension
long EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension
long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension
long FermionEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension
long FermionShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);


int main(int argc, char** argv)
{
  OptionManager Manager ("GetBosonsDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta", 20);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "ground-only", "get the dimension only for the largest subspace");
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_torus_n_nbrparticles_q_nbrfluxquanta.dim");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereGetDimension -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrFluxQuanta = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger(); 
  int  LzMin = 0;
  if (((NbrParticles * NbrFluxQuanta) & 1) != 0)
    LzMin = 1;
  if (((BooleanOption*) Manager["ground-only"])->GetBoolean() == true)
    {
      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	cout << EvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, LzMin) << endl;
      else
	cout << FermionEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, LzMin) << endl;
    }
  else
    {
      int LzMax = 0;
      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	LzMax = (NbrParticles * NbrFluxQuanta);
      else
	LzMax = ((NbrFluxQuanta - NbrParticles + 1) * NbrParticles);
      long* LzDimensions = new long [1 + ((LzMax - LzMin) >> 1)];
      long* LDimensions = new long [1 + ((LzMax - LzMin) >> 1)];
      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	for (int x = LzMin; x <= LzMax; x += 2)
	  LzDimensions[(x - LzMin) >> 1] = EvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, x);
      else
	for (int x = LzMin; x <= LzMax; x += 2)
	  LzDimensions[(x - LzMin) >> 1] =  FermionEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, x);
      LDimensions[(LzMax - LzMin) >> 1] = LzDimensions[(LzMax - LzMin) >> 1];
      long TotalDimension = LzDimensions[(LzMax - LzMin) >> 1];
      for (int x = LzMax - 2; x >= LzMin; x -= 2)
	{
	  LDimensions[(x - LzMin) >> 1] =  LzDimensions[(x - LzMin) >> 1] - LzDimensions[((x - LzMin) >> 1) + 1];
	  TotalDimension += LzDimensions[(x - LzMin) >> 1];
	}
            
      if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
	{
	  char* OutputFileName = 0;
	  if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
	    {
	      OutputFileName = new char[256];
	      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
		sprintf (OutputFileName, "bosons_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
	      else
		sprintf (OutputFileName, "fermions_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
	    }
	  else
	    {
	      OutputFileName = new char[strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
	      strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
	    }
	  ofstream File;
	  File.open(OutputFileName, ios::binary | ios::out);
	  File << "# Hilbert space dimension in each L and Lz sector for " << NbrParticles << " ";
	  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	    File << "bosons";
	  else
	    File << "femions";
	  File << " on the sphere geometry with " << NbrFluxQuanta << " flux quanta" << endl;
	  File << "# total Hilbert space dimension = " << TotalDimension << endl << endl 
	       << "N = " << NbrParticles << endl
	       << "2S = " << NbrFluxQuanta << endl << endl
	       << "#  dimensions for the Lz subspaces (starting from 2Lz = " << LzMin << " " << (LzMin + 2) << " ..." << endl
	       << "Lz =";
	  for (int x = LzMin; x <= LzMax; x += 2)
	    File << " " << LzDimensions[(x - LzMin) >> 1];
	  File << endl
	       << "#  dimensions for the L subspaces (starting from 2L = " << LzMin << " " << (LzMin + 2) << " ..." << endl
	       << "L =";
	  for (int x = LzMin; x <= LzMax; x += 2)
	    File << " " << LDimensions[(x - LzMin) >> 1];
	  File << endl;
	  File.close();
	  delete[] OutputFileName;
	}
      else
	{
	  cout << "Lz =";
	  for (int x = LzMin; x <= LzMax; x += 2)
	    cout << " " << LzDimensions[(x - LzMin) >> 1];
	  cout << endl << "L =";
	  for (int x = LzMin; x <= LzMax; x += 2)
	    cout << " " << LDimensions[(x - LzMin) >> 1];	  
	  cout << endl;
	}
      delete[] LzDimensions;
      delete[] LDimensions;
    }
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

long EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  return ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, (totalLz + lzMax * nbrBosons) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons == 0) || ((nbrBosons * lzMax) < totalLz))
    return 0l;
  if (((nbrBosons * lzMax) == totalLz) || (lzMax == 0) || (totalLz == 0))
    {
      return 1l;
    }
  long TmpDim = 0;
  while ((totalLz >= 0) && (nbrBosons > 0))
    {
      TmpDim += ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz);
      --nbrBosons;
      totalLz -= lzMax;
    }
  return TmpDim;
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  return FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + nbrFermions * lzMax) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long FermionShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return (long) 0;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return (long) 0;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return (long) 1;
  if (LzTotalMax == totalLz)
    return (long) 1;
  return  (FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax)
	   +  FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz));
}
