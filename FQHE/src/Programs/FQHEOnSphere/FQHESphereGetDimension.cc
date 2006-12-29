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

// evaluate Hilbert space dimension for bosons
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension
long BosonEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for bosons
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension
long BosonShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

// evaluate Hilbert space dimension for fermions
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension
long FermionEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for fermions
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension
long FermionShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

// evaluate Hilbert space dimension for fermions with SU(2) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value (with shift nbrFermions * lzMax)
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension
long FermionSU2ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin);

// evaluate Hilbert space dimension for fermions with SU(2)xSU(2) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value (with shift nbrFermions * lzMax)
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// return value = Hilbert space dimension
long FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin);

// evaluate Hilbert space dimension for fermions with SU(4) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// totalEntanglement = number of particles with entanglement plus
// return value = Hilbert space dimension
long FermionSU4ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin, int totalEntanglement);

// save dimensions in a given file
//
// outputFileName = output file name
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true for bosons, false for fermions
// lzDimensions = array taht contains dimension of each Lz sector
// lDimensions = array that contains dimension of each L sector (without taking into account the Lz degeneracy)
// lzMin = twice the minimum Lz value
// lzMax = twice the maximum Lz value
// totalDimension = total Hilbert space dimension
bool WriteDimensionToDisk(char* outputFileName, int nbrParticles, int nbrFluxQuanta, bool statistics,
			  long* lzDimensions, long* lDimensions, int lzMin, int lzMax, long totalDimension);


// save dimensions in a given file (for particles with SU(2) spin)
//
// outputFileName = output file name
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true for bosons, false for fermions
// lzDimensions = array that contains dimension of each Lz sector per Sz
// lDimensions = array that contains dimension of each L sector per Sz (without taking into account the Lz degeneracy)
// lzMin = twice the minimum Lz value for each value of Sz
// lzMax = twice the maximum Lz value for each value of Sz
// sz = minal number of particles with spin up
// totalDimension = total Hilbert space dimension
bool SU2WriteDimensionToDisk(char* outputFileName, int nbrParticles, int nbrFluxQuanta, bool statistics,
			     long** lzDimensions, long** lDimensions, int* lzMin, int* lzMax, long totalDimension, int sz);

// save dimensions in a given file (for particles with SU(4) spin)
//
// outputFileName = output file name
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true for bosons, false for fermions
// lzDimensions = array that contains dimension of each Lz sector per Sz and per Iz
// lDimensions = array that contains dimension of each L sector per Sz and per Iz (without taking into account the Lz degeneracy)
// lzMin = twice the minimum Lz value for each value of Sz and Iz
// lzMax = twice the maximum Lz value for each value of Sz and Iz
// sz = minal number of particles with spin up
// iz = minal number of particles with spin plus
// totalDimension = total Hilbert space dimension
bool SU4WriteDimensionToDisk(char* outputFileName, int nbrParticles, int nbrFluxQuanta, bool statistics,
			     long*** lzDimensions, long*** lDimensions, int** lzMin, int** lzMax, long totalDimension, int sz, int iz);



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
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "su2su2-spin", "consider particles with SU(2)xSU(2) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "su4-spin", "consider particles with SU(4) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "ground-only", "get the dimension only for the largest subspace");
  (*SystemGroup) += new BooleanOption ('\n', "use-files", "use dimension files that have been previously generated to increase speed. Files must be in current directory and obey the statistics_sphere_n_nbrparticles_q_nbrfluxquanta.dim naming convention");
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_sphere_n_nbrparticles_q_nbrfluxquanta.dim");
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
  if ((((BooleanOption*) Manager["su4-spin"])->GetBoolean() == false) && (((BooleanOption*) Manager["su2-spin"])->GetBoolean() == false) && (((BooleanOption*) Manager["su2su2-spin"])->GetBoolean() == false))
    if (((BooleanOption*) Manager["ground-only"])->GetBoolean() == true)
      {
	if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	  cout << BosonEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, LzMin) << endl;
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
	    LzDimensions[(x - LzMin) >> 1] = BosonEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, x);
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
	    WriteDimensionToDisk (OutputFileName, NbrParticles, NbrParticles, ((BooleanOption*) Manager["boson"])->GetBoolean(),
				  LzDimensions, LDimensions, LzMin, LzMax, TotalDimension);
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
  else
    {
      int Sz = 0;
      if (NbrParticles & 1)
	Sz = 1;
      int NbrSz = NbrParticles;//(NbrParticles + 2 - Sz) >> 1;
      if (((BooleanOption*) Manager["su4-spin"])->GetBoolean() == true)
	{
	  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	    {
	      cout << "SU(4) mode for bosons not yet available" << endl;	
	      return -1;
	    }
	  if (((BooleanOption*) Manager["ground-only"])->GetBoolean() == true)
	    cout << FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, (LzMin + (NbrFluxQuanta * NbrParticles)) >> 1, 
								   (Sz + NbrParticles) >> 1, (Sz + NbrParticles) >> 1)<< endl;
	  else
	    {
	      	      
	      int** FullLzMax = new int* [NbrSz];
	      int** FullLzMin = new int* [NbrSz];
	      long*** LzDimensions = new long**[NbrSz];
	      long*** LDimensions = new long**[NbrSz];
	      for (int i = 0; i < NbrSz; ++i)
		{
		  FullLzMax[i] = new int[NbrSz];
		  FullLzMin[i] = new int[NbrSz];
		  LzDimensions[i] = new long*[NbrSz];
		  LDimensions[i] = new long*[NbrSz];
		  for (int j = 0; j < NbrSz; ++j)
		    {
		      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
			{
			  FullLzMin[i][j] = 0;
			  FullLzMax[i][j] = (NbrParticles * NbrFluxQuanta);
			}
		      else
			{
			  FullLzMin[i][j] = ((NbrParticles - 1) * NbrParticles);
			  FullLzMax[i][j] = ((NbrFluxQuanta - NbrParticles + 1) * NbrParticles);
			}
		      LzDimensions[i][j] = new long [1 + ((FullLzMax[i][j] - FullLzMin[i][j]) >> 1)];
		      LDimensions[i][j] = new long [1 + ((FullLzMax[i][j] - FullLzMin[i][j]) >> 1)];
		    }
		}
	      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
		for (int MinSz = 0; MinSz < NbrSz; ++MinSz)
		  for (int MinIz = 0; MinIz < NbrSz; ++MinIz)
		    for (int x = FullLzMin[MinSz][MinIz]; x <= FullLzMax[MinSz][MinIz]; x += 2)
		      LzDimensions[MinSz][MinIz][(x - FullLzMin[MinSz][MinIz]) >> 1] = BosonEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, x);
	      else
		for (int MinSz = 0; MinSz < NbrSz; ++MinSz)
		  for (int MinIz = 0; MinIz < NbrSz; ++MinIz)
		    {
		      cout << MinSz << " " << MinIz << " " << NbrSz << endl;
		      for (int x = FullLzMin[MinSz][MinIz]; x <= FullLzMax[MinSz][MinIz]; x += 2)
			LzDimensions[MinSz][MinIz][(x - FullLzMin[MinSz][MinIz]) >> 1] =  FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, x, MinSz, MinIz);
		    }
	      long TotalDimension = 0l;
	      for (int MinSz = 0; MinSz < NbrSz; ++MinSz)
		for (int MinIz = 0; MinIz < NbrSz; ++MinIz)
		  {
		    LDimensions[MinSz][MinIz][(FullLzMax[MinSz][MinIz] - FullLzMin[MinSz][MinIz]) >> 1] = LzDimensions[MinSz][MinIz][(FullLzMax[MinSz][MinIz] - FullLzMin[MinSz][MinIz]) >> 1];
		    TotalDimension += LzDimensions[MinSz][MinIz][(FullLzMax[MinSz][MinIz] - FullLzMin[MinSz][MinIz]) >> 1];
		  }
	      for (int MinSz = 0; MinSz < NbrSz; ++MinSz)
		for (int MinIz = 0; MinIz < NbrSz; ++MinIz)
		  for (int x = FullLzMax[MinSz][MinIz] - 2; x >= FullLzMin[MinSz][MinIz]; x -= 2)
		    {
		      LDimensions[MinSz][MinIz][(x - FullLzMin[MinSz][MinIz]) >> 1] =  (LzDimensions[MinSz][MinIz][(x - FullLzMin[MinSz][MinIz]) >> 1] 
											- LzDimensions[MinSz][MinIz][((x - FullLzMin[MinSz][MinIz]) >> 1) + 1]);
		      TotalDimension += LzDimensions[MinSz][MinIz][(x - FullLzMin[MinSz][MinIz]) >> 1];
		    }	      
	      if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
		{
		  char* OutputFileName = 0;
		  if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
		    {
		      OutputFileName = new char[256];
		      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
			sprintf (OutputFileName, "bosons_sphere_su4_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		      else
			sprintf (OutputFileName, "fermions_sphere_su4_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		    }
		  else
		    {
		      OutputFileName = new char[strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
		      strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
		    }
		  SU4WriteDimensionToDisk (OutputFileName, NbrParticles, NbrFluxQuanta, ((BooleanOption*) Manager["boson"])->GetBoolean(),
					   LzDimensions, LDimensions, FullLzMin, FullLzMax, TotalDimension, Sz, Sz);
		  delete[] OutputFileName;
		}
	      else
		{
		  for (int MinSz = Sz; MinSz <= NbrParticles; ++MinSz)
		    for (int MinIz = Sz; MinIz <= NbrParticles; ++MinIz)
		      for (int x = FullLzMin[MinSz][MinIz]; x <= FullLzMax[MinSz][MinIz]; x += 2)
			cout << x << " " << ((2 * MinSz) - NbrParticles)<< " " << ((2 * MinIz) - NbrParticles) << " " 
			     << LzDimensions[MinSz][MinIz][(x - FullLzMin[MinSz][MinIz]) >> 1] << " " << LDimensions[MinSz][MinIz][(x - FullLzMin[MinSz][MinIz]) >> 1] << endl;
		}
	      for (int i = 0; i < NbrSz; ++i)
		{		 
		  for (int j = 0; j < NbrSz; ++j)
		    {
		      delete[] LzDimensions[i][j];
		      delete[] LDimensions[i][j];
		    }
		  delete[] FullLzMin[i];
		  delete[] FullLzMax[i];
		  delete[] LzDimensions[i];
		  delete[] LDimensions[i];
		}
	      delete[] FullLzMin;
	      delete[] FullLzMax;
	      delete[] LzDimensions;
	      delete[] LDimensions;	      
	    }
	}
      else
	{
	  cout << "SU(2) mode not yet available" << endl;	
	  return -1;
	}
    }

}

// evaluate Hilbert space dimension for bosons
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

long BosonEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  return BosonShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, (totalLz + lzMax * nbrBosons) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for bosons
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long BosonShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
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
      TmpDim += BosonShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz);
      --nbrBosons;
      totalLz -= lzMax;
    }
  return TmpDim;
}

// evaluate Hilbert space dimension for fermions
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  return FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + nbrFermions * lzMax) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for fermions
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

// evaluate Hilbert space dimension for fermions with SU(2) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value (with shift nbrFermions * lzMax)
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension

long FermionSU2ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalSpin > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz))
    return 0l;
    
  if (nbrFermions == 1) 
    if (lzMax >= totalLz)
      return 1l;
    else
      return 0l;

  if ((lzMax == 0)  && (totalLz != 0))
    return 0l;

  if ((nbrFermions == 2) && (totalLz == (2 * lzMax)))
    if (totalSpin == 1)
      return 1l;
    else
      return 0l;
  else
    return  (FermionSU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1)
	     + FermionSU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin)
	     + FermionSU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, totalSpin));
}

// evaluate Hilbert space dimension for fermions with SU(2)xSU(2) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value (with shift nbrFermions * lzMax)
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// return value = Hilbert space dimension

long FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalIsospin < 0) ||  
      (totalSpin > nbrFermions) || (totalIsospin > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < totalIsospin) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalIsospin)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz) || ((((2 * lzMax + nbrFermions + 1 - totalIsospin) * nbrFermions) >> 1) < totalLz))
    return 0l;
    
  if (nbrFermions == 1) 
    if (lzMax >= totalLz)
      return 1l;
    else
      return 0l;

  if ((lzMax == 0)  && (totalLz != 0))
    return 0l;

  unsigned long Tmp = 0l;
  if (nbrFermions >= 3)    
    {
      Tmp += (FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 2, totalIsospin - 1)
	      + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 2)
	      + (2l * FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 1))
	      + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin)
	      + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin, totalIsospin - 1));

      if (nbrFermions > 3)
	{
	  Tmp += (FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 2)
		  + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 1)
		  + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 2)
		  + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 1));
	  if (nbrFermions == 4)
	    {
	      if ((totalLz == (4 * lzMax)) && (totalSpin == 2) && (totalIsospin == 2))
		++Tmp;      
	    }
	  else
	    Tmp += FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 4, lzMax - 1, totalLz - (4 * lzMax), totalSpin - 2, totalIsospin - 2);
	}
      else
	if ((totalLz == (3 * lzMax)) && (((totalSpin == 2) || (totalSpin == 1)) && ((totalIsospin == 2) || (totalIsospin == 1))))
	  ++Tmp;
    }
  else
    if (totalLz == (2 * lzMax))
      {
 	switch (totalSpin)
 	  {
 	  case 2:
	    if (totalIsospin == 1)
	      ++Tmp;
 	    break;
 	  case 1:
	    switch (totalIsospin)
	      {
	      case 2:
		++Tmp;
		break;
	      case 1:
		Tmp += 2l;
		break;
	      case 0:
		++Tmp;
		break;
	      }
	    break;
 	  case 0:
	    if (totalIsospin == 1) 
	      ++Tmp;
	    break; 
 	  }
      }

  return  (Tmp + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin)
	   + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin - 1)
	   + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin - 1)
	   + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin)
	   + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, totalSpin, totalIsospin));
}

// evaluate Hilbert space dimension for fermions with SU(4) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// totalEntanglement = number of particles with entanglement plus
// return value = Hilbert space dimension

long FermionSU4ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin, int totalEntanglement)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalIsospin < 0) || (totalEntanglement < 0) ||
      (totalSpin > nbrFermions) || (totalIsospin > nbrFermions) || (totalEntanglement > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < totalIsospin) || ((2 * (lzMax + 1)) < totalEntanglement) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalIsospin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalEntanglement)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz) 
      || ((((2 * lzMax + nbrFermions + 1 - totalIsospin) * nbrFermions) >> 1) < totalLz)
      || ((((2 * lzMax + nbrFermions + 1 - totalEntanglement) * nbrFermions) >> 1) < totalLz))
    return 0l;
    
  if ((nbrFermions == 0) && (totalLz == 0) && (totalSpin == 0) && (totalIsospin == 0) && (totalEntanglement == 0))
    return 1l;
  if (nbrFermions == 1) 
    if ((lzMax >= totalLz) && (totalEntanglement != (totalSpin ^ totalIsospin)))
      return 1l;
    else
      return 0l;

  if ((lzMax == 0)  && (totalLz != 0))
    return 0l;

  unsigned long Tmp = 0l;
  if (nbrFermions >= 3)    
    {
      Tmp += (FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 2, totalIsospin - 1, totalEntanglement - 1)
	      + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 2, totalEntanglement - 1)
	      + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 1, totalEntanglement - 2)
	      + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 1, totalEntanglement)
	      + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin, totalEntanglement - 1)
	      + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin, totalIsospin - 1, totalEntanglement - 1));
      
      if (nbrFermions > 3)
	{
 	  Tmp += (FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 2, totalEntanglement - 1)
 		  + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 1, totalEntanglement - 1)
 		  + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 2, totalEntanglement - 2)
 		  + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 1, totalEntanglement - 2));
 	  if (nbrFermions == 4)
 	    {
 	      if ((totalLz == (4 * lzMax)) && (totalSpin == 2) && (totalIsospin == 2) && (totalEntanglement == 2))
 		++Tmp;      
 	    }
 	  else
 	    Tmp += FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 4, lzMax - 1, totalLz - (4 * lzMax), totalSpin - 2, totalIsospin - 2, totalEntanglement - 2);
 	}
      else
	if ((totalLz == (3 * lzMax)) && (((totalSpin * totalIsospin * totalEntanglement) == 4) || ((totalSpin * totalIsospin * totalEntanglement) == 1)))
 	  ++Tmp;
    }
  else
    if (totalLz == (2 * lzMax))
      {
 	switch (totalSpin)
 	  {
 	  case 2:
	    if ((totalIsospin == 1) && (totalEntanglement == 1))
	      ++Tmp;
 	    break;
 	  case 1:
	    switch (totalIsospin)
	      {
	      case 2:
		if (totalEntanglement == 1)
		  ++Tmp;
		break;
	      case 1:
		if (totalEntanglement != 1)
		  ++Tmp;
		break;
	      case 0:
		if (totalEntanglement == 1)
		  ++Tmp;
		break;
	      }
	    break;
 	  case 0:
	    if ((totalIsospin == 1)  && (totalEntanglement == 1))
	      ++Tmp;
	    break; 
 	  }
      }

  return  (Tmp + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin, totalEntanglement)
	   + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin - 1, totalEntanglement)
	   + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin - 1, totalEntanglement - 1)
	   + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin, totalEntanglement - 1)
	   + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, totalSpin, totalIsospin, totalEntanglement));

}

// evaluate Hilbert space dimension using previously generated Hilbert space dimension files (or compute them if they don't exist)
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

// long EvaluateHilbertSpaceDimensionWithDiskStorage(int nbrParticles, int nbrFluxQuanta, bool statistics,
// 						  long* lzDimensions, long* lDimensions, int lzMin, int lzMax)
// {
//   OutputFileName = new char[256];
//   if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
//     sprintf (OutputFileName, "bosons_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
//   else
//     sprintf (OutputFileName, "fermions_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
  
//   return FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + nbrFermions * lzMax) >> 1);
// }

// save dimensions in a given file
//
// outputFileName = output file name
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true for bosons, false for fermions
// lzDimensions = array taht contains dimension of each Lz sector
// lDimensions = array that contains dimension of each L sector (without taking into account the Lz degeneracy)
// lzMin = twice the minimum Lz value
// lzMax = twice the maximum Lz value
// totalDimension = total Hilbert space dimension

bool WriteDimensionToDisk(char* outputFileName, int nbrParticles, int nbrFluxQuanta, bool statistics,
			  long* lzDimensions, long* lDimensions, int lzMin, int lzMax, long totalDimension)
{
  ofstream File;
  File.open(outputFileName, ios::binary | ios::out);
  File << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " ";
  if (statistics == true)
    File << "bosons";
  else
    File << "femions";
  File << " on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  File << "# total Hilbert space dimension = " << totalDimension << endl << endl 
       << "N = " << nbrParticles << endl
       << "2S = " << nbrFluxQuanta << endl << endl
       << "#  dimensions for the Lz subspaces (starting from 2Lz = " << lzMin << " " << (lzMin + 2) << " ..." << endl
       << "Lz =";
  for (int x = lzMin; x <= lzMax; x += 2)
    File << " " << lzDimensions[(x - lzMin) >> 1];
  File << endl
       << "#  dimensions for the L subspaces (starting from 2L = " << lzMin << " " << (lzMin + 2) << " ..." << endl
       << "L =";
  for (int x = lzMin; x <= lzMax; x += 2)
    File << " " << lDimensions[(x - lzMin) >> 1];
  File << endl;
  File.close();
  return true;
}

// save dimensions in a given file (for particles with SU(2) spin)
//
// outputFileName = output file name
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true for bosons, false for fermions
// lzDimensions = array that contains dimension of each Lz sector per Sz
// lDimensions = array that contains dimension of each L sector per Sz (without taking into account the Lz degeneracy)
// lzMin = twice the minimum Lz value for each value of Sz
// lzMax = twice the maximum Lz value for each value of Sz
// sz = minal number of particles with spin up
// totalDimension = total Hilbert space dimension

bool SU2WriteDimensionToDisk(char* outputFileName, int nbrParticles, int nbrFluxQuanta, bool statistics,
			     long** lzDimensions, long** lDimensions, int* lzMin, int* lzMax, long totalDimension, int sz)
{
  ofstream File;
  File.open(outputFileName, ios::binary | ios::out);
  File << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " ";
  if (statistics == true)
    File << "bosons";
  else
    File << "femions";
  File << " with SU(2) spin on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  File << "# total Hilbert space dimension = " << totalDimension << endl << endl 
       << "N = " << nbrParticles << endl
       << "2S = " << nbrFluxQuanta << endl << endl
       << "#  dimensions for each subspaces with the following convention " << endl 
       << "(twice the total Lz value) (twice the total Sz value) (dimension of the subspace with fixed Lz and Sz) (dimension of the subspace with fixed L, Lz=L and Sz)" << endl;
  for (int MinSz = sz; MinSz <= nbrParticles; ++MinSz)
    for (int x = lzMin[MinSz]; x <= lzMax[MinSz]; x += 2)
      File << x << " " << ((2 * MinSz) - nbrParticles) << " " << lzDimensions[MinSz][(x - lzMin[MinSz]) >> 1] 
	   << " " << lDimensions[MinSz][(x - lzMin[MinSz]) >> 1] << endl;
  File << endl;
  File.close();
  return true;
}

// save dimensions in a given file (for particles with SU(4) spin)
//
// outputFileName = output file name
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true for bosons, false for fermions
// lzDimensions = array that contains dimension of each Lz sector per Sz and per Iz
// lDimensions = array that contains dimension of each L sector per Sz and per Iz (without taking into account the Lz degeneracy)
// lzMin = twice the minimum Lz value for each value of Sz and Iz
// lzMax = twice the maximum Lz value for each value of Sz and Iz
// sz = minal number of particles with spin up
// iz = minal number of particles with spin plus
// totalDimension = total Hilbert space dimension

bool SU4WriteDimensionToDisk(char* outputFileName, int nbrParticles, int nbrFluxQuanta, bool statistics,
			  long*** lzDimensions, long*** lDimensions, int** lzMin, int** lzMax, long totalDimension, int sz, int iz)
{
  sz = 0;
  iz = 0;
  ofstream File;
  File.open(outputFileName, ios::binary | ios::out);
  File << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " ";
  if (statistics == true)
    File << "bosons";
  else
    File << "femions";
  File << " with SU(4) spin on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  File << "# total Hilbert space dimension = " << totalDimension << endl << endl 
       << "N = " << nbrParticles << endl
       << "2S = " << nbrFluxQuanta << endl << endl
       << "#  dimensions for each subspaces with the following convention " << endl 
       << "(twice the total Lz value) (twice the total Sz value) (twice the total Iz value) (dimension of the subspace with fixed Lz, Sz and Iz) (dimension of the subspace with fixed L, Lz=L, Sz and Iz)" << endl;
  for (int MinSz = sz; MinSz < nbrParticles; ++MinSz)
    for (int MinIz = iz; MinIz < nbrParticles; ++MinIz)
      for (int x = lzMin[MinSz][MinIz]; x <= lzMax[MinSz][MinIz]; x += 2)
	File << x << " " << ((2 * MinSz) - nbrParticles)<< " " << ((2 * MinIz) - nbrParticles) << " " 
	     << lzDimensions[MinSz][MinIz][(x - lzMin[MinSz][MinIz]) >> 1] << " " << lDimensions[MinSz][MinIz][(x - lzMin[MinSz][MinIz]) >> 1] << endl;
  File << endl;
  File.close();
  return true;
}

