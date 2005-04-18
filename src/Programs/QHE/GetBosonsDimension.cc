#include "HilbertSpace/QHEHilbertSpace/BosonOnSphere.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>

using std::cout;
using std::endl;

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

// fake run to generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// currentLzMax = momentum maximum value for bosons that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// memory = reference on amount of memory needed
// return value = position from which new states have to be stored
int FakeGenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int pos, int& memory);


int main(int argc, char** argv)
{
  OptionManager Manager ("GetBosonsDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-particles", "number of particles", 20);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-lz", "number of particles", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-lz", "number of particles", 20);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "l-dimension", "get dimensions of all subspaces with fixed total l value");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((BooleanOption*) Manager["l-dimension"])->GetBoolean() == false)
    {
      for (int NbrBosons = ((SingleIntegerOption*) Manager["min-particles"])->GetInteger(); 
	   NbrBosons <= ((SingleIntegerOption*) Manager["max-particles"])->GetInteger(); ++NbrBosons)
	{
	  for (int LzMax = ((SingleIntegerOption*) Manager["min-lz"])->GetInteger(); 
	       LzMax <= ((SingleIntegerOption*) Manager["max-lz"])->GetInteger(); ++LzMax)
	    {
	      if (((BooleanOption*) Manager["fermion"])->GetBoolean() == true)
		{
		  int Max = ((LzMax - NbrBosons + 1) * NbrBosons);
		  int  L = 0;
		  if ((abs(Max) & 1) != 0)
		    L = 1;
		  cout << FermionEvaluateHilbertSpaceDimension(NbrBosons, LzMax, L) << " ";
		}
	      else
		{
		  int Max = (LzMax * NbrBosons);
		  int  L = 0;
		  if ((abs(Max) & 1) != 0)
		    L = 1;
		  cout << EvaluateHilbertSpaceDimension(NbrBosons, LzMax, L) << " ";
		}
	    }
	  cout << endl;
	}
    }
  else
    {
     int NbrBosons = ((SingleIntegerOption*) Manager["min-particles"])->GetInteger();
     int LzMax = ((SingleIntegerOption*) Manager["min-lz"])->GetInteger();
     int TotalLzMax = 0;
     int StartL = 0;
     if (((BooleanOption*) Manager["fermion"])->GetBoolean() == true)
       {
	 TotalLzMax = ((LzMax - NbrBosons + 1) * NbrBosons);
       }
     else
       {
	 TotalLzMax = (LzMax * NbrBosons);
       }
     StartL = 0;
     if ((abs(TotalLzMax) & 1) != 0)
       StartL = 1;
     long* LzDimensions = new long [1 + ((TotalLzMax - StartL) >> 1)];
     for (int L = StartL; L <= TotalLzMax; L += 2)
       {
	 long TmpDimension = 0;
	 if (((BooleanOption*) Manager["fermion"])->GetBoolean() == true)
	   {
	     TmpDimension = FermionEvaluateHilbertSpaceDimension(NbrBosons, LzMax, L);
	   }
	 else
	   {
	     TmpDimension = EvaluateHilbertSpaceDimension(NbrBosons, LzMax, L);
	   }
	 LzDimensions[(L - StartL) >> 1] = TmpDimension;
       }
     cout << "N = " << NbrBosons << endl;
     cout << "2S = " << LzMax << endl;
     cout << "Lz = ";
     for (int L = StartL; L < TotalLzMax; L += 2)
       {
	 cout << LzDimensions[(L - StartL) >> 1] << " ";
       }
     cout << LzDimensions[( TotalLzMax- StartL) >> 1] << endl;
     cout << "L = ";
     for (int L = StartL; L < TotalLzMax; L += 2)
       {
	 cout << (LzDimensions[(L - StartL) >> 1] - LzDimensions[1 + ((L - StartL) >> 1)]) << " ";
       }
     cout << LzDimensions[( TotalLzMax- StartL) >> 1] << endl;    
     delete[] LzDimensions;
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
    return 0;
  if (((nbrBosons * lzMax) == totalLz) || (lzMax == 0) || (totalLz == 0))
    {
      return 1;
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

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// currentLzMax = momentum maximum value for bosons that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// memory = reference on amount of memory needed
// return value = position from which new states have to be stored

int FakeGenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int pos, int& memory)
{
  if ((nbrBosons == 0) || ((nbrBosons * currentLzMax) < totalLz))
    {
      return pos;
    }
  if ((nbrBosons * currentLzMax) == totalLz)
    {
      if (memory > 0)
	memory += sizeof(int) * (lzMax + 1);
      return pos + 1;
    }
  if ((currentLzMax == 0) || (totalLz == 0))
    {
      memory += sizeof(int) * (lzMax + 1);
      return pos + 1;
    }

  int TmpTotalLz = totalLz / currentLzMax;
  int TmpNbrBosons = nbrBosons - TmpTotalLz;
  TmpTotalLz = totalLz - TmpTotalLz * currentLzMax;
  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      TmpPos = FakeGenerateStates(TmpNbrBosons, lzMax, ReducedCurrentLzMax, TmpTotalLz, pos, memory);
      ++TmpNbrBosons;
      pos = TmpPos;
      TmpTotalLz += currentLzMax;
    }
  if (lzMax == currentLzMax)
    return FakeGenerateStates(nbrBosons, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, pos, memory);
  else
    return FakeGenerateStates(nbrBosons, lzMax, ReducedCurrentLzMax, totalLz, pos, memory);
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
