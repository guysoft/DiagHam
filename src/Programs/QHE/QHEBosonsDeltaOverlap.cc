#include "Vector/RealVector.h"

#include "HilbertSpace/QHEHilbertSpace/BosonOnSphere.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereDeltaHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereDeltaModifiedHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereCoulombDeltaHamiltonian.h"
#include "FunctionBasis/QHEFunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MainTask/QHEMainTask/QHEOnSphereMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


Complex LaughlinWaveFunction(RealVector& position, int nbrBosons);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEBosonsDelta" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleStringOption  ('\n', "exact-state", "name of the file containing the vector obtained using exact diagonalization");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsDeltaOverlap -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();

  if (((SingleStringOption*) Manager["exact-state"])->GetString() == 0)
    {
      cout << "QHEBosonsDeltaOverlap requires an exact state" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["exact-state"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["exact-state"])->GetString() << endl;
      return -1;      
    }
  BosonOnSphere Space (NbrBosons, 0, LzMax);
  RealVector Location(2 * NbrBosons, true);
  srand48(29457);
  for (int k = 0; k < 100; ++k)
    {
      for (int i = 0; i < NbrBosons; ++i)
	{
	  Location[i << 1] = M_PI * drand48();
	  Location[(i << 1) + 1] = 2.0 * M_PI * drand48();
	  cout << Location[i << 1] << " " << Location[(i << 1) + 1] << endl;
	}
      //  Location[4] = Location[0];
      //  Location[5] = Location[1];
      ParticleOnSphereFunctionBasis Basis(LzMax);
      QHEParticleWaveFunctionOperation Operation(&Space, &State, &Location, &Basis);
      Architecture.GetArchitecture()->ExecuteOperation(&Operation);      
      Complex ValueExact (Operation.GetScalar());
      //      Complex ValueExact = Space.EvaluateWaveFunction(State, Location, Basis);
      Complex ValueLaughlin = LaughlinWaveFunction(Location, NbrBosons) * 0.36563112422012;
      cout << ValueExact  << endl; 
      cout << "-------------------------------------" << endl;
//      cout << ValueExact  << " " << ValueLaughlin << " " << (Norm(ValueExact) / Norm(ValueLaughlin)) << endl;  
      
    }
  return 0;
}


Complex LaughlinWaveFunction(RealVector& position, int nbrBosons)
{
  Complex Value (1.0, 0.0);
  Complex Tmp;
  for (int i = 0; i < (nbrBosons - 1); ++i)
    for (int j = i + 1; j < nbrBosons; ++j)
      {
	Tmp.Re = sin(0.5 * (position[i << 1] - position[j << 1])) * cos(0.5 * (position[1 + (i << 1)] - position[1 + (j << 1)]));
	Tmp.Im = sin(0.5 * (position[i << 1] + position[j << 1])) * sin(0.5 * (position[1 + (i << 1)] - position[1 + (j << 1)]));
	Value *= Tmp;
	Value *= Tmp;
      }
  return Value;
}
