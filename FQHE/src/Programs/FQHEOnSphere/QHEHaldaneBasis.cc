#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  
  // some running options and help
  OptionManager Manager ("QHEHaldaneSphere" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  
  ArchitectureManager Architecture;
  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be taken from reference state file definition)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be taken from reference state file definition)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the lz value corresponding to the eigenvector", 0);
  (*SystemGroup) += new SingleStringOption  ('s', "state", "name of the file containing the eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n', "huge", "use huge Hilbert space support");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "file-size", "maximum file size (in MBytes) when using huge mode", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldabe algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEHaldaneSphere -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["lz-value"])->GetInteger();

  if (((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "QHEHaldaneSphere requires a state" << endl;
      return -1;
    }
//   RealVector State;
//   if (State.ReadVector (((SingleStringOption*) Manager["state"])->GetString()) == false)
//     {
//       cout << "can't open vector file " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
//       return -1;      
//     }


//   ParticleOnSphere* Space;
// #ifdef __64_BITS__
//   if (LzMax <= 63)
//     {
//       Space = new FermionOnSphere(NbrParticles, Lz, LzMax);
//     }
//   else
//     {
//       Space = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax);
//     }
// #else
//   if (LzMax <= 31)
//     {
//       Space = new FermionOnSphere(NbrParticles, Lz, LzMax);
//     }
//   else
//     {
//       Space = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax);
//     }
// #endif


//   int NbrNonZeroComponents = 0;
//   for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
//     {
//       if (fabs(State[i]) > 1e-12)
// 	{
// // 	  Space->PrintState(cout, i) << " : " << (State[i]) << endl;
// 	  ++NbrNonZeroComponents;
// 	}
//     }
//   cout << NbrNonZeroComponents << " / " << Space->GetHilbertSpaceDimension() << endl;

  int* ReferenceState = 0;
  if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
    {
      ReferenceState = new int[LzMax + 1];
      for (int i = 0; i <= LzMax; ++i)
	ReferenceState[i] = 0;
      if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "laughlin") == 0)
	for (int i = 0; i <= LzMax; i += 3)
	  ReferenceState[i] = 1;
      else
	if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "pfaffian") == 0)
	  for (int i = 0; i <= LzMax; i += 4)
	    {
	      ReferenceState[i] = 1;
	      ReferenceState[i + 1] = 1;
	    }
	else
	  if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "readrezayi3") == 0)
	    for (int i = 0; i <= LzMax; i += 5)
	      {
		ReferenceState[i] = 1;
		ReferenceState[i + 1] = 1;
		ReferenceState[i + 2] = 1;
	      }
	  else
	    {
	      cout << "unknown reference state " << ((SingleStringOption*) Manager["reference-state"])->GetString() << endl;
	      return -1;
	    }
    }
  else
    {
      ConfigurationParser ReferenceStateDefinition;
      if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
	{
	  ReferenceStateDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	{
	  cout << "NbrParticles is not defined or as a wrong value" << endl;
	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
	{
	  cout << "LzMax is not defined or as a wrong value" << endl;
	  return -1;
	}
      int MaxNbrLz;
      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	{
	  cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
	  return -1;     
	}
      if (MaxNbrLz != (LzMax + 1))
	{
	  cout << "wrong LzMax value in ReferenceState" << endl;
	  return -1;     
	}
    }
  if (((BooleanOption*) Manager["huge"])->GetBoolean() == true)
    {
      FermionOnSphereHaldaneHugeBasis ReducedSpace(NbrParticles, Lz, LzMax, ((SingleIntegerOption*) Manager["file-size"])->GetInteger(),
						   ReferenceState, ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20);
      cout << ReducedSpace.GetHilbertSpaceDimension() << endl;
    }
  else
    {
      FermionOnSphereHaldaneBasis ReducedSpace(NbrParticles, Lz, LzMax, ReferenceState);
      //  for (int i = 0; i < ReducedSpace.GetHilbertSpaceDimension(); ++i)
      //    ReducedSpace.PrintState(cout, i) << endl;
      cout << ReducedSpace.GetHilbertSpaceDimension() << endl;
    }
  return 0;
}


