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
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the lz value corresponding to the eigenvector", 0);
  (*SystemGroup) += new SingleStringOption  ('s', "state", "name of the file containing the eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n', "huge", "use huge Hilbert space support");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "file-size", "maximum file size (in MBytes) when using huge mode", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  
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

  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
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
//       Space = new FermionOnSphere(NbrFermions, Lz, LzMax);
//     }
//   else
//     {
//       Space = new FermionOnSphereUnlimited(NbrFermions, Lz, LzMax);
//     }
// #else
//   if (LzMax <= 31)
//     {
//       Space = new FermionOnSphere(NbrFermions, Lz, LzMax);
//     }
//   else
//     {
//       Space = new FermionOnSphereUnlimited(NbrFermions, Lz, LzMax);
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

  if (((BooleanOption*) Manager["huge"])->GetBoolean() == true)
    {
      FermionOnSphereHaldaneHugeBasis ReducedSpace(NbrFermions, Lz, LzMax, ((SingleIntegerOption*) Manager["file-size"])->GetInteger(),
						   ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20);
      cout << ReducedSpace.GetHilbertSpaceDimension() << endl;
    }
  else
    {
      FermionOnSphereHaldaneBasis ReducedSpace(NbrFermions, Lz, LzMax);
      //  for (int i = 0; i < ReducedSpace.GetHilbertSpaceDimension(); ++i)
      //    ReducedSpace.PrintState(cout, i) << endl;
      cout << ReducedSpace.GetHilbertSpaceDimension() << endl;
    }
  return 0;
}


