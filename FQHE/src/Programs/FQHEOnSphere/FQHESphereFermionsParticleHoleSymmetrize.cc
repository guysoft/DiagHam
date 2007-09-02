#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereFermionsParticleHoleSymmetrize.cc" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "input state file name");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new BooleanOption  ('r', "symmetrize", "symmetrize state (instead of unsymmetrizing it)");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (changing N and 2S)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsParticleHoleSymmetrize.cc -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((SingleStringOption*) Manager["input-file"])->GetString() == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereFermionsParticleHoleSymmetrize -h" << endl;
      return -1;
    }
  if (IsFile(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrFluxQuanta = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger(); 
  int NbrHoles = NbrFluxQuanta + 1 - NbrParticles;

  bool SymmetrizeFlag = ((BooleanOption*) Manager["symmetrize"])->GetBoolean();
  bool Statistics = true;
  int TotalLz = 0;
  if (QHEOnSphereFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["input-file"])->GetString(),
						  NbrParticles, NbrFluxQuanta, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;
    }
  if ((((BooleanOption*) Manager["boson"])->GetBoolean() == true) || (((BooleanOption*) Manager["fermion"])->GetBoolean() == true))
    {
      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	Statistics = false;
      else
	Statistics = true;
    }
  if (((NbrParticles * NbrFluxQuanta) & 1) != (TotalLz != 0))
    {
      cout << "incompatible values for nbr-particles, nbr-flux and total-lz" << endl;
      return -1;
    }

  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;      
    }


  if (Statistics == true)
    {
      if (((BooleanOption*) Manager["haldane"])->GetBoolean() == false)
	{
	  RealVector OutputState;
	  FermionOnSphereSymmetricBasis InitialSpace(NbrParticles, NbrFluxQuanta);
	  FermionOnSphere TargetSpace(NbrParticles, TotalLz, NbrFluxQuanta);
	  if (SymmetrizeFlag)
	    {
	      if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;
 		}
	      OutputState = InitialSpace.ConvertToSymmetricNbodyBasis(State, TargetSpace);
	    }
	  else
	    {
	      if (InitialSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;
 		}
	      OutputState = InitialSpace.ConvertToNbodyBasis(State, TargetSpace);
	    }
	  if (OutputState.WriteVector(((SingleStringOption*) Manager["output-file"])->GetString()) == false)
	    {
	      cout << "error while writing output state " << ((SingleStringOption*) Manager["output-file"])->GetString() << endl;
	      return -1;
	    }
	}
      else
	{
	  int* ReferenceState = 0;
	  if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
	    {
	      ReferenceState = new int[NbrFluxQuanta + 1];
	      for (int i = 0; i <= NbrFluxQuanta; ++i)
		ReferenceState[i] = 0;
	      if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "laughlin") == 0)
		for (int i = 0; i <= NbrFluxQuanta; i += 3)
		  ReferenceState[i] = 1;
	      else
		if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "pfaffian") == 0)
		  for (int i = 0; i <= NbrFluxQuanta; i += 4)
		    {
		      ReferenceState[i] = 1;
		      ReferenceState[i + 1] = 1;
		    }
		else
		  if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "readrezayi3") == 0)
		    for (int i = 0; i <= NbrFluxQuanta; i += 5)
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
	      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", NbrFluxQuanta) == false) || (NbrFluxQuanta <= 0))
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
	      if (MaxNbrLz != (NbrFluxQuanta + 1))
		{
		  cout << "wrong LzMax value in ReferenceState" << endl;
		  return -1;     
		}
	    }
	  RealVector OutputState;
	  FermionOnSphereHaldaneSymmetricBasis InitialSpace(NbrParticles, NbrFluxQuanta, ReferenceState);	  
	  FermionOnSphereHaldaneBasis TargetSpace(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);

	  if (SymmetrizeFlag)
	    {
	      if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;
 		}
	      OutputState = InitialSpace.ConvertToSymmetricHaldaneNbodyBasis(State, TargetSpace);
	    }
	  else
	    {
	      if (InitialSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		{
		  cout << "dimension mismatch between Hilbert space and input state" << endl;
		  return -1;
 		}
	      OutputState = InitialSpace.ConvertToHaldaneNbodyBasis(State, TargetSpace);
	    }
 	  if (OutputState.WriteVector(((SingleStringOption*) Manager["output-file"])->GetString()) == false)
 	    {
 	      cout << "error while writing output state " << ((SingleStringOption*) Manager["output-file"])->GetString() << endl;
 	      return -1;
 	    }
	}
    }
}

