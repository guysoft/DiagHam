#include "config.h"

#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHESphereSymmetrizeU1U1StateOperation.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereMultipleU1ToU1" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
	
  ArchitectureManager Architecture;
	
  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleStringOption  ('1', "state-1", "vector file that corresponds to the first component");
  (*SystemGroup) += new SingleStringOption  ('2', "state-2", "vector file that corresponds to the second component");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane-1", "use squeezed basis instead of the usual n-body basis for the first component");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane-2", "use squeezed basis instead of the usual n-body basis for the second component");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane-output", "use squeezed basis instead of the usual n-body basis for the symmetrized state");
  (*SystemGroup) += new BooleanOption  ('\n', "single-state", "symmetrize a unique state with even number of orbitals");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file1", "use a file as the definition of the reference state for the first component squeezed basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file2", "use a file as the definition of the reference state for the second component squeezed basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-fileoutput", "use a file as the definition of the reference state for the symmetrized state squeezed basis");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file names");
  (*OutputGroup) += new BooleanOption  ('u', "unnormalized-basis", "indicates that calculations and data are in the unnormalized basis");
  (*SystemGroup) += new BooleanOption  ('\n', "rational", "use rational numbers instead of double precision floating point numbers");
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the target system (in single-state mode)", 0);
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("state-1") == 0)
    {
      cout << "error, an input file should be provided for the first component. See man page for option syntax or type FQHESphereMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-1")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-1") << endl;
    }
  if (Manager.GetString("state-2") == 0 && Manager.GetBoolean("single-state") == false)
    {
      cout << "error, an input file should be provided for the second component. See man page for option syntax or type FQHESphereMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("single-state") == false && IsFile(Manager.GetString("state-2")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-2") << endl;
    }

  int NbrParticles1 = 0; 
  int NbrFluxQuanta1 = 0; 
  int TotalLz1 = 0;
  bool Statistics = true;
  
  int NbrParticles2 = 0; 
  int NbrFluxQuanta2 = 0; 
  int TotalLz2 = 0;
  
    
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("state-1"),
						   NbrParticles1, NbrFluxQuanta1, TotalLz1, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-1") << endl;
      return -1;
    }
  
  RealVector State1;
  LongRationalVector RationalState1;
  RealVector State2;
  if (Manager.GetBoolean("rational") == false)
  {
    if (State1.ReadVector (Manager.GetString("state-1")) == false)
      {
	cout << "can't open vector file " << Manager.GetString("state-1") << endl;
	return -1;      
      }
  }
  else
  {
    if (RationalState1.ReadVector (Manager.GetString("state-1")) == false)
      {
	cout << "can't open vector file " << Manager.GetString("state-1") << endl;
	return -1;      
      }
  }
  

  BosonOnSphereShort* Space1 = 0;
  BosonOnSphereShort* Space2 = 0;
  if (Manager.GetBoolean("haldane-1") == true)
    {
      int* ReferenceState = 0;
      if (FQHEGetRootPartition(Manager.GetString("reference-file1"), NbrParticles1, NbrFluxQuanta1, ReferenceState) == false)
	return -1;
      Space1 = new BosonOnSphereHaldaneBasisShort(NbrParticles1, TotalLz1, NbrFluxQuanta1, ReferenceState);	       
    }
  else
    {
      Space1 = new BosonOnSphereShort(NbrParticles1, TotalLz1, NbrFluxQuanta1);	       
    }
  
  if (Manager.GetBoolean("single-state") == false)
  {
    Statistics = true;
    if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("state-2"),
						   NbrParticles2, NbrFluxQuanta2, TotalLz2, Statistics) == false)
      {
	cout << "error while retrieving system parameters from file name " << Manager.GetString("state-2") << endl;
	return -1;
      }
  /*  if (NbrParticles1 != NbrParticles2)
      {
      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of particles" << endl;
      return -1;
      }*/
  if (NbrFluxQuanta2 != NbrFluxQuanta1)
    {
      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of flux quanta" << endl;
      return -1;
    }
    
    
  if (State2.ReadVector (Manager.GetString("state-2")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-2") << endl;
      return -1;      
    }
  
  if (Manager.GetBoolean("haldane-2") == true)
    {
      int* ReferenceState = 0;
      if (FQHEGetRootPartition(Manager.GetString("reference-file2"), NbrParticles2, NbrFluxQuanta2, ReferenceState) == false)
	return -1;
      Space2 = new BosonOnSphereHaldaneBasisShort(NbrParticles2, TotalLz2, NbrFluxQuanta2, ReferenceState);	       
    }
  else
    {
      Space2 = new BosonOnSphereShort(NbrParticles2, TotalLz2, NbrFluxQuanta2);	       
    }
  
  }
  else
  {
   if ((NbrFluxQuanta1 % 2) == 0) 
   {
    cout << "Error: number of orbitals should be even for symmetrization of a single state" << endl;
    return -1;
   }
  }
  
  
  
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = new char [512];
      if (Manager.GetBoolean("single-state") == false)
      {
	if (Manager.GetBoolean("unnormalized-basis") == false)
	  {
	    if (Manager.GetBoolean("haldane-output") == true)
	      sprintf (OutputFileName, "bosons_haldane_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
	    else
	      sprintf (OutputFileName, "bosons_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
	  }
	else
	  {
	    if (Manager.GetBoolean("haldane-output") == true)
	      sprintf (OutputFileName, "bosons_haldane_unnormalized_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
	    else
	      sprintf (OutputFileName, "bosons_unnormalized_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
	  }
      }
      else
      {
	if (Manager.GetBoolean("unnormalized-basis") == false)
	 {
	   cout << "Normalized basis not implemented" << endl;
	   return -1;
	 }
	else
	{
	  sprintf(OutputFileName, "bosons_unnormalized_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1, (NbrFluxQuanta1 - 1)/2, Manager.GetInteger("total-lz"));
	}
      }
    }

/*  if (NbrParticles1 != NbrParticles2)
    {
      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of particles" << endl;
      return -1;
    }*/
  
  BosonOnSphereShort* TargetSpace = 0;
  if (Manager.GetBoolean("single-state") == false)
  {
    if (Manager.GetBoolean("haldane-output") == true)
      {
	int* ReferenceState = 0;
	int NbrParticles = 0;
	int NbrFluxQuanta = 0;
	if (FQHEGetRootPartition(Manager.GetString("reference-fileoutput"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
	  return -1;
	if (NbrParticles != (NbrParticles1 + NbrParticles2))
	  {
	    cout << "error : number of particles for the output state ( " << NbrParticles << " ) is different from the one computed from the input states (" << (NbrParticles1 + NbrParticles2) << ")" << endl;
	    return -1;
	  }
	if (NbrFluxQuanta != NbrFluxQuanta2)
	  {
	    cout << "error : number of flux quanta for the output state ( " << NbrFluxQuanta << " ) is different from the one computed from the input states (" << NbrFluxQuanta2 << ")" << endl;
	    return -1;
	  }
	int TotalLz = 0;
	TargetSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, NbrFluxQuanta2, ReferenceState);	       
	if (TotalLz != (TotalLz1 + TotalLz2))
	  {
	    cout << "error : total Lz for the output state ( " << TotalLz << " ) is different from the one computed from the input states (" << (TotalLz1 + TotalLz2) << ")" << endl;
	    return -1;
	  }
      }
    else
    {
      TargetSpace = new BosonOnSphereShort(NbrParticles1 + NbrParticles2, TotalLz1 + TotalLz2, NbrFluxQuanta2);	       
    }
    RealVector OutputState = TargetSpace->SymmetrizeU1U1State (State1 , State2, Space1 , Space2 , Manager.GetBoolean("unnormalized-basis") , Architecture.GetArchitecture());
    if (OutputState.WriteVector(OutputFileName) == false)
    {
      cout << "error while writing output state " << OutputFileName << endl;
      return -1;
    }
  }
  else
  {
   int NbrFluxQuanta = (NbrFluxQuanta1 - 1)/2;
   int TotalLz = Manager.GetInteger("total-lz");
   TargetSpace = new BosonOnSphereShort(NbrParticles1, TotalLz, NbrFluxQuanta);
   RealVector OutputState;
   LongRationalVector RationalOutputState;
   if (Manager.GetBoolean("rational") == false)
   {
    OutputState = TargetSpace->SymmetrizeU1U1SingleState(State1, Space1, Manager.GetBoolean("unnormalized-basis"));
    if (OutputState.Norm() == 0)
    {
      cout << "Symmetrized state is zero. No output." << endl;
      return -1;
    }
    else
    {
      if (OutputState.WriteVector(OutputFileName) == false)
      {
	cout << "error while writing output state " << OutputFileName << endl;
	return -1;
      }
      if (Manager.GetBoolean("unnormalized-basis") == true);
      {
	int RootPosition = 0;
	while (abs(OutputState[RootPosition]) < 1.0e-12)
	{
	  RootPosition += 1;
	}
	double RootCoef = OutputState[RootPosition];
	for (int TmpPos = 0; TmpPos < TargetSpace->GetHilbertSpaceDimension(); ++TmpPos)
	{
	 OutputState[TmpPos] /= RootCoef; 
	}
      }
    }
   }
   else
   {
     RationalOutputState = TargetSpace->SymmetrizeU1U1SingleState(RationalState1, Space1, Manager.GetBoolean("unnormalized-basis"));
     bool zeroFlag = true;
     for (int TmpPos = 0; TmpPos < TargetSpace->GetHilbertSpaceDimension(); TmpPos ++)
      {
	if (RationalOutputState[TmpPos].IsZero() == false)
	  zeroFlag = false;
      }
      if (zeroFlag)
      {
	cout << "Symmetrized state is zero. No output." << endl;
	return -1;
      }
      
	if (RationalOutputState.WriteVector(OutputFileName) == false)
	{
	  cout << "error while writing output state " << OutputFileName << endl;
	  return -1;
	}
	if (Manager.GetBoolean("unnormalized-basis") == true);
	{
	  int RootPosition = 0;
	  while (RationalOutputState[RootPosition].IsZero())
	  {
	    RootPosition += 1;
	  }
	  LongRational RootCoef = RationalOutputState[RootPosition];
	  RationalOutputState /= RootCoef;
	} 
      }
    }
  
  delete Space1;
  delete Space2;
  delete TargetSpace;
}

