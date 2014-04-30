#include "config.h"

#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnTorusShort.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
// #include "Architecture/ArchitectureOperation/FQHESphereSymmetrizeU1U1StateOperation.h"

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
  
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file names");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("state-1") == 0)
    {
      cout << "error, an input file should be provided for the first component. See man page for option syntax or type FQHETorusMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-1")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-1") << endl;
    }
  if (Manager.GetString("state-2") == 0)
    {
      cout << "error, an input file should be provided for the second component. See man page for option syntax or type FQHETorusMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-2")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-2") << endl;
    }

  int NbrParticles1 = 0; 
  int NbrFluxQuanta1 = 0; 
  int TotalKy1 = 0;
  bool Statistics = true;
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("state-1"),
						   NbrParticles1, NbrFluxQuanta1, TotalKy1, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-1") << endl;
      return -1;
    }
  
  int NbrParticles2 = 0; 
  int NbrFluxQuanta2 = 0; 
  int TotalKy2 = 0;
  Statistics = true;
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("state-2"),
						   NbrParticles2, NbrFluxQuanta2, TotalKy2, Statistics) == false)
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
  

  RealVector State1;
  if (State1.ReadVector (Manager.GetString("state-1")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-1") << endl;
      return -1;      
    }
  RealVector State2;
  if (State2.ReadVector (Manager.GetString("state-2")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-2") << endl;
      return -1;      
    }

  BosonOnTorusShort* Space1 = 0;
  Space1 = new BosonOnTorusShort(NbrParticles1, NbrFluxQuanta1, TotalKy1);	       
  
  BosonOnTorusShort* Space2 = 0;
  Space2 = new BosonOnTorusShort(NbrParticles2, NbrFluxQuanta2, TotalKy2);	       
   
  bool UnnormalizedBasisFlag = false;
  
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = new char [512];
      sprintf (OutputFileName, "bosons_torus_kysym_symmetrized_n_%d_2s_%d_ky_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalKy1 + TotalKy2));
	
    }

  if (NbrFluxQuanta2 != NbrFluxQuanta1)
    {
      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of flux quanta" << endl;
      return -1;
    }

  BosonOnTorusShort* TargetSpace = 0;
  TargetSpace = new BosonOnTorusShort(NbrParticles1 + NbrParticles2, NbrFluxQuanta2, TotalKy1 + TotalKy2);	       
 
  RealVector OutputState = TargetSpace->SymmetrizeU1U1State (State1 , State2, Space1 , Space2 , UnnormalizedBasisFlag , Architecture.GetArchitecture());
    
	
  if (OutputState.WriteVector(OutputFileName) == false)
    {
      cout << "error while writing output state " << OutputFileName << endl;
      return -1;
    }
  delete Space1;
  delete Space2;
  delete TargetSpace;
}

