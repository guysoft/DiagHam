#include "config.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnTorusShort.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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
  (*SystemGroup) += new SingleStringOption  ('\n', "multiple-states", "text file that give the list of vectors to symmetrize");
  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new BooleanOption  ('s', "single-state", "vector file that corresponds to the second component");
  (*SystemGroup) += new BooleanOption  ('a', "sym-y", "apply antiperiodic conditions with respect to Ly before symmetrizing");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-orbitals", "number of orbitals to group together when using the single-state option", 2);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ky-momentum", "compute the vector with given ky in mode sym-y", 0);
//   (*SystemGroup) += new SingleBooleanOption  ('x', "magnetic-translation", "vector file that corresponds to the second component");
  
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

  if ((Manager.GetString("state-1") == 0) && (Manager.GetString("multiple-states") == 0))
    {
      cout << "error, an input file should be provided for the first component. See man page for option syntax or type FQHETorusMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if ((Manager.GetString("state-1") != 0) && (IsFile(Manager.GetString("state-1")) == false))
    {
      cout << "can't open file " << Manager.GetString("state-1") << endl;
    }
  

  bool HaveComplex = Manager.GetBoolean("complex");
  int NbrParticles1 = 0; 
  int NbrFluxQuanta1 = 0; 
  int TotalKy1 = 0;
  bool Statistics = true;
  BosonOnTorusShort* Space1 = 0;
  bool UnnormalizedBasisFlag = false;
  bool SingleStateFlag = Manager.GetBoolean("single-state");
  int NbrStates = 1;

  if (Manager.GetString("state-1") != 0)
    {  
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("state-1"),
						      NbrParticles1, NbrFluxQuanta1, TotalKy1, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("state-1") << endl;
	  return -1;
	}
      Space1 = new BosonOnTorusShort(NbrParticles1, NbrFluxQuanta1, TotalKy1);
    }
  else
    {
      MultiColumnASCIIFile MultipleStateFile;
      if (MultipleStateFile.Parse(Manager.GetString("multiple-states")) == false)
	{
	  MultipleStateFile.DumpErrors(cout);
	  return -1;
	}
      NbrStates = MultipleStateFile.GetNbrLines();
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(MultipleStateFile(0, 0),
						      NbrParticles1, NbrFluxQuanta1, TotalKy1, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << MultipleStateFile(0, 0) << endl;
	  return -1;
	}
      Space1 = new BosonOnTorusShort(NbrParticles1, NbrFluxQuanta1, TotalKy1);      
    }

  if ((SingleStateFlag == true) && ((NbrFluxQuanta1 % Manager.GetInteger("nbr-orbitals")) != 0))
    {
      cout << "error: number of flux quanta should be a multiple of " << Manager.GetInteger("nbr-orbitals") << " to perform symmetrization operation on single state"<< endl;
      return -1;
    }
	  
  
  BosonOnTorusShort** InputSpaces = 0;
  BosonOnTorusShort* TargetSpace = 0; 
  int TotalNbrParticles = NbrParticles1;
  int TotalKy = TotalKy1;
  if (SingleStateFlag == false)    
    {      
      if ((Manager.GetString("state-2") == 0) && (Manager.GetString("multiple-states") == 0))
	{
	  cout << "error, an input file should be provided for the second component. See man page for option syntax or type FQHETorusMultipleU1ToU1 -h" << endl;
	  return -1;
	}
      if (Manager.GetString("state-2") != 0)
	{
	  NbrStates = 2;
	  if (IsFile(Manager.GetString("state-2")) == false)
	    {
	      cout << "can't open file " << Manager.GetString("state-2") << endl;
	    }
	  InputSpaces = new BosonOnTorusShort*[NbrStates];
	  InputSpaces[0] = Space1;
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
	  if (NbrFluxQuanta2 != NbrFluxQuanta1)
	    {
	      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of flux quanta" << endl;
	      return -1;
	    }
	  TotalNbrParticles += NbrParticles2;
	  TotalKy += TotalKy2;
	  InputSpaces[1] = new BosonOnTorusShort(NbrParticles2, NbrFluxQuanta1, TotalKy2);
	} 
      else
	{
	  MultiColumnASCIIFile MultipleStateFile;
	  if (MultipleStateFile.Parse(Manager.GetString("multiple-states")) == false)
	    {
	      MultipleStateFile.DumpErrors(cout);
	      return -1;
	    }
	  InputSpaces = new BosonOnTorusShort*[NbrStates];
	  InputSpaces[0] = Space1;
	  for (int i = 1; i < NbrStates; ++i)
	    {
	      int NbrParticles2 = 0; 
	      int NbrFluxQuanta2 = 0; 
	      int TotalKy2 = 0;
	      Statistics = true;
	      if (FQHEOnTorusFindSystemInfoFromVectorFileName(MultipleStateFile(0, i),
							      NbrParticles2, NbrFluxQuanta2, TotalKy2, Statistics) == false)
		{
		  cout << "error while retrieving system parameters from file name " << Manager.GetString("state-2") << endl;
		  return -1;
		}
	      if (NbrFluxQuanta2 != NbrFluxQuanta1)
		{
		  cout << "error, " << Manager.GetString("state-1") << " and " << MultipleStateFile(0, i) << " don't have the same number of flux quanta" << endl;
		  return -1;
		}
	      TotalNbrParticles += NbrParticles2;
	      TotalKy += TotalKy2;
	      InputSpaces[i] = new BosonOnTorusShort(NbrParticles2, NbrFluxQuanta1, TotalKy2);
	    }
	}
      TotalKy %= NbrFluxQuanta1;
      cout << "target space N="<< TotalNbrParticles << " Nphi=" << NbrFluxQuanta1 << " Ky=" << TotalKy << endl;
      TargetSpace = new BosonOnTorusShort(TotalNbrParticles, NbrFluxQuanta1, TotalKy);	       
    }

  if (HaveComplex == false)
    {
      if (SingleStateFlag == false)
	{	      
	  RealVector* States = new RealVector[NbrStates];
	  if ((Manager.GetString("state-1") != 0) && (Manager.GetString("state-2") != 0))
	    {
	      if (States[0].ReadVector (Manager.GetString("state-1")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("state-1") << endl;
		  return -1;      
		}
	      if (States[1].ReadVector (Manager.GetString("state-2")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("state-2") << endl;
		  return -1;      
		}
	    }
	  else
	    {
	      MultiColumnASCIIFile MultipleStateFile;
	      if (MultipleStateFile.Parse(Manager.GetString("multiple-states")) == false)
		{
		  MultipleStateFile.DumpErrors(cout);
		  return -1;
		}
	      for (int i = 0; i < NbrStates; ++i)
		{
		  if (States[i].ReadVector (MultipleStateFile(0, i)) == false)
		    {
		      cout << "can't open vector file " << MultipleStateFile(0, i) << endl;
		      return -1;      
		    }
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
	      sprintf (OutputFileName, "bosons_torus_kysym_symmetrized_n_%d_2s_%d_ky_%d.0.vec", TotalNbrParticles, NbrFluxQuanta1, TotalKy);
	    }
	  
      	  RealVector OutputState = TargetSpace->SymmetrizeU1U1State (States, InputSpaces, NbrStates, Architecture.GetArchitecture());
	  if (OutputState.WriteVector(OutputFileName) == false)
	    {
	      cout << "error while writing output state " << OutputFileName << endl;
	      return -1;
	    }
	  
	  delete TargetSpace;
	}
      else
	{
	  RealVector State1;
	  if (State1.ReadVector (Manager.GetString("state-1")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("state-1") << endl;
	      return -1;      
	    }      
	  int NbrFluxQuanta = NbrFluxQuanta1 / Manager.GetInteger("nbr-orbitals");
	  char* OutputFileName = 0;
	  char* FullOutputFileName = 0;
	  RealVector* OutputStates = 0;
	  int* KySectors = 0;
	  int NbrKySectors = 0;
	  if (Manager.GetBoolean("sym-y") == false)
	    {
	      NbrKySectors = Space1->SymmetrizeSingleStateGroupingDistantOrbitals(State1, Manager.GetInteger("nbr-orbitals"), OutputStates, KySectors, Architecture.GetArchitecture()); 
	    }
	  else
	    {
	      NbrKySectors = Space1->SymmetrizeSingleStateGroupingNeighbouringOrbitals(State1, Manager.GetInteger("nbr-orbitals"), OutputStates, KySectors, Architecture.GetArchitecture()); 	      
	    }
	  int NbrGeneratedStates = 0;
	  if (Manager.GetString("output-file") != 0)
	    {
	      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
	      strcpy (OutputFileName, Manager.GetString("output-file"));
	    }
	  else
	    {
	      OutputFileName = new char [512];
	      if (Manager.GetBoolean("sym-y") == false)
		{
		  sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_xsymmetrized_n_%d_2s_%d", TotalKy1, NbrParticles1, NbrFluxQuanta);
		}
	      else
		{
		  sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_ysymmetrized_n_%d_2s_%d", TotalKy1, NbrParticles1, NbrFluxQuanta);
		}
	    }
	  for (int i = 0; i < NbrKySectors; ++i)
	    {
	      cout << "state generated in the Ky=" << KySectors[i] << " sector" << endl;
	      char* FullOutputFileName = new char [256 + strlen(OutputFileName)];
	      sprintf (FullOutputFileName , "%s_ky_%d.0.vec", OutputFileName, KySectors[i]);	      
	      if (OutputStates[i].WriteVector(FullOutputFileName) == false)
		{
		  cout << "error while writing output state " << FullOutputFileName << endl;
		  return -1;
		}
	      delete[] FullOutputFileName;
	    }
	  cout << "Symmetrization has generated " << NbrKySectors << " state(s)" << endl;
	}
    }
  else
    {
      if (SingleStateFlag == false)
	{	  
	  ComplexVector* States = new ComplexVector[NbrStates];
	  if ((Manager.GetString("state-1") != 0) && (Manager.GetString("state-2") != 0))
	    {
	      if (States[0].ReadVector (Manager.GetString("state-1")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("state-1") << endl;
		  return -1;      
		}
	      if (States[1].ReadVector (Manager.GetString("state-2")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("state-2") << endl;
		  return -1;      
		}
	    }
	  else
	    {
	      MultiColumnASCIIFile MultipleStateFile;
	      if (MultipleStateFile.Parse(Manager.GetString("multiple-states")) == false)
		{
		  MultipleStateFile.DumpErrors(cout);
		  return -1;
		}
	      for (int i = 0; i < NbrStates; ++i)
		{
		  if (States[i].ReadVector (MultipleStateFile(0, i)) == false)
		    {
		      cout << "can't open vector file " << MultipleStateFile(0, i) << endl;
		      return -1;      
		    }
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
	      sprintf (OutputFileName, "bosons_torus_kysym_symmetrized_n_%d_2s_%d_ky_%d.0.vec", TotalNbrParticles, NbrFluxQuanta1, TotalKy);
	    }
	  
	  	  
	  ComplexVector OutputState = TargetSpace->SymmetrizeU1U1State (States, InputSpaces, NbrStates, Architecture.GetArchitecture());

	  if (OutputState.WriteVector(OutputFileName) == false)
	    {
	      cout << "error while writing output state " << OutputFileName << endl;
	      return -1;
	    }
	  
	  delete TargetSpace;
	}
      else
	{
	  ComplexVector State1;
	  if (State1.ReadVector (Manager.GetString("state-1")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("state-1") << endl;
	      return -1;      
	    }      
	  int NbrFluxQuanta = NbrFluxQuanta1 / Manager.GetInteger("nbr-orbitals");
	  char* OutputFileName = 0;
	  char* FullOutputFileName = 0;
	  ComplexVector* OutputStates = 0;
	  int* KySectors = 0;
	  int NbrKySectors = 0;
	  if (Manager.GetBoolean("sym-y") == false)
	    {
	      NbrKySectors = Space1->SymmetrizeSingleStateGroupingDistantOrbitals(State1, Manager.GetInteger("nbr-orbitals"), OutputStates, KySectors, Architecture.GetArchitecture()); 
	    }
	  else
	    {
	      NbrKySectors = Space1->SymmetrizeSingleStateGroupingNeighbouringOrbitals(State1, Manager.GetInteger("nbr-orbitals"), OutputStates, KySectors, Architecture.GetArchitecture()); 	      
	    }
	  int NbrGeneratedStates = 0;
	  if (Manager.GetString("output-file") != 0)
	    {
	      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
	      strcpy (OutputFileName, Manager.GetString("output-file"));
	    }
	  else
	    {
	      OutputFileName = new char [512];
	      if (Manager.GetBoolean("sym-y") == false)
		{
		  sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_xsymmetrized_n_%d_2s_%d", TotalKy1, NbrParticles1, NbrFluxQuanta);
		}
	      else
		{
		  sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_ysymmetrized_n_%d_2s_%d", TotalKy1, NbrParticles1, NbrFluxQuanta);
		}
	    }
	  for (int i = 0; i < NbrKySectors; ++i)
	    {
	      cout << "state generated in the Ky=" << KySectors[i] << " sector" << endl;
	      char* FullOutputFileName = new char [256 + strlen(OutputFileName)];
	      sprintf (FullOutputFileName , "%s_ky_%d.0.vec", OutputFileName, KySectors[i]);	      
	      if (OutputStates[i].WriteVector(FullOutputFileName) == false)
		{
		  cout << "error while writing output state " << FullOutputFileName << endl;
		  return -1;
		}
	      delete[] FullOutputFileName;
	    }
	  cout << "Symmetrization has generated " << NbrKySectors << " state(s)" << endl;
	}
    }
  
  delete Space1;
}

