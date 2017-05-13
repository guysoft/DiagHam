#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"
 
#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "MathTools/BinomialCoefficients.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereBosonicStateTimesFermionicStateOperation.h"


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
  OptionManager Manager ("FQHESphereLLLProjectedLLLTimesLLL" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  ArchitectureManager Architecture;
  Architecture.AddOptionGroup(&Manager);
  
  (*SystemGroup) += new SingleStringOption  ('1', "state-1", "vector file that corresponds to the first state (should be bosonic)");
  (*SystemGroup) += new SingleStringOption  ('2', "state-2", "vector file that corresponds to the second state");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file1", "use a squeezed basis for the first state described by a given root configuration");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file2", "use a squeezed basis for the first state described by a given root configuration");
  (*SystemGroup) += new BooleanOption ('\n', "reverse-flux", "the fluxes bind to each particle are in the opposite direction than the magnetic field");
  (*SystemGroup) += new BooleanOption ('\n',"disable-lzsymmetry","disable the Lz<->-Lz symmetry for the Lz=0 sector");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "normalize the projected state assuming the sphere geometry");
  (*OutputGroup) += new BooleanOption  ('\n', "rational" , "use rational numbers instead of double precision floating point numbers");
  (*OutputGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*OutputGroup) += new SingleIntegerOption ('\n', "outputvector-index", "set the index of the output vector (i.e. the integer in the extention *.xxx.vec)", 0);
  (*OutputGroup) += new BooleanOption  ('\n', "save-partial" , "save the partial results");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereLLLProjectedLLLTimesLLL -h" << endl;
      return -1;
    }	
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrParticles1 = 0; 
  int NbrFluxQuanta1 = 0; 
  int TotalLz1 = 0;
  bool Statistics1 = true;
  
  int NbrParticles2 = 0; 
  int NbrFluxQuanta2 = 0; 
  int TotalLz2 = 0;
  bool Statistics2 = true;

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("state-1"),
						   NbrParticles1, NbrFluxQuanta1, TotalLz1, Statistics1) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-1") << endl;
      return -1;
    }
  if (Statistics1 == true)
    {
      cout << Manager.GetString("state-1") << " should be bosonic" << endl;
      return -1;
    }
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("state-2"),
						   NbrParticles2, NbrFluxQuanta2, TotalLz2, Statistics2) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-2") << endl;
      return -1;
    }

  if (NbrParticles1 != NbrParticles2)
    {
      cout << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " have a different number of particles" << endl;
      return -1;
    }

  int NbrFluxQuantumOutputState = 0;
  int TotalLzOutputState = 0;
  if (Manager.GetBoolean("reverse-flux") == true)
    {
      NbrFluxQuantumOutputState = NbrFluxQuanta2 - NbrFluxQuanta1; 
      TotalLzOutputState = TotalLz2 - TotalLz1;
      if (NbrFluxQuantumOutputState < 0)
	{
	  cout << "the number of flux quanta of " << Manager.GetString("state-1") << " should be lower than " << Manager.GetString("state-2") << endl;
	  return -1;
	}      
    }
  else
    {
      NbrFluxQuantumOutputState = NbrFluxQuanta2 + NbrFluxQuanta1; 
      TotalLzOutputState = TotalLz2 + TotalLz1;
    }

  bool LzSymmmetryFlag = true;
  if ((TotalLzOutputState != 0) || (Manager.GetBoolean("disable-lzsymmetry") == true))
    {
      LzSymmmetryFlag = false;
    }
  char* DiscreteSymmetryName = new char[128];
  if (LzSymmmetryFlag == false)
    {
      sprintf (DiscreteSymmetryName, "");
    }
  else
    {
      if (Manager.GetBoolean("minus-lzparity") == false)
	{
	  sprintf (DiscreteSymmetryName, "_lzsym_1");
	}
      else
	{
	  sprintf (DiscreteSymmetryName, "_lzsym_-1");
	}
    }

  RealVector BosonicInputVector;
  if (BosonicInputVector.ReadVector (Manager.GetString("state-1")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-1") << endl;
      return -1;      
    }
  BosonOnSphereShort* BosonicInputSpace = 0;
  if (Manager.GetString("reference-file1") == 0)
    {
      BosonicInputSpace = new BosonOnSphereShort(NbrParticles1, TotalLz1, NbrFluxQuanta1);  
    }
  else
    {
      int* ReferenceState = 0;
      if (FQHEGetRootPartition(Manager.GetString("reference-file1"), NbrParticles1, NbrFluxQuanta1, ReferenceState) == false)
	return -1;
      BosonicInputSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles1, TotalLz1, NbrFluxQuanta1, ReferenceState);        
    }

  RealVector InputVector;
  if (InputVector.ReadVector (Manager.GetString("state-2")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-2") << endl;
      return -1;      
    }
  ParticleOnSphere* InputSpace = 0;  
  if (Manager.GetString("reference-file2") == 0)
    {
      if (Statistics2 == true)
	{
	  InputSpace = new FermionOnSphere(NbrParticles2, TotalLz2, NbrFluxQuanta2);  
	}
      else
	{
	  InputSpace = new BosonOnSphereShort(NbrParticles2, TotalLz2, NbrFluxQuanta2);  
	}
    }
  else
    {
      int* ReferenceState = 0;
      if (FQHEGetRootPartition(Manager.GetString("reference-file2"), NbrParticles2, NbrFluxQuanta2, ReferenceState) == false)
	return -1;
      if (Statistics2 == true)
	{
	  InputSpace = new FermionOnSphereHaldaneBasis(NbrParticles2, TotalLz2, NbrFluxQuanta2, ReferenceState);              
	}
      else
	{
	  InputSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles2, TotalLz2, NbrFluxQuanta2, ReferenceState);              
	}
    }
  ParticleOnSphere* OutputSpace = 0;
  if (Statistics2 == true)
    {
      if (LzSymmmetryFlag == true)
	{
	  OutputSpace = new FermionOnSphere(NbrParticles1, TotalLzOutputState, NbrFluxQuantumOutputState);
	}
      else
	{
	  OutputSpace = new FermionOnSphere(NbrParticles1, TotalLzOutputState, NbrFluxQuantumOutputState);
	}
    }
  else
    {
      if (LzSymmmetryFlag == true)
	{
	  OutputSpace = new BosonOnSphereShort(NbrParticles1, TotalLzOutputState, NbrFluxQuantumOutputState);
	}
      else
	{
	  OutputSpace = new BosonOnSphereShort(NbrParticles1, TotalLzOutputState, NbrFluxQuantumOutputState);
	}
    }

  char* GeometryName = new char[128];
  if (Manager.GetBoolean("normalize") == true)
    {
      sprintf (GeometryName, "sphere");
    }
  else
    {
      sprintf (GeometryName, "unnormalized");
    }

  char* OutputName = new char [512  + strlen(DiscreteSymmetryName)+ strlen(Manager.GetString("interaction-name")) + strlen(GeometryName)];
  if (Statistics2 == true)
    {
      sprintf (OutputName, "fermions_%s%s_%s_n_%d_2s_%d_lz_%d.%ld.vec", GeometryName, DiscreteSymmetryName, Manager.GetString("interaction-name"), 
	       NbrParticles1, NbrFluxQuantumOutputState, TotalLzOutputState, Manager.GetInteger("outputvector-index"));
    }
  else
    {
      sprintf (OutputName, "bosons_%s%s_%s_n_%d_2s_%d_lz_%d.%ld.vec", GeometryName, DiscreteSymmetryName, Manager.GetString("interaction-name"), 
	       NbrParticles1, NbrFluxQuantumOutputState, TotalLzOutputState, Manager.GetInteger("outputvector-index"));
    }

  if (Architecture.GetArchitecture()->CanWriteOnDisk())
    {
      cout << "generating state " << OutputName << endl;
    }


  char* PartialProductPrefix = 0;
  if (Manager.GetBoolean("save-partial") == true)
    {
      PartialProductPrefix = ReplaceString(Manager.GetString("state-2"), ".vec", "_partialllltimeslll");
    }

  RealVector OutputVector;
  if (Statistics2 == true)
    {
      FQHESphereBosonicStateTimesFermionicStateOperation Operation (BosonicInputVector, InputVector, 
								    BosonicInputSpace, (FermionOnSphere*) InputSpace, (FermionOnSphere*) OutputSpace,
								    !(Manager.GetBoolean("normalize")), PartialProductPrefix);
      Operation.ApplyOperation(Architecture.GetArchitecture());
      OutputVector = Operation.GetState();
    }
  else
    {
      FQHESphereBosonicStateTimesFermionicStateOperation Operation (BosonicInputVector, InputVector, 
								    BosonicInputSpace, (BosonOnSphereShort*) InputSpace, (BosonOnSphereShort*) OutputSpace,
								    !(Manager.GetBoolean("normalize")), PartialProductPrefix);
      Operation.ApplyOperation(Architecture.GetArchitecture());
      OutputVector = Operation.GetState();
    }
  delete InputSpace;
  delete BosonicInputSpace;

  if ((Architecture.GetArchitecture()->CanWriteOnDisk()) && (Manager.GetBoolean( "save-partial") == false))
    {
      if (Manager.GetBoolean("normalize") == true)
	OutputVector.Normalize();
      if (OutputVector.WriteVector(OutputName) == false)
	{
	  cout << "can't write " << OutputName << endl;
	}
    }

  delete[] OutputName;
  
  delete OutputSpace;

  return 0;
}

