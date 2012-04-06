#include "Vector/RealVector.h"

#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/GroupedPermutations.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"


#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereBosonsWithSpinLandauLevelLiftOperation.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  
  OptionManager Manager ("FQHESphereBosonicStateLandauLevelLift" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  ArchitectureManager Architecture;
	
  Manager += SystemGroup;
  Manager += MiscGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += OutputGroup;
	
  (*SystemGroup) += new SingleStringOption ('\0', "state", "name of a spinful bosonic state");
  (*SystemGroup) += new SingleStringOption ('p', "polarized-state", "name of a polarized bosonic state");

  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files, default is nbody)");
	// (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  // (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  // (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "resume-file", "use this file as the partial vector to resume from");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "resume-idx", "use this file as the partial vector to resume from", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "mpi-stages", "the number of stages divide into when using MPI  (default is 20)", 20);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "smp-stages", "the number of stages divide into when using SMP  (default is 20)", 20);
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "the output vector will be written in the normalized basis",false);
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the Jack polynomial decomposition into a binary file");
  //(*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the Jack polynomial decomposition into a text file");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);


  // testing:
  GroupedPermutations TestPer(2,2,false);
  SmallIntegerArray* Perms = TestPer.GetPermutations();
  cout << "total number of grouped permutations (of 4+4 particles) : "<<TestPer.GetNbrPermutations()<<endl;
  
  for (int i=0; i<TestPer.GetNbrPermutations(); ++i)
    cout << i<<"\t"<<Perms[i]<<endl;
  exit(1);
 	
  RealVector InitialState;
  if (InitialState.ReadVector(Manager.GetString("state")) == false)
    {
      cout << "error while reading " << Manager.GetString("state") << endl;
      return -1;
    }

  RealVector PolarizedState;
  if (PolarizedState.ReadVector(Manager.GetString("polarized-state")) == false)
    {
      cout << "error while reading " << Manager.GetString("polarized-state") << endl;
      return -1;
    }

 
  int NbrParticles = 0, NbrParticles2=0;
  int LzMaxIn = 0;
  int TotalLzIn = 0;
  int TotalSzIn = 0;
  bool FermionFlag = false;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMaxIn, TotalLzIn, TotalSzIn, FermionFlag) == false)
   {
     cout << "Error: Could not read system information from input file name"<<endl;
     return -1;
   }
  if (FermionFlag==true)
    {
      cout << "Error: bosonic input file required!"<<endl;
      return -1;
    }
  
  int LzMaxPo = 0;
  int TotalLzPo = 0;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("polarized-state"), NbrParticles2, LzMaxPo, TotalLzPo, FermionFlag) == false)
    {
      cout << "Error: Could not read system information from polarized file name"<<endl;
      return -1;
    }
  if (FermionFlag==true)
    {
      cout << "Error: bosonic input file required!"<<endl;
      return -1;
    }
  if (NbrParticles2!=NbrParticles)
    {
      cout << "Error: Initial state and polarized state need to have the same number of particles"<<endl;
      return -1;
    }
  
  BosonOnSphereWithSpin * InitialSpace;
  InitialSpace = new BosonOnSphereWithSpin (NbrParticles,TotalLzIn,LzMaxIn,TotalSzIn);
  
  if (InitialSpace->GetHilbertSpaceDimension() != InitialState.GetVectorDimension())
    {
      cout << "dimension mismatch between the initial state (" << InitialState.GetVectorDimension() << ") and the Hilbert space (" << InitialSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  
  BosonOnSphereShort * PolarizedSpace;
  PolarizedSpace = new BosonOnSphereShort (NbrParticles,TotalLzPo,LzMaxPo);
  
  if (PolarizedSpace->GetHilbertSpaceDimension() != PolarizedState.GetVectorDimension())
    {
      cout << "dimension mismatch between the polarized state (" << PolarizedState.GetVectorDimension() << ") and the Hilbert space (" << PolarizedSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }


  int TotalLzFi = TotalLzIn + TotalLzPo;
  int LzMaxFi = LzMaxIn + LzMaxPo;
  BosonOnSphereWithSpin * FinalSpace;

  FinalSpace = new BosonOnSphereWithSpin (NbrParticles,TotalLzFi,LzMaxFi,TotalSzIn);

    
  RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
  int ResumeIdx = Manager.GetInteger("resume-idx");
		
  if (Manager.GetString("resume-file") != 0 )
    {
      if ( Architecture.GetArchitecture()->ReadVector(OutputVector, Manager.GetString("resume-file")) == false )
	{
	  cout << "error while reading " << Manager.GetString("resume-file") << endl;
	  return -1;
	}		
    }
  

  // process operation tasks:
  
  FQHESphereBosonsWithSpinLandauLevelLiftOperation MainOperation(InitialSpace, PolarizedSpace, FinalSpace, &InitialState, &PolarizedState, &OutputVector, Manager.GetInteger("mpi-stages"), Manager.GetInteger("smp-stages"), ResumeIdx);
  MainOperation.ApplyOperation(Architecture.GetArchitecture());    
  
      
  if(Manager.GetBoolean("normalize"))
    FinalSpace->ConvertFromUnnormalizedMonomial(OutputVector,0l,true);
      
  Architecture.GetArchitecture()->WriteVector(OutputVector, Manager.GetString("bin-output"));

  delete InitialSpace;
  delete PolarizedSpace;
  delete FinalSpace;
}


