#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereMonomialsTimesSlaterProjectionOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereMultipleMonomialsTimesSlaterProjectionOperation.h"


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
  OptionManager Manager ("FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  ArchitectureManager Architecture;
  Architecture.AddOptionGroup(&Manager);
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new BooleanOption ('\n',"reverse-flux","the fluxes bind to each particle are in the opposite direction than the magnetic field");

  (*OutputGroup) += new BooleanOption ('\n', "normalize", "normalize the projected state assuming the sphere geometry");
  (*OutputGroup) += new BooleanOption  ('\n', "rational" , "use rational numbers instead of double precision floating point numbers");
  (*OutputGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*OutputGroup) += new SingleIntegerOption ('\n', "outputvector-index", "set the index of the output vector (i.e. the integer in the extention *.xxx.vec)", 0);
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL -h" << endl;
      return -1;
    }	
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrParticles  = Manager.GetInteger("nbr-particles");
  int TotalSz = 0;
  int NbrParticleUp = (NbrParticles + TotalSz) / 2;
  int NbrParticleDown = (NbrParticles - TotalSz) / 2;  
  int NbrLandauLevel = 1;
  int NbrFluxQuantumLambdaLevels = NbrParticleUp - (NbrLandauLevel * (NbrLandauLevel - 1));
  if ((NbrFluxQuantumLambdaLevels % NbrLandauLevel) != 0)
    {
      cout << "error, the number of particles is not compatible with " <<  NbrLandauLevel << " filled Landau levels" << endl;
      return -1;
    }
  NbrFluxQuantumLambdaLevels /= NbrLandauLevel;
  --NbrFluxQuantumLambdaLevels;
  int NbrFluxQuanta = (NbrParticles - 1);
  if (Manager.GetBoolean("reverse-flux") == false)
    NbrFluxQuanta += NbrFluxQuantumLambdaLevels;
  else
    NbrFluxQuanta -= NbrFluxQuantumLambdaLevels;
  int TotalLz = 0;
    
  FermionOnSphereWithSpin* InputSpace = new FermionOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuantumLambdaLevels, TotalSz);

  BosonOnSphereWithSU2Spin* OutputSpace = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);

  char* OutputName = new char [512 + strlen(Manager.GetString("interaction-name"))];
  sprintf (OutputName, "bosons_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d.%ld.vec", Manager.GetString("interaction-name"), NbrParticles, NbrFluxQuanta, TotalSz, TotalLz, Manager.GetInteger("outputvector-index"));

  cout << "generating state " << OutputName << endl;
  RealVector InputVector (InputSpace->GetHilbertSpaceDimension(), true);
  RealVector OutputVector (OutputSpace->GetHilbertSpaceDimension(), true);
  
  InputVector[0] = 1.0;
  OutputSpace->SlaterTimeSpinfulFermionicState(InputVector, OutputVector, InputSpace, 0, InputSpace->GetHilbertSpaceDimension(),
					       !(Manager.GetBoolean("normalize")));

  if (Manager.GetBoolean("normalize") == true)
    OutputVector.Normalize();
  if (OutputVector.WriteVector(OutputName) == false)
    {
      cout << "can't write " << OutputName << endl;
    }
  delete[] OutputName;
  
  delete OutputSpace;
  delete InputSpace;

  return 0;
}

