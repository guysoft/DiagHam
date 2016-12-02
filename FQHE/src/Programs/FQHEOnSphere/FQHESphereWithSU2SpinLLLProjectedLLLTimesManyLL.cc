#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSzSymmetry.h"

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
  (*SystemGroup) += new BooleanOption ('\n',"disable-szsymmetry","disable the Sz<->-Sz symmetry for the Sz=0 sector");
  (*SystemGroup) += new BooleanOption ('\n',"disable-lzsymmetry","disable the Lz<->-Lz symmetry for the Lz=0 sector");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");

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
    
  bool LzSymmmetryFlag = true;
  if ((TotalLz != 0) || (Manager.GetBoolean("disable-lzsymmetry") == true))
    {
      LzSymmmetryFlag = false;
    }
  bool SzSymmmetryFlag = true;
  if ((TotalSz != 0) || (Manager.GetBoolean("disable-szsymmetry") == true))
    {
      SzSymmmetryFlag = false;
    }
  char* DiscreteSymmetryName = new char[128];
  if (SzSymmmetryFlag == false)
    {
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
    }
  else
    {
      if (LzSymmmetryFlag == false)
	{
	  if (Manager.GetBoolean("minus-szparity") == false)
	    {
	      sprintf (DiscreteSymmetryName, "_szsym_1");
	    }
	  else
	    {
	      sprintf (DiscreteSymmetryName, "_szsym_-1");
	    }
	}
      else
	{
	  if (Manager.GetBoolean("minus-szparity") == false)
	    {
	      if (Manager.GetBoolean("minus-lzparity") == false)
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_1_szsym_1");
		}
	      else
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_-1_szsym_1");
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("minus-lzparity") == false)
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_1_szsym_-1");
		}
	      else
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_-1_szsym_-1");
		}
	    }
	}
    }


  FermionOnSphereWithSpin* InputSpace = new FermionOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuantumLambdaLevels, TotalSz);

  BosonOnSphereWithSU2Spin* OutputSpace = 0;
  if (LzSymmmetryFlag == true)
    {
      if (SzSymmmetryFlag == true)
	{
	  OutputSpace = new BosonOnSphereWithSU2SpinLzSzSymmetry(NbrParticles, NbrFluxQuanta, TotalSz, Manager.GetBoolean("minus-szparity"), 
								  Manager.GetBoolean("minus-lzparity"));
	}
      else
	{
	  OutputSpace = new BosonOnSphereWithSU2SpinLzSymmetry(NbrParticles, NbrFluxQuanta, TotalSz, Manager.GetBoolean("minus-lzparity"));
	}
    }
  else
    {
      if (SzSymmmetryFlag == true)
	{
	  OutputSpace = new BosonOnSphereWithSU2SpinSzSymmetry(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, Manager.GetBoolean("minus-szparity"));
	}
      else
	{
	  OutputSpace = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
	}
    }

  char* OutputName = new char [512  + strlen(DiscreteSymmetryName)+ strlen(Manager.GetString("interaction-name"))];
  sprintf (OutputName, "bosons_sphere_su2%s_%s_n_%d_2s_%d_sz_%d_lz_%d.%ld.vec", DiscreteSymmetryName, Manager.GetString("interaction-name"), NbrParticles, NbrFluxQuanta, TotalSz, TotalLz, Manager.GetInteger("outputvector-index"));

  if (Architecture.GetArchitecture()->CanWriteOnDisk())
    {
      cout << "generating state " << OutputName << endl;
    }
  RealVector InputVector (InputSpace->GetHilbertSpaceDimension(), true);
  RealVector OutputVector (OutputSpace->GetHilbertSpaceDimension(), true);
  
  InputVector[0] = 1.0;
  OutputSpace->SlaterTimeSpinfulFermionicState(InputVector, OutputVector, InputSpace, 0, InputSpace->GetHilbertSpaceDimension(),
					       !(Manager.GetBoolean("normalize")), Architecture.GetArchitecture());

  if (Architecture.GetArchitecture()->CanWriteOnDisk())
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
  delete InputSpace;

  return 0;
}

