#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"

#include "Hamiltonian/ParticleOnSphereGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereTwoLandauLevelDeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnSphereMainTask.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <cstring>
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
  OptionManager Manager ("FQHESphereBosonsTwoLandauLevelTwoBodyGeneric" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  
  ArchitectureManager Architecture;  
  LanczosManager Lanczos(false);
  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 10);
  (*SystemGroup) += new SingleIntegerOption  ('z', "initial-lz", "initial lz value", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "potential-lll", "one particle potential lll", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "potential-sll", "one particle potential sll", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file-lll", "file describing the 2-body interaction in terms of the pseudo-potential on LLL");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file-sll", "file describing the 2-body interaction in terms of the pseudo-potential on SLL");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
#endif

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");  

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsTwoLandauLevelTwoBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
   {
      Manager.DisplayHelp (cout);
      return 0;
   }
  
  // Get options for the options manager.  
  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");
  int InitialLz = Manager.GetInteger("initial-lz");
  long Memory = Manager.GetInteger("memory") << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  double cyclotron_energies[2];
  cyclotron_energies[0] = (double)Manager.GetInteger("potential-lll");
  cyclotron_energies[1] = (double)Manager.GetInteger("potential-sll");
  
  bool FirstRun = true;
  
  // Read in Pseudopotentials from the interaction file. 
  double** PseudoPotentials ;
  PseudoPotentials = new double*[4];
  for ( int i = 0 ; i < 4 ;  i++ ) {
      PseudoPotentials[i] = new double[NbrFluxQuanta+3]; //make enough space for each level
      for ( int j = 0 ; j < NbrFluxQuanta+3 ; j++ ) 
        {
	  PseudoPotentials[i][j] = 0.0;
        }
  }
  
  if ( ( Manager.GetString("interaction-file-lll") == 0 ) && ( Manager.GetString("interaction-file-sll") == 0 ) ) 
    {
      cout << "No interaction file supplied. Will use delta interaction." << endl;
    }
  else
    {
      if ( Manager.GetString("interaction-file-lll") != 0 )
        {
	  ConfigurationParser InteractionDefinition;
	  if (InteractionDefinition.Parse(Manager.GetString("interaction-file-lll")) == false)
	    {
	      InteractionDefinition.DumpErrors(cout) << endl;
	      return -1;
	    }
	  int TmpNbrPseudoPotentials;
	  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials[0], TmpNbrPseudoPotentials) == false)
	    {
	      cout << "Weights is not defined or as a wrong value in " << Manager.GetString("interaction-file-lll") << endl;
	      return -1;
	    }
	  if (TmpNbrPseudoPotentials != (NbrFluxQuanta +1))
	    {
	      cout << "Invalid number of pseudo-potentials: " << TmpNbrPseudoPotentials <<  endl;
	      return -1;	  
	    }
	}
	
      if ( Manager.GetString("interaction-file-sll") != 0 )
        {
	  ConfigurationParser InteractionDefinition;
	  if (InteractionDefinition.Parse(Manager.GetString("interaction-file-sll")) == false)
	    {
	      InteractionDefinition.DumpErrors(cout) << endl;
	      return -1;
	    }
	  int TmpNbrPseudoPotentials;
	  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials[1], TmpNbrPseudoPotentials) == false)
	    {
	      cout << "Weights is not defined or as a wrong value in " << Manager.GetString("interaction-file-sll") << endl;
	      return -1;
	    }
	  if (TmpNbrPseudoPotentials != (NbrFluxQuanta +3))
	    {
	      cout << "Invalid number of pseudo-potentials: " << TmpNbrPseudoPotentials <<  endl;
	      return -1;	  
	    }
	}
	
    }

  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  sprintf (OutputNameLz, "bosons_%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), NbrBosons, NbrFluxQuanta);
  int Max = ((NbrFluxQuanta+2) * NbrBosons); //the maximum possible TotalLz value is when all the bosons are placed in the SLL on the far right.
  int  L = InitialLz;

  
  if ((abs(Max) & 1) != (InitialLz & 1))
    L += 1;
  
  if ((abs(Max) & 1) != (InitialLz & 1)) //since total Lz goes in steps of 2 make sure we start at the right parity.
    L += 1;
  
  if (GroundFlag == true)
    Max = L;

  
  for (; L <= Max; L += 2)
    {
      // create space.
      ParticleOnSphereWithSpin* Space = new BosonOnSphereTwoLandauLevels(NbrBosons, L, NbrFluxQuanta+2, NbrFluxQuanta);    
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	  Memory = Architecture.GetArchitecture()->GetLocalMemory();
      
      // Create data structure for the Hamiltonian.
      AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
      
      if ( ( Manager.GetString("interaction-file-lll") == 0 ) && ( Manager.GetString("interaction-file-sll") == 0 ) ) 
	{
	  Hamiltonian = new ParticleOnSphereTwoLandauLevelDeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta+2, NULL, cyclotron_energies,  
									   Architecture.GetArchitecture(), 
									   Memory, DiskCacheFlag,
									   LoadPrecalculationFileName);
	}
      else 
	{
	  Hamiltonian = new ParticleOnSphereTwoLandauLevelDeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta+2, PseudoPotentials, cyclotron_energies,  
									   Architecture.GetArchitecture(), 
									   Memory, DiskCacheFlag,
									   LoadPrecalculationFileName);
	}

      
      double Shift = 0;
      int LzMax = NbrFluxQuanta + 2;
      
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "bosons_2ll_%s_n_%d_2s_%d_lz_%d", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrBosons, NbrFluxQuanta, L);
	}
      
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
      MainTaskOperation TaskOperation (&Task); TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      
      if (FirstRun == true)
	FirstRun = false;
    }

  return 0;
}
