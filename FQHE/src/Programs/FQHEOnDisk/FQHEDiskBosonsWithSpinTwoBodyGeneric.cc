#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnDiskWithSU2Spin.h"

#include "Hamiltonian/ParticleOnSphereWithSpinGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#ifdef __MPI__
#include <mpi.h>
#endif


using std::ios;
using std::cout;
using std::endl;

using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);


  OptionManager Manager ("FQHEDiskBosonsWithSpinTwoBodyGeneric" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('l', "maximum-momentum", "maximum total angular momentum to study", 10, true, 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "minimum-momentum", "minimum total angular momentum to study", 1, true, 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "force-maxmomentum", "force the maximum single particle momentum to a particular value (negative from the number of particles and the state total angular momentum)", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-exponent", "start the Haldane algorithm from Laughlin state with exponent m)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new BooleanOption  ('\n', "allow-disk-storage", "expand memory for fast multiplication using disk storage",false);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the Hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskBosonsWithSpinTwoBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int MMin = Manager.GetInteger("minimum-momentum");
  int MMax = Manager.GetInteger("maximum-momentum");
  if (MMin < (((NbrParticles - 1) * (NbrParticles)) / 2))
    MMin = (((NbrParticles - 1) * (NbrParticles)) / 2);
  if (MMax < MMin)
    MMax = MMin;
  int ForceMaxMomentum = Manager.GetInteger("force-maxmomentum");
  int SzTotal = Manager.GetInteger("total-sz");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;  
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  bool onDiskCacheFlag = Manager.GetBoolean("allow-disk-storage");
  bool FirstRun = true;
  double** PseudoPotentials  = new double*[10];
  for (int i = 0; i < 3; ++i)
    {
      PseudoPotentials[i] = new double[LzMax + 1];
      for (int j = 0; j <= LzMax; ++j)
	PseudoPotentials[i][j] = 0.0;
    };
  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;

  int NbrUp = (NbrBosons + SzTotal) >> 1;
  int NbrDown = (NbrBosons - SzTotal) >> 1;
  if ((NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }

  if (Manager.GetString("interaction-file") == 0)
    {
      if (!Manager.GetBoolean("l2-s2-only"))
	{
	  cout << "an interaction file has to be provided" << endl;
	  return -1;
	}
    }
  else
    {
      if (FQHEDiskSU2GetPseudopotentials(Manager.GetString("interaction-file"), LzMax, PseudoPotentials,
					 OneBodyPotentialUpUp, OneBodyPotentialDownDown) == false)
	return -1;
    }

  char* OutputNameLz = new char [512 + strlen(Manager.GetString("interaction-name"))];
  if (ForceMaxMomentum >= 0)
    sprintf (OutputNameLz, "bosons_disk_su2_%s_n_%d_lzmax_%d_lz_%d_sz_%d.dat", Manager.GetString("interaction-name"),
	     NbrBosons, ForceMaxMomentum, MMax, SzTotal);
  else
    sprintf (OutputNameLz, "bosons_disk_su2_%s_n_%d_lzmax_%d_lz_%d_sz_%d.dat", Manager.GetString("interaction-name"),
	     NbrBosons, MMax, SzTotal);

  for (int  L = MMin; L <= MMax; ++L)
    {
      int TmpMaxMomentum = L;
      if ((ForceMaxMomentum >= 0) && (ForceMaxMomentum < TmpMaxMomentum))
	TmpMaxMomentum = ForceMaxMomentum;
      ParticleOnSphereWithSpin* Space = 0;
      Space = new FermionOnDisk(NbrParticles, L, TmpMaxMomentum, MemorySpace);

      double Shift = 0.0;
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
        Memory = Architecture.GetArchitecture()->GetLocalMemory();

      AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian;

      Hamiltonian = new ParticleOnDiskWithSpinGenericHamiltonian(Space, NbrBosons, TmpMaxMomentum, PseudoPotentials, OneBodyPotentialUpUp, OneBodyPotentialDownDown, NULL, 
								 Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);
       
      Hamiltonian->ShiftHamiltonian(Shift);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [120];
	  sprintf (EigenvectorName, "bosons_disk_su2_%s_n_%d_2s_%d_sz_%d_lz_%d",
		   Manager.GetString("interaction-name"),
		   NbrBosons, LzMax, SzTotal, L*LSign);
	}
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      delete Space;      
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	  EigenvectorName = 0;
	}
      if (FirstRun == true)
	FirstRun = false;
      if (HaldaneBasisFlag) return 0; // only one subspace defined...
    }
  delete[] OutputNameLz;
  delete[] ExtraTerms;
  return 0;
}


