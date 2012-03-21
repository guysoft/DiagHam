#include "HilbertSpace/BosonOnSphereWithSU4Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU4SpinAllEntanglement.h"

#include "Hamiltonian/ParticleOnSphereWithSU4SpinS2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSU4SpinL2Hamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"

#include "LanczosAlgorithm/LanczosManager.h"
#include "HilbertSpace/ParticleOnSphereManager.h"

#include <iostream>
#include <stdlib.h>
#include <cstring>
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
  OptionManager Manager ("FQHESphereWithSU4SpinL2Diagonalize" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  ParticleOnSphereManager ParticleManager(true, true, 4);
  ParticleManager.AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");  
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  
  Manager += LanczosGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "l2");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "energy-shift", "if non zero, override energy shift using the indicated value ", -10.0);

  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag",
						"maximum Hilbert space dimension for which full diagonalization is applied",
						500, true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);  
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);
  (*LanczosGroup) += new SingleDoubleOption ('\n', "lanczos-precision", "define Lanczos precision for eigenvalues (0 if automatically defined by the program)", 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "fast-disk", "use disk storage to increase speed of ground state calculation and decrease memory footprint when using Lanczos algorithm");
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the haldane or symmetrized bases)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the haldane or symmetrized bases)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz  = Manager.GetInteger("total-lz");
  int TotalSz  = Manager.GetInteger("total-sz");
  int IsoSzTotal = Manager.GetInteger("total-isosz");
  int TotalEntanglement = Manager.GetInteger("total-entanglement");

//   bool LzSymmetrizedBasis = ((BooleanOption*) Manager["lzsymmetrized-basis"])->GetBoolean();
//   bool TzSymmetrizedBasis = ((BooleanOption*) Manager["tzsymmetrized-basis"])->GetBoolean();
//   bool Z3SymmetrizedBasis = ((BooleanOption*) Manager["z3symmetrized-basis"])->GetBoolean();
//   bool TzMinusParity = ((BooleanOption*) Manager["minus-tzparity"])->GetBoolean();
  unsigned long Memory = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];

  if (Manager.GetBoolean("use-entanglement"))
    sprintf (OutputNameLz, "fermions_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_pz_%d_lz.dat", Manager.GetString("interaction-name"), 
	     NbrParticles, LzMax, TotalSz, IsoSzTotal, TotalEntanglement);
  else
    sprintf (OutputNameLz, "fermions_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_lz.dat", Manager.GetString("interaction-name"), 
	     NbrParticles, LzMax, TotalSz, IsoSzTotal);



  ParticleOnSphereWithSU4Spin* Space=(ParticleOnSphereWithSU4Spin*) ParticleManager.GetHilbertSpace(TotalLz);

  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
  AbstractQHEOnSphereHamiltonian* Hamiltonian = new ParticleOnSphereWithSU4SpinL2Hamiltonian(Space, NbrParticles, LzMax, TotalLz,
											     Architecture.GetArchitecture(), 
											     ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(),
											     Memory, DiskCacheFlag,
											     LoadPrecalculationFileName);

  double Shift = ((SingleDoubleOption*) Manager["energy-shift"])->GetDouble();
  Hamiltonian->ShiftHamiltonian(Shift);
  char* EigenvectorName = 0;
  if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
    {
      EigenvectorName = new char [128 + strlen(Manager.GetString("interaction-name"))];
      if (Manager.GetBoolean("use-entanglement"))
	sprintf (EigenvectorName, "fermions_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_pz_%d_lz", Manager.GetString("interaction-name"), 
		 NbrParticles, LzMax, TotalSz, IsoSzTotal, TotalEntanglement);
      else
	sprintf (EigenvectorName, "fermions_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_lz", Manager.GetString("interaction-name"), 
	     NbrParticles, LzMax, TotalSz, IsoSzTotal);
    }
  QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, TotalLz, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
  MainTaskOperation TaskOperation (&Task);
  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
  if (EigenvectorName != 0)
    delete[] EigenvectorName;
  delete Hamiltonian;
  if (FirstRun == true)
    FirstRun = false;

  return 0;
}
