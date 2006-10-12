#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"

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

#include <iostream>
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

  OptionManager Manager ("FQHESphereFermionsGraphene" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += ToolsGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", 
						"maximum Hilbert space dimension for which full diagonalization is applied", 
						500, true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "limit-time", "use limit in time instead of a number of lanczos iteration (0 if none, time in seconds)", 0);
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 15);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin (i.e valley SU(2) degeneracy) of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new BooleanOption  ('\n', "allow-disk-storage", "expand memory for fast multiplication using disk storage",false);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsSphereWithSpin -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int SzTotal = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  int IsoSzTotal = ((SingleIntegerOption*) Manager["total-isosz"])->GetInteger();

  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  char* SavePrecalculationFileName = ((SingleStringOption*) Manager["save-precalculation"])->GetString();
  bool onDiskCacheFlag = ((BooleanOption*) Manager["allow-disk-storage"])->GetBoolean();
  bool FirstRun = true;
  

  int NbrUp = (NbrFermions + SzTotal) >> 1;
  int NbrDown = (NbrFermions - SzTotal) >> 1;
  if ((NbrUp+NbrDown != NbrFermions) || ( NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }
  int NbrPlus = (NbrFermions + IsoSzTotal) >> 1;
  int NbrMinus = (NbrFermions - IsoSzTotal) >> 1;
  if ((NbrPlus + NbrMinus != NbrFermions) || (NbrPlus < 0) || (NbrMinus < 0))
    {
      cout << "This value of the isospin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }
  int Max = ((LzMax - NbrUp + 1) * NbrUp) + ((LzMax - NbrDown + 1) * NbrDown);

  int  L = 0;
  if ((abs(Max) & 1) != 0)
     L = 1;
  if (InitialLz >= 0)
    {
      L = InitialLz;
      if ((abs(Max) & 1) != 0)
	L |= 1;
      else
	L &= ~0x1;
    }
  if (NbrLz > 0)
    {
      if (L + (2 * (NbrLz - 1)) < Max)
	Max = L + (2 * (NbrLz - 1));
    }
  for (; L <= Max; L += 2)
    {
      double Shift = -10.0;
      ParticleOnSphereWithSpin* Space;
#ifdef __64_BITS__
      if (LzMax <= 15)
#else
      if (LzMax <= 7)
#endif
        {
          Space = new FermionOnSphereWithSU4Spin(NbrFermions, L, LzMax, SzTotal, IsoSzTotal, MemorySpace);
        }
      else
	{
	  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	  return -1;
	}	
    }
      
//       Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
//       AbstractQHEHamiltonian* Hamiltonian;

//       Hamiltonian = new ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian(Space, NbrFermions, LzMax, V0, V1, 
// 									       Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);
//       Hamiltonian->ShiftHamiltonian(Shift);
//       if (SavePrecalculationFileName != 0)
// 	{
// 	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
// 	}
//       char* EigenvectorName = 0;
//       if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
// 	{
// 	  EigenvectorName = new char [120];
// 	  sprintf (EigenvectorName, "fermions_sphere_spin_n_%d_2S_%d_Sz_%d_lz_%d_V_%g_W_%g.ev",
// 		   NbrFermions, LzMax, SzTotal, L, V0, V1);
// 	}
//       QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName);
//       MainTaskOperation TaskOperation (&Task);
//       TaskOperation.ApplyOperation(Architecture.GetArchitecture());
//       delete Hamiltonian;
//       delete Space;
//       if (EigenvectorName != 0)
// 	{
// 	  delete[] EigenvectorName;
// 	  EigenvectorName = 0;
// 	}
//       if (FirstRun == true)
// 	FirstRun = false; 
//     }
//   delete [] OutputNameLz;
  return 0;
}
