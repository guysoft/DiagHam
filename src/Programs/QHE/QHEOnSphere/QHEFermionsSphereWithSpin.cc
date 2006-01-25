#include "HilbertSpace/QHEHilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/QHEHilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnSphere.h"
#include "Hamiltonian/QHEHamiltonian/AbstractQHEHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEMainTask/QHEOnSphereMainTask.h"

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


using std::ios;
using std::cout;
using std::endl;

using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("QHEFermionsSphereWithSpin" , "0.01");
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
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);  

  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 15);
  (*SystemGroup) += new SingleIntegerOption  ('s', "SzTotal", "twice the z component of the total spin of the system", 0);
  (*SystemGroup) += new SingleDoubleOption  ('v', "V0-Interaction", "Interaction in s-wave channel", 0);
  (*SystemGroup) += new SingleDoubleOption  ('w', "V1-Interaction", "Interaction in p-wave channel", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
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
  
  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int SzTotal = ((SingleIntegerOption*) Manager["SzTotal"])->GetInteger();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  double V0 = ((SingleDoubleOption*) Manager["V0-Interaction"])->GetDouble();
  double V1 = ((SingleDoubleOption*) Manager["V1-Interaction"])->GetDouble();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  char* SavePrecalculationFileName = ((SingleStringOption*) Manager["save-precalculation"])->GetString();
  bool onDiskCacheFlag = ((BooleanOption*) Manager["allow-disk-storage"])->GetBoolean();
  bool FirstRun = true;

  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "fermions_sphere_spin_n_%d_2S_%d_Sz_%d_lz_V_%g_W_%g.dat",
	   NbrFermions, LzMax,SzTotal,V0,V1);

  int NbrUp = (NbrFermions + SzTotal)/2;
  int NbrDown = (NbrFermions - SzTotal)/2;
  if ((NbrUp+NbrDown != NbrFermions) || ( NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of Sz cannot be achieved with this particle number!" << endl;
      exit(5);
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
  if (GroundFlag == true)
      Max = L;
  else
    {
      if (NbrLz > 0)
	{
	  if (L + (2 * (NbrLz - 1)) < Max)
	    Max = L + (2 * (NbrLz - 1));
	}
    }
  for (; L <= Max; L += 2)
    {
      double Shift = -10.0;
      ParticleOnSphereWithSpin* Space;
#ifdef __64_BITS__
      if (LzMax <= 31)
        {
          Space = new FermionOnSphereWithSpin(NbrFermions, L, LzMax, SzTotal, MemorySpace);
        }
      else
	{
	  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	  return -1;
	}	
#else
      if (LzMax <= 15)
        {
          Space = new FermionOnSphereWithSpin(NbrFermions, L, LzMax, SzTotal, MemorySpace);
	}
      else
	{
	  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	  return -1;
	}	
#endif
      
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      /*
      // Some testing of the Hilbert Space...
      if (Space->GetHilbertSpaceDimension() < 100)
	for (int i=0; i<Space->GetHilbertSpaceDimension(); i++)
	  { Space->PrintState(cout, i); cout << endl; }
      else
	for (int i=0; i<100; i++)
	  { Space->PrintState(cout, i); cout << endl; }

      // test Four-Point Operators:
      ParticleOnSphere* Space2=0;
      if (NbrFermions==abs(SzTotal))
	{
	  Space2 = new FermionOnSphere(NbrFermions, L, LzMax, MemorySpace);
	  cout << "For comparison: space without spin:" <<endl;
	        if (Space->GetHilbertSpaceDimension() < 100)
		  {
		    if (Space->GetHilbertSpaceDimension() < 100)
		      for (int i=0; i<Space->GetHilbertSpaceDimension(); i++)
			{ Space2->PrintState(cout, i); cout << endl; }
		    else
		      for (int i=0; i<100; i++)
			{ Space2->PrintState(cout, i); cout << endl; }
		  }
	}
      for (int i=0; i<Space->GetHilbertSpaceDimension(); i++)
	{
	  cout << "Applying operators to: "; Space->PrintState(cout, i); cout << endl;
	  for (int m1=0; m1 <= LzMax; ++m1)
	    for (int m2=0; m2 <= m1; ++m2)
	      for (int m3=0; m3 <= LzMax; ++m3)
		{
		  int m4= m1+m2-m3, ret;
		  double coeff=0.0;
		  ret=Space->AddAddAdAd(i,m1,m2,m3,m4,coeff);
		  printf("Pair (m1=%d,d; m2=%d,d; m3=%d,d; m4=%d,d) : ",m1,m2,m3,m4);
		  if (ret < Space->GetHilbertSpaceDimension())
		    {
		      cout << coeff <<"*";
		      Space->PrintState(cout,ret);
		      cout << endl;
		    }
		  else cout << "void" << endl;
		}
	}
      */
      // introduce interaction and test Hamiltonian!
      AbstractQHEHamiltonian* Hamiltonian;

      Hamiltonian = new ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian(Space, NbrFermions, LzMax, V0, V1, Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);
      Hamiltonian->ShiftHamiltonian(Shift);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "fermions_sphere_spin_n_%d_2S_%d_Sz_%d_lz_V_%g_W_%g.ev",
		   NbrFermions, LzMax,SzTotal,V0,V1);
	}
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      delete Space;
      if (FirstRun == true)
	FirstRun = false; 
    }
  return 0;
}
