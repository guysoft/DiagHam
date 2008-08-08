#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "Hamiltonian/ParticleOnSphereGenericHamiltonian.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

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

  // some running options and help
  OptionManager Manager ("FQHEFermionsTwoBodyGeneric" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 8);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "get-lvalue", "compute mean l value from <L^2> for each eigenvalue");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

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
  (*LanczosGroup) += new  BooleanOption ('\n', "resume-fastdisk", "resume the fast-disk mode Lanczos algorithm from a stopped one (for example due to computer crash)");
  (*LanczosGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEFermionsTwoBodyGeneric -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  bool HaldaneBasisFlag = ((BooleanOption*) Manager["haldane"])->GetBoolean();
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();  
  bool DiskCacheFlag = ((BooleanOption*) Manager["disk-cache"])->GetBoolean();
  bool FirstRun = true;
  double* PseudoPotentials = 0;
  double* OneBodyPotentials = 0;
  if (((SingleStringOption*) Manager["interaction-file"])->GetString() == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(((SingleStringOption*) Manager["interaction-file"])->GetString()) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      int TmpNbrPseudoPotentials;
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials, TmpNbrPseudoPotentials) == false)
	{
	  cout << "Weights is not defined or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	  return -1;
	}
      cout << "LzMax= " << TmpNbrPseudoPotentials << " " << LzMax << endl;
      if (TmpNbrPseudoPotentials != (LzMax +1))
	{
	  cout << "Invalid number of pseudo-potentials" << endl;
	  return -1;	  
	}
      if (InteractionDefinition.GetAsDoubleArray("Onebodypotentials", ' ', OneBodyPotentials, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials" << endl;
	      return -1;	  
	    }
	}
    }

  char* OutputNameLz = new char [256 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
  sprintf (OutputNameLz, "fermions_%s_n_%d_2s_%d_lz.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax);

  int Max = ((LzMax - NbrParticles + 1) * NbrParticles);

  int  L = InitialLz;
  if (L < -Max)
    L = -Max;
  else
    if (L > Max)
      L = Max;
  if ((abs(Max) & 1) != (abs(InitialLz) & 1))
    L += 1;
  if (GroundFlag == true)
      Max = L;
  else
    {
      if (NbrLz > 0)
	{
	  Max = L + (2 * (NbrLz - 1));
	}
    }
  for (; L <= Max; L += 2)
    {
      ParticleOnSphere* Space = 0;
      if (HaldaneBasisFlag == false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	  if (LzMax <= 30)
#endif
	    if ((SymmetrizedBasis == false) || (L != 0))
	      Space = new FermionOnSphere(NbrParticles, L, LzMax, MemorySpace);
	    else
	      {
		if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		  Space = new FermionOnSphereSymmetricBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		else
		  Space = new FermionOnSphereSymmetricBasis(NbrParticles, LzMax, MemorySpace);
		if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
		  {
		    ((FermionOnSphereSymmetricBasis*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
		    return 0;
		  }
	      }
	  else
#ifdef __128_BIT_LONGLONG__
	    if (LzMax <= 126)
#else
	      if (LzMax <= 62)
#endif
		{
		  if ((SymmetrizedBasis == false) || (L != 0))
		    Space = new FermionOnSphereLong(NbrParticles, L, LzMax, MemorySpace);
		  else
		    {
		      if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
			Space = new FermionOnSphereSymmetricBasisLong(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		      else
			Space = new FermionOnSphereSymmetricBasisLong(NbrParticles, LzMax, MemorySpace);
		      if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
			{
			  ((FermionOnSphereSymmetricBasisLong*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
			  return 0;
			}
		    }
		}
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, L, LzMax, MemorySpace);
	}
      else
	{
	  int* ReferenceState = 0;
	  if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
	    {
	      ReferenceState = new int[LzMax + 1];
	      for (int i = 0; i <= LzMax; ++i)
		ReferenceState[i] = 0;
	      if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "laughlin") == 0)
		for (int i = 0; i <= LzMax; i += 3)
		  ReferenceState[i] = 1;
	      else
		if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "pfaffian") == 0)
		  for (int i = 0; i <= LzMax; i += 4)
		    {
		      ReferenceState[i] = 1;
		      ReferenceState[i + 1] = 1;
		    }
		else
		  if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "readrezayi3") == 0)
		    for (int i = 0; i <= LzMax; i += 5)
		      {
			ReferenceState[i] = 1;
			ReferenceState[i + 1] = 1;
			ReferenceState[i + 2] = 1;
		      }
		  else
		    {
		      cout << "unknown reference state " << ((SingleStringOption*) Manager["reference-state"])->GetString() << endl;
		      return -1;
		    }
	    }
	  else
	    {
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
		{
		  ReferenceStateDefinition.DumpErrors(cout) << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
		{
		  cout << "NbrParticles is not defined or as a wrong value" << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
		{
		  cout << "LzMax is not defined or as a wrong value" << endl;
		  return -1;
		}
	      int MaxNbrLz;
	      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		{
		  cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
		  return -1;     
		}
	      if (MaxNbrLz != (LzMax + 1))
		{
		  cout << "wrong LzMax value in ReferenceState" << endl;
		  return -1;     
		}
	    }
	  if (SymmetrizedBasis == false)
	     {
#ifdef __64_BITS__
	       if (LzMax <= 62)
#else
		 if (LzMax <= 30)
#endif
		   {
		     if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		       Space = new FermionOnSphereHaldaneBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		     else
		       Space = new FermionOnSphereHaldaneBasis(NbrParticles, L, LzMax, ReferenceState, MemorySpace);
		     if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
		       {
			 ((FermionOnSphereHaldaneBasis*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
			 return 0;
		       }
		   }
	       else
#ifdef __128_BIT_LONGLONG__
		 if (LzMax <= 126)
#else
		   if (LzMax <= 62)
#endif
		     {
		       if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
			 Space = new FermionOnSphereHaldaneBasisLong(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		       else
			 Space = new FermionOnSphereHaldaneBasisLong(NbrParticles, L, LzMax, ReferenceState, MemorySpace);
		       if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
			 {
			   ((FermionOnSphereHaldaneBasisLong*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
			   return 0;
			 }
		     }	       
	     }
	  else
	    {
#ifdef __64_BITS__
	       if (LzMax <= 62)
#else
		 if (LzMax <= 30)
#endif
		   {
		     if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		       Space = new FermionOnSphereHaldaneSymmetricBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		     else
		       Space = new FermionOnSphereHaldaneSymmetricBasis(NbrParticles, LzMax, ReferenceState, MemorySpace);
		     if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
		       {
			 ((FermionOnSphereHaldaneSymmetricBasis*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
			 return 0;
		       }
		   }
		 else
#ifdef __128_BIT_LONGLONG__
		   if (LzMax <= 126)
#else
		     if (LzMax <= 62)
#endif
		       {
			 if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
			   Space = new FermionOnSphereHaldaneSymmetricBasisLong(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
			 else
			   Space = new FermionOnSphereHaldaneSymmetricBasisLong(NbrParticles, LzMax, ReferenceState, MemorySpace);
			 if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
			   {
			     ((FermionOnSphereHaldaneSymmetricBasisLong*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
			     return 0;
			   }
		       }
	    }
	}
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      if (OneBodyPotentials == 0)
	Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials,
							     ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(),
							     Architecture.GetArchitecture(), 
							     Memory, DiskCacheFlag,
							     LoadPrecalculationFileName);
      else
	Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials, OneBodyPotentials,
							     ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(),
							     Architecture.GetArchitecture(), 
							     Memory, DiskCacheFlag,
							     LoadPrecalculationFileName);
      if (Manager.GetString("energy-expectation") != 0 )
	{
	  char* StateFileName = Manager.GetString("energy-expectation");
	  if (IsFile(StateFileName) == false)
	    {
	      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
	      return -1;           
	    }
	  RealVector State;
	  if (State.ReadVector(StateFileName) == false)
	    {
	      cout << "error while reading " << StateFileName << endl;
	      return -1;
	    }
	  if (State.GetVectorDimension()!=Space->GetHilbertSpaceDimension()) break;
	  RealVector TmpState(Space->GetHilbertSpaceDimension());
	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  double EnergyValue = TmpState * State;
	  cout << "< Energy > = "<<EnergyValue<<endl;
	}
      else // normal procedure: diagonalize Hamiltonian
	{
	  double Shift = - 0.5 * ((double) (NbrParticles * NbrParticles)) / (0.5 * ((double) LzMax));
	  Hamiltonian->ShiftHamiltonian(Shift);
	  char* EigenvectorName = 0;
	  if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	    {
	      EigenvectorName = new char [64];
	      sprintf (EigenvectorName, "fermions_%s_n_%d_2s_%d_lz_%d", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax, L);
	    }
	  
	  QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  if (EigenvectorName != 0)
	    {
	      delete[] EigenvectorName;
	    }
	}
      delete Hamiltonian;
      if (FirstRun == true)
	FirstRun = false;
    }

  return 0;
}
