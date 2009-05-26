#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"

#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereGetAngularMomentumDecomposition" , "0.01");

  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
 
  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the file that contains the state whose average L value has to be evaluated");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleStringOption  ('s', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleDoubleOption ('\n', "accuracy", "set the accuracy that defines a non-zero projection on a given L subspace", 1e-10);
  (*OutputGroup) += new SingleStringOption ('o', "output-prefix", "prefix to use for the output files");

  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereGetAngularMomentumDecomposition -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int TotalLz = ((SingleIntegerOption*) Manager["total-lz"])->GetInteger();
  bool FermionFlag = false;
  bool HaldaneBasisFlag = ((BooleanOption*) Manager["haldane"])->GetBoolean();
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  double Accuracy = Manager.GetDouble("accuracy");
  if (((SingleStringOption*) Manager["statistics"])->GetString() == 0)
    FermionFlag = true;
  if (NbrParticles == 0)
    if (FQHEOnSphereFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["state"])->GetString(), NbrParticles, LzMax, TotalLz, FermionFlag) == false)
      {
	return -1;
      }


  long MemorySpace = 9l << 20;
  ParticleOnSphere* Space=0;
  if (FermionFlag == true)
    {
      if (HaldaneBasisFlag == false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	    if (LzMax <= 30)
#endif
	      if ((SymmetrizedBasis == false) || (TotalLz != 0))
		Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax, MemorySpace);
	      else
		Space = new FermionOnSphereSymmetricBasis(NbrParticles, LzMax, MemorySpace);
	    else
#ifdef __128_BIT_LONGLONG__
	      if (LzMax <= 126)
#else
		if (LzMax <= 62)
#endif
		  {
		    Space = new FermionOnSphereLong(NbrParticles, TotalLz, LzMax, MemorySpace);
		  }
		else
		  Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, LzMax, MemorySpace);
	}
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceState) == false)
	    return -1;
	  if (SymmetrizedBasis == false)
	    {
	      if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		Space = new FermionOnSphereHaldaneBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
	      else
		Space = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
	      if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
		{
		  ((FermionOnSphereHaldaneBasis*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
		  return 0;
		}
	    }
	  else
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
	}
    }
  else
    {
#ifdef  __64_BITS__
      if ((LzMax + NbrParticles - 1) < 63)
#else
	if ((LzMax + NbrParticles - 1) < 31)	
#endif
	  {
	    if ((SymmetrizedBasis == false) || (TotalLz != 0))
	      Space = new BosonOnSphereShort(NbrParticles, TotalLz, LzMax);
	    else
	      Space = new BosonOnSphereSymmetricBasisShort(NbrParticles, LzMax);
	  }
	else
	  {
	    if ((SymmetrizedBasis == false) || (TotalLz != 0))
	      Space = new BosonOnSphere (NbrParticles, TotalLz, LzMax);
	    else
	      Space = new BosonOnSphereSymmetricBasis(NbrParticles, LzMax);
	  }
    }
    
  
  RealVector State;
  if (State.ReadVector(Manager.GetString("state")) == false)
    {
      cout << "error while reading " << Manager.GetString("state") << endl;
      return -1;
    }
  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }

  int TotalMaxLz = LzMax * NbrParticles;
  if (FermionFlag == true)
    {
      TotalMaxLz = (LzMax - NbrParticles + 1) * NbrParticles;
    }
  cout << LzMax << " " << NbrParticles << " " << TotalMaxLz << endl;
  cout << "dim = " << Space->GetHilbertSpaceDimension() << endl;
  //  ParticleOnSphereL2Hamiltonian Hamiltonian (Space, NbrParticles, LzMax, TotalLz, Architecture.GetArchitecture(), 1.0, 0, true);

  RealVector TmpState(Space->GetHilbertSpaceDimension());
  RealVector TmpState2(Space->GetHilbertSpaceDimension());
  int l = 0;
  if (((NbrParticles * LzMax) & 1) != 0)
    l = 1;
  int* LValues = new int [TotalMaxLz + 1];
  double* Components = new double [TotalMaxLz + 1];
  char**  VectorFileNames = new char* [TotalMaxLz + 1];
  int NbrLValues = 0;
  //  RealVector InitialState;
  //  InitialState.Copy(State);
  double Sum = 0.0;
  cout.precision(14);

  for (; (l <= TotalMaxLz) && (fabs(1.0 - Sum) > Accuracy); l += 2)
    {
      TmpState.Copy(State);
      double TmpNorm = TmpState.Norm();
      cout << "truc " << TmpNorm << endl;
      for (int CurrentL = l + 2; (CurrentL <= TotalMaxLz) && (TmpNorm != Accuracy); CurrentL +=2)
	{
<<<<<<< .mine
	  ParticleOnSphereL2Hamiltonian Hamiltonian (Space, NbrParticles, LzMax, TotalLz, Architecture.GetArchitecture(), 1.0, 0, true);// / (-0.25 * ((double) (CurrentL * (CurrentL + 2)))), 0, true);
	  //	  Hamiltonian.ShiftHamiltonian(1.0);
	  //	  Hamiltonian.ShiftHamiltonian(-0.25 * (CurrentL * (CurrentL + 2)));
=======
//	  TmpState /= TmpNorm;
//	  Hamiltonian.ShiftHamiltonian(- 0.25 * (CurrentL * (CurrentL + 2)));
>>>>>>> .r1123
	  VectorHamiltonianMultiplyOperation Operation (&Hamiltonian, &TmpState, &TmpState2);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
<<<<<<< .mine
	  //cout << (TmpState * TmpState2) << endl;
	  //	  TmpNorm = TmpState2.Norm();
=======
//	  TmpNorm = TmpState2.Norm();
>>>>>>> .r1123
	  //	  TmpState2 /=  TmpNorm;
<<<<<<< .mine
	  //	  cout << CurrentL << " " << TmpNorm << endl;
	  //	  TmpState2 /= -0.25 * ((double) (CurrentL * (CurrentL + 2)));
=======
//	  cout << CurrentL << " " << TmpNorm << endl;
	  TmpState2 /= -0.25 * ((double) (CurrentL * (CurrentL + 2)));
//	  TmpState2 += TmpState;
>>>>>>> .r1123
	  RealVector TmpVector = TmpState2;
	  TmpState2 = TmpState;
	  TmpState = TmpVector;
	  TmpNorm = TmpState.Norm();
	  cout << CurrentL << " " << TmpNorm << endl;
	}
      if (TmpNorm > Accuracy)
	{
	  cout << "norm = " << TmpNorm << endl;
	  cout << Manager.GetString("state") << " has a non-zero projection on the L=" << (0.5 * ((double) l)) << " subspace" << endl;
	  //	  TmpState /= TmpNorm;
	  ParticleOnSphereL2Hamiltonian Hamiltonian (Space, NbrParticles, LzMax, TotalLz, Architecture.GetArchitecture(), 1.0, 0, true);
	  Hamiltonian.ShiftHamiltonian(0.0);
	  VectorHamiltonianMultiplyOperation Operation (&Hamiltonian, &TmpState, &TmpState2);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  cout << (TmpState * TmpState) << " " << (TmpState * TmpState2) << endl;
	  char* Extension = new char[32];
	  sprintf (Extension, "l_%d.vec", l);
	  if (Manager.GetString("output-prefix") != 0)
	    {
	      VectorFileNames[NbrLValues] = AddExtensionToFileName(Manager.GetString("output-prefix"), Extension);
	    }
	  else
	    {
	      VectorFileNames[NbrLValues] = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", Extension);
	    }
	  LValues[NbrLValues] = l;
	  Components[NbrLValues] = (TmpState * State);
	  Sum += Components[NbrLValues] * Components[NbrLValues];
	  cout << Components[NbrLValues] << " " << (1.0 - Sum) << endl;
	  TmpState.WriteVector(VectorFileNames[NbrLValues]);	  
	  delete[] Extension;
	  double Tmp = TmpState * State;	  
	  TmpState *= Tmp;
	  State -= TmpState;
	  ++NbrLValues;
	}
    }
  for (int i = 0; i < NbrLValues; ++i)
    {
      cout << LValues[i] << " " << Components[i] << " " << VectorFileNames[i] << endl;
    }
  cout << "sum = " << Sum << endl;
}
