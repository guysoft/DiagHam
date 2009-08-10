#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Operator/ParticleOnSphereWithSpinLMinusOperator.h"

#include "Tools/FQHESpectrum/QHEOnSphereLzSortedSpectrum.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
//#include "HilbertSpace/FermionOnSphereUnlimited.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereWithSpinLMinus" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += DataGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "sz", "twice the total sz value of the system for the initial state", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "twice the total lz value of the system for the initial state", 0);
  (*SystemGroup) += new SingleStringOption  ('S', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lm", "number of time the L- operator has to be applied", 1);

  (*DataGroup) += new SingleStringOption  ('i', "input-file", "input vector file name");
  (*DataGroup) += new SingleStringOption  ('o', "output-file", "output vector file name");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type LMinus -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int Lz = Manager.GetInteger("lz");
  int TotalSz = Manager.GetInteger("sz");
  int NbrLMinus = Manager.GetInteger("nbr-lm");
  bool FermionFlag = false;
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-file"), NbrParticles, LzMax, Lz, TotalSz, FermionFlag) == false)
    {
      return -1;
    }
  if (Manager.GetString("statistics") != 0)
    if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
      {
	FermionFlag = true;
      }
    else
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = false;
	}
      else
	{
	  cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	}  
  int Parity = Lz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }
  if (NbrParticles&1 != TotalSz&1)
    {
      cout << "Sz and NbrParticles must have the parity" << endl;
      return -1;
    }

  RealVector InitialVector; 
  RealVector TargetVector; 
  if (InitialVector.ReadVector(Manager.GetString("input-file")) == false)
    {
      cout << "error while reading " << Manager.GetString("input-file") << endl;
      return -1;
    }
	
  long MemorySpace = 9l << 20;
  ParticleOnSphereWithSpin* IntialSpace;
  ParticleOnSphereWithSpin* TargetSpace;
  if (FermionFlag == true)
    {
#ifdef __64_BITS__
      if (LzMax <= 31)
	{
	  IntialSpace = new FermionOnSphereWithSpin(NbrParticles, Lz, LzMax, TotalSz, MemorySpace);
	}
      else
	{
	  // IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	  cout << "Fermions with Spin not defined yet for LzMax > 31"<<endl;
	  exit(-1);
	}
#else
      if (LzMax <= 15)
	{
	  IntialSpace = new FermionOnSphereWithSpin(NbrParticles, Lz, LzMax, TotalSz, MemorySpace);
	}
      else
	{
	  // IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	  cout << "Fermions with Spin not defined yet for LzMax > 15, consider using a 64 bit machine!"<<endl;
	  exit(-1);
	}
#endif
    }
  else
    {
      IntialSpace = new BosonOnSphereWithSpin(NbrParticles, Lz, LzMax, TotalSz);
    }
  for (int i = 1; i <= NbrLMinus; ++i)
    {
      if (FermionFlag == true)
	{
#ifdef __64_BITS__
	  if (LzMax <= 31)
	    {
	      TargetSpace = new FermionOnSphereWithSpin(NbrParticles, Lz, LzMax - (2 * i), TotalSz, MemorySpace);
	    }
	  else
	    {
	      // TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax - (2 * i), TotalSz, MemorySpace);
	      cout << "Fermions with Spin not defined yet for LzMax > 31"<<endl;
	      exit(-1);
	    }
#else
	  if (LzMax <= 15)
	    {
	      TargetSpace = new FermionOnSphereWithSpin(NbrParticles, Lz, LzMax - (2 * i), TotalSz, MemorySpace);
	    }
	  else
	    {
	      // TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax - (2 * i), TotalSz, MemorySpace);
	      cout << "Fermions with Spin not defined yet for LzMax > 15, consider using a 64 bit machine!"<<endl;
	      exit(-1);
	    }
#endif
	}
      else
	{
	  TargetSpace = new BosonOnSphereWithSpin(NbrParticles, Lz, LzMax - (2 * i), TotalSz);
	}
      IntialSpace->SetTargetSpace(TargetSpace);
      TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension());
      ParticleOnSphereWithSpinLMinusOperator LMinus(IntialSpace, Lz, LzMax  - (2 * i) + 2);
      LMinus.Multiply(InitialVector, TargetVector);
      delete IntialSpace;
      IntialSpace = TargetSpace;
      InitialVector = TargetVector;
    }  
  if (InitialVector.WriteVector(Manager.GetString("output-file")) == false)
    {
      cout << "error while writing " << Manager.GetString("output-file") << endl;
      return -1;
    }
  delete IntialSpace;
  return 0;
}

