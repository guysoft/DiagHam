#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Operator/ParticleOnSphereLMinusOperator.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("LMinus" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += DataGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "twice the total lz value of the system for the initial state", 0);
  (*SystemGroup) += new SingleStringOption  ('s', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lm", "number of time the L- operator has to be applied", 1);

  (*SystemGroup) += new SingleStringOption  ('\n', "input-reference", "use a haldane basis with the given reference file for the input file");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-reference", "use a haldane basis with the given reference file for the output file");


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

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["lz"])->GetInteger();
  int NbrLMinus = ((SingleIntegerOption*) Manager["nbr-lm"])->GetInteger();
  bool FermionFlag = false;
  if (((SingleStringOption*) Manager["statistics"])->GetString() == 0)
    FermionFlag = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["input-file"])->GetString(), NbrParticles, LzMax, Lz, FermionFlag) == false)
    {
      return -1;
    }
  if ((((SingleStringOption*) Manager["statistics"])->GetString()) != 0)
    if ((strcmp ("fermions", ((SingleStringOption*) Manager["statistics"])->GetString()) == 0))
      {
	FermionFlag = true;
      }
    else
      if ((strcmp ("fermions", ((SingleStringOption*) Manager["statistics"])->GetString()) == 0))
	{
	  FermionFlag = false;
	}
      else
	{
	  cout << ((SingleStringOption*) Manager["statistics"])->GetString() << " is an undefined statistics" << endl;
	}  
  int Parity = Lz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the same parity" << endl;
      return -1;           
    }

  RealVector InitialVector; 
  RealVector TargetVector; 
  if (InitialVector.ReadVector(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "error while reading " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;
    }
	
  long MemorySpace = 9l << 20;
  ParticleOnSphere* IntialSpace;
  ParticleOnSphere* TargetSpace;
  cout << "Creating space for Lz="<<Lz<<"/2"<<endl;
  if (FermionFlag == true)
    {
#ifdef __64_BITS__
      if (LzMax <= 63)
	{
	  IntialSpace = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#else
      if (LzMax <= 31)
	{
	  IntialSpace = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#endif
    }
  else
    {
      IntialSpace = new BosonOnSphereShort(NbrParticles, Lz, LzMax);
    }
  for (int i = 1; i <= NbrLMinus; ++i)
    {
//       if (Lz - (2 * i) < - LzMax)
// 	{
// 	  cout << "Cannot apply LMinus more than "<<i-1<<" times"<<endl;
// 	  exit(-1);
// 	}
      cout << "Creating space for Lz="<<Lz - (2 * i)<<"/2"<<endl;
      if (FermionFlag == true)
	{
#ifdef __64_BITS__
	  if (LzMax <= 63)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles, Lz - (2 * i), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz - (2 * i), LzMax, MemorySpace);
	    }
#else
	  if (LzMax <= 31)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles, Lz - (2 * i), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz - (2 * i), LzMax, MemorySpace);
	    }
#endif
	}
      else
	{
	  TargetSpace = new BosonOnSphereShort(NbrParticles, Lz - (2 * i), LzMax);
	}
      IntialSpace->SetTargetSpace(TargetSpace);
      TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension(), true);
      if (TargetSpace->GetHilbertSpaceDimension()!=IntialSpace->GetTargetHilbertSpaceDimension())
	{
	  cout << "Problem with setting target space"<<endl;
	  exit(-1);
	}
      ParticleOnSphereLMinusOperator LMinus(IntialSpace, Lz - (2 * i) + 2, LzMax);
      LMinus.Multiply(InitialVector, TargetVector);
      delete IntialSpace;
      IntialSpace = TargetSpace;
      InitialVector = TargetVector;
      InitialVector/=InitialVector.Norm();
    }
  char *OutputName;
  if (Manager.GetString("output-file")!=NULL)
    {
      OutputName = new char[strlen(Manager.GetString("output-file"))+1];
      strcpy(OutputName,Manager.GetString("output-file"));
    }
  else
    {
      OutputName = new char[strlen(Manager.GetString("input-file"))+10];
      sprintf(OutputName,"%s_L-",Manager.GetString("input-file"));
    }  
  if (InitialVector.WriteVector(OutputName) == false)
    {
      cout << "error while writing " << OutputName << endl;
      return -1;
    }
  delete IntialSpace;
  delete [] OutputName;
  return 0;
}

