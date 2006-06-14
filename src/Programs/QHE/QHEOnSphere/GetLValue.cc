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

#include "Operator/QHEOperator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Tools/QHE/QHESpectrum/QHEOnSphereLzSortedSpectrum.h"

#include "HilbertSpace/QHEHilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnSphereUnlimited.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FindSystemInfoFromFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, bool& statistics);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("GetLValue" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += DataGroup;
  Manager += MiscGroup;
 
  (*DataGroup) += new SingleStringOption  ('\0', "state", "name of the file that contains the state whose average L value has to be evaluated");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleStringOption  ('s', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");
  (*DataGroup) += new SingleDoubleOption  ('\n', "l-error", "error above which a vector is no more considerated as an eigenvector of L^2", 1e-12);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHENBodyQuasiHoleOverlap -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "no input state" << endl << "see man page for option syntax or type QHENBodyQuasiHoleOverlap -h" << endl;
      return -1;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["lz"])->GetInteger();
  bool FermionFlag = false;
  if (((SingleStringOption*) Manager["statistics"])->GetString() == 0)
    FermionFlag = true;
  if (FindSystemInfoFromFileName(((SingleStringOption*) Manager["state"])->GetString(), NbrParticles, LzMax, Lz, FermionFlag) == false)
    {
      return -1;
    }
  if (((SingleStringOption*) Manager["statistics"])->GetString() != 0)
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
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  char* StateFileName = ((SingleStringOption*) Manager["state"])->GetString();
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


  long MemorySpace = 9l << 20;
  ParticleOnSphere* Space;
  if (FermionFlag == true)
    {
#ifdef __64_BITS__
      if (LzMax <= 63)
	{
	  Space = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  Space = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#else
      if (LzMax <= 31)
	{
	  Space = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  Space = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#endif
    }
  else
    {
      Space = new BosonOnSphere(NbrParticles, Lz, LzMax);
    }

  int TotalMaxLz = LzMax * NbrParticles;
  if (FermionFlag == true)
    {
      TotalMaxLz = (LzMax - NbrParticles + 1) * NbrParticles;
    }

   ParticleOnSphereSquareTotalMomentumOperator oper(Space, Lz, LzMax);
   double RawTmpAngularMomentum = 0.5 * (sqrt ((4.0 * oper.MatrixElement(State, State).Re) + 1.0) - 1.0);
   cout << "<L> = " << RawTmpAngularMomentum << endl;
   return 0;
}


// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured
 
bool FindSystemInfoFromFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, bool& statistics)
{
  char* StrNbrParticles;
  if (nbrParticles == 0)
    {
      StrNbrParticles = strstr(filename, "_n_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 3;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      nbrParticles = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess number of particles from file name " << filename << endl
	       << "use --nbr-particles option" << endl;
	  return false;            
	}
    }
  if (lzMax == 0)
    {
      StrNbrParticles = strstr(filename, "_2s_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 4;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      lzMax = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess maximum momentum from file name " << filename << endl
	       << "use --lzmax option" << endl;
	  return false;            
	}
    }
  if (lz == 0)
    {
      StrNbrParticles = strstr(filename, "_lz_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 4;
	  int SizeString = 0;
	  if (StrNbrParticles[SizeString] == '-')
	    ++SizeString;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '.') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      lz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '.';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the angular momentum from file name " << filename << endl
	       << "use --lz option" << endl;
	  return false;            
	}
    }
  if (statistics == true)
    {
      if (strstr(filename, "fermion") == 0)
	{
	  if (strstr(filename, "boson") == 0)
	    {
	      cout << "can't guess particle statistics from file name " << filename << endl
		   << "use --statistics option" << endl;
	      return false;	  
	    }
	  else
	    {
	      statistics = false;
	    }
	}
      else
	{
	  statistics = true;
	}
    }
  return true;
}
