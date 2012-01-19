#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorHamiltonian.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/GenericComplexMainTask.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// evaluate Hilbert space dimension for fermions
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension
long EvaluateHilbertSpaceDimension(int nbrFermions, int currentPx, int currentPy, int currentTotalPx, int currentTotalPy, int pxConstraint, int pyConstraint);

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereGetDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "Hubbard potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "band-parameter", "band structure parameter", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "kinetic-factor", "multiplicative factor in front of the kinetic term", 1.0);
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereGetDimension -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n#");
  char* EigenvalueOutputFile = new char [512];
  sprintf (EigenvalueOutputFile, "fermions_cherninsulator_n_%d_x_%d_y_%d_u_%f_m_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("band-parameter"));

  int MinKx = 0;
  int MaxKx = NbrSitesX - 1;
  if (Manager.GetBoolean("full-momentum") == false)
    {
      MaxKx = NbrSitesX / 2;
    }
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrSitesY - 1;
  if (Manager.GetBoolean("full-momentum") == false)
    {
      MaxKy = NbrSitesY / 2;
    }
  if (Manager.GetInteger("only-ky") >= 0)
    {						
      MinKy = Manager.GetInteger("only-ky");
      MaxKy = MinKy;
    }
  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
	  FermionOnSquareLatticeWithSpinMomentumSpace Space(NbrParticles, NbrSitesX, NbrSitesY, i, j);
 	  cout << "dim = " << Space.GetHilbertSpaceDimension()  << endl;
	  Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
	  AbstractQHEHamiltonian* Hamiltonian = new ParticleOnLatticeWithSpinChernInsulatorHamiltonian(&Space, NbrParticles, NbrSitesX, NbrSitesY,
												       Manager.GetDouble("kinetic-factor"), Manager.GetDouble("u-potential"), Manager.GetDouble("band-parameter"),
												       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d %d", i, j);
	  char* EigenstateOutputFile = new char [512];
	  sprintf (EigenstateOutputFile, "fermions_cherninsulator_n_%d_x_%d_y_%d_u_%f_m_%f_kx_%d_ky_%d", NbrParticles, NbrSitesX, NbrSitesY, 
		   Manager.GetDouble("u-potential"), Manager.GetDouble("band-parameter"), i, j);
	  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	  FirstRunFlag = false;
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  delete Hamiltonian;
	  delete[] EigenstateOutputFile;
	  delete[] ContentPrefix;
	}
    }
  return 0;
}

// evaluate Hilbert space dimension for fermions
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long EvaluateHilbertSpaceDimension(int nbrFermions, int currentPx, int currentPy, int currentTotalPx, int currentTotalPy, int pxConstraint, int pyConstraint)
{
  if (currentPy >= pyConstraint)
    {
      currentPy = 0;
      currentPx++;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalPx % pxConstraint) == 0) && ((currentTotalPy % pyConstraint) == 0))
	return 1l;
      else	
	return 0l;
    }
  if (currentPx >= pxConstraint)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int i = currentPx; i < pxConstraint; ++i)
	{
	  for (int j = currentPy; j < pyConstraint; ++j)
	    {
	      if ((((i + currentTotalPx) % pxConstraint) == 0) && (((j + currentTotalPy) % pyConstraint) == 0))
		++Count;
	    }
	}
      return Count;
    }
  Count += EvaluateHilbertSpaceDimension(nbrFermions - 2, currentPx, currentPy + 1, currentTotalPx + (2 * currentPx), currentTotalPy + (2 * currentPy), pxConstraint, pyConstraint);
  Count += EvaluateHilbertSpaceDimension(nbrFermions - 1, currentPx, currentPy + 1, currentTotalPx + currentPx, currentTotalPy + currentPy, pxConstraint, pyConstraint);
  Count += EvaluateHilbertSpaceDimension(nbrFermions, currentPx, currentPy + 1, currentTotalPx, currentTotalPy, pxConstraint, pyConstraint);
  return Count;
}
