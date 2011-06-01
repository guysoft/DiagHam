#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandCheckerboardHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

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


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEQuantumSpinHallCheckerboardModelTwoBands" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive nearest next neighbor potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mixing-norm", "norm of the mixing term coupling the two copies of the checkerboard lattice", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mixing-arg", "argument of the mixing term coupling the two copies of the checkerboard lattice (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new BooleanOption ('\n', "decoupled", "assume two decoupled copies of the checkerboard lattice");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEQuantumSpinHallCheckerboardModelTwoBands -h" << endl;
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
  if (Manager.GetBoolean("flat-band") == true)
    {
      sprintf (EigenvalueOutputFile, "fermions_twoband_quantumspinhall_checkerboardlattice_n_%d_x_%d_y_%d_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
    }
  else
    {
      sprintf (EigenvalueOutputFile, "fermions_twoband_quantumspinhall_checkerboardlattice_n_%d_x_%d_y_%d_u_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
    }

  int MinKx = 0;
  int MaxKx = NbrSitesX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrSitesY - 1;
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
 	  AbstractQHEHamiltonian* Hamiltonian = 0;
	  if (Manager.GetBoolean("decoupled") == false)
	    {
	      Hamiltonian = new ParticleOnLatticeQuantumSpinHallTwoBandCheckerboardHamiltonian(&Space, NbrParticles, NbrSitesX, NbrSitesY,
											       Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
											       Manager.GetDouble("tpp"), Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg") * 2.0 * M_PI, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
											       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	    }
	  else
	    {
	      Hamiltonian = new ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian(&Space, NbrParticles, NbrSitesX, NbrSitesY,
													Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
													Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
													Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	    }
	      char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d %d", i, j);
	  char* EigenstateOutputFile = new char [512];
	  if (Manager.GetBoolean("flat-band") == true)
	    {
	      sprintf (EigenstateOutputFile, "fermions_twoband_quantumspinhall_checkerboardlattice_n_%d_x_%d_y_%d_t1_%f_t2_%f_gx_%f_gy_%f_kx_%d_ky_%d", NbrParticles, NbrSitesX, NbrSitesY, 
		       Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
	    }
	  else
	    {
	      sprintf (EigenstateOutputFile, "fermions_twoband_quantumspinhall_checkerboardlattice_n_%d_x_%d_y_%d_u_%f_t1_%f_t2_%f_gx_%f_gy_%f_kx_%d_ky_%d", NbrParticles, NbrSitesX, NbrSitesY, 
		       Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
	    }
	  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	  FirstRunFlag = false;
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  cout << "------------------------------------" << endl;
	  delete Hamiltonian;
	  delete[] EigenstateOutputFile;
	  delete[] ContentPrefix;
	}
    }
  return 0;
}

