#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
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

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEChernInsulatorSingleBand" , "0.01");
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
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive nearest next neighbor potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "band-parameter", "band structure parameter", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEChernInsulatorSingleBand -h" << endl;
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

  char* FilePrefix = new char [256];
  if (Manager.GetBoolean("single-band") == false)
    {
      sprintf (FilePrefix, "fermions_checkerboardlattice_n_%d_x_%d_y_%d",  NbrParticles, NbrSitesX, NbrSitesY);
    }
  else
    {
      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
	{ 
	  sprintf (FilePrefix, "fermions_singleband_checkerboardlattice_n_%d_x_%d_y_%d",  NbrParticles, NbrSitesX, NbrSitesY);
	}
      else
	{
	  if (Manager.GetBoolean("three-body") == true)
	    {
	      sprintf (FilePrefix, "fermions_singleband_threebody_checkerboardlattice_n_%d_x_%d_y_%d",  NbrParticles, NbrSitesX, NbrSitesY);
	    }
	  else
	    {
	      if (Manager.GetBoolean("four-body") == true)
		{
		  sprintf (FilePrefix, "fermions_singleband_fourbody_checkerboardlattice_n_%d_x_%d_y_%d",  NbrParticles, NbrSitesX, NbrSitesY);
		}
	      else
		{
		  sprintf (FilePrefix, "fermions_singleband_fivebody_checkerboardlattice_n_%d_x_%d_y_%d",  NbrParticles, NbrSitesX, NbrSitesY);
		}
	    }
	}
    }

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  sprintf (EigenvalueOutputFile, "fermions_singleband_twoorbitalcherninsulator_n_%d_x_%d_y_%d_m_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("band-parameter"));

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("mu-s"));
      return 0;
    }

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
	  FermionOnSquareLatticeMomentumSpace Space(NbrParticles, NbrSitesX, NbrSitesY, i, j);
 	  cout << "dim = " << Space.GetHilbertSpaceDimension()  << endl;
 	  Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
 	  AbstractQHEHamiltonian* Hamiltonian = new ParticleOnLatticeChernInsulatorSingleBandHamiltonian(&Space, NbrParticles, NbrSitesX, NbrSitesY, 
													 Manager.GetDouble("band-parameter"),
													 Architecture.GetArchitecture(), Memory);
 	  char* ContentPrefix = new char[256];
 	  sprintf (ContentPrefix, "%d %d", i, j);
	  char* EigenstateOutputFile = new char [512];
	  sprintf (EigenstateOutputFile, "fermions_singleband_twoorbitalcherninsulator_n_%d_x_%d_y_%d_m_%f_kx_%d_ky_%d", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("band-parameter"), i, j);
 	  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, 
				      FirstRunFlag, EigenstateOutputFile);
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

