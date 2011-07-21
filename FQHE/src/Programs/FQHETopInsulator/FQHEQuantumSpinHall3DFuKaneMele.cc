#include "Options/Options.h"

#include "HilbertSpace/FermionOnCubicLatticeWithSpinMomentumSpace.h"

#include "Hamiltonian/ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian.h"

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


// compute the single particle spectrum 
//
// outputFileName = name of the output file
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// nnHoping = nearest neighbor hoping amplitude
// nnnHoping =  next nearest neighbor hoping amplitude
// nnnnHoping =  second next nearest neighbor hoping amplitude
// mus = sublattice staggered chemical potential 
// mixingTermNorm = norm of the mixing term coupling the two copies of the checkerboard lattice
// mixingTermArgv = argument of the mixing term coupling the two copies of the checkerboard lattice
void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, int nbrSitesZ, double nnHoping, double nnnHoping, double nnnnHoping, double mus, double mixingTermNorm, double mixingTermArg);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEQuantumSpinHall3DFuKaneMele" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('z', "nbr-sitez", "number of sites along the z direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kz", "only evalute a given y momentum sector (negative if all kz sectors have to be computed)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive on-site potential strength between opposite spins", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive nearest neighbor potential strength between opposite spins", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mixing-norm", "norm of the mixing term coupling the two copies of the checkerboard lattice", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mixing-arg", "argument of the mixing term coupling the two copies of the checkerboard lattice (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption ('\n', "flat-band", "use flat band model");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEQuantumSpinHall3DFuKaneMele -h" << endl;
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
  int NbrSitesZ = Manager.GetInteger("nbr-sitez"); 
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky kz ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetBoolean("flat-band") == true)
    {
      if (Manager.GetDouble("mu-s") == 0.0)
	sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_v_%f_w_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
      else
	sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_v_%f_w_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f_mus_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
    }
  else
    {
      if (Manager.GetDouble("mu-s") == 0.0)
	sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_u_%f_v_%f_w_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
      else
	sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_u_%f_v_%f_w_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f_mus_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
    }

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"), Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg") * 2.0 * M_PI);
      return 0;
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
  int MinKz = 0;
  int MaxKz = NbrSitesZ - 1;
  if (Manager.GetInteger("only-kz") >= 0)
    {						
      MinKz = Manager.GetInteger("only-kz");
      MaxKz = MinKz;
    }
  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  for (int k = MinKz; k <= MaxKz; ++k)
	    {
	      cout << "(kx=" << i << ",ky=" << j << ",kz=" << k << ") " << endl;
	      FermionOnCubicLatticeWithSpinMomentumSpace Space(NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
	      cout << "dim = " << Space.GetHilbertSpaceDimension()  << endl;
	      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      Hamiltonian = new ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian(&Space, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ,
											       Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
											       Manager.GetDouble("tpp"), Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg") * 2.0 * M_PI, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
											       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d %d", i, j, k);
	      char* EigenstateOutputFile = new char [512];
	      if (Manager.GetBoolean("flat-band") == true)
		{
		  sprintf (EigenstateOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_t1_%f_t2_%f_gx_%f_gy_%f_kx_%d_ky_%d_kz_%d", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, 
			   Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j, k);
		}
	      else
		{
		  sprintf (EigenstateOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_u_%f_t1_%f_t2_%f_gx_%f_gy_%f_kx_%d_ky_%d_kz_%d", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, 
			   Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j, k);
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
    }
  return 0;
}

// compute the single particle spectrum 
//
// outputFileName = name of the output file
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// nnHoping = nearest neighbor hoping amplitude
// nnnHoping =  next nearest neighbor hoping amplitude
// nnnnHoping =  second next nearest neighbor hoping amplitude
// mus = sublattice staggered chemical potential 
// mixingTermNorm = norm of the mixing term coupling the two copies of the checkerboard lattice
// mixingTermArgv = argument of the mixing term coupling the two copies of the checkerboard lattice

void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, int nbrSitesZ, double nnHoping, double nnnHoping, double nnnnHoping, double mus, double mixingTermNorm, double mixingTermArg)
{
  ofstream File;
  File.open(outputFileName);
  File << "# kx    ky     kz     E_{-,1}    E_{-,2}    E_{+,1}    E_{+,2}" << endl;
  double MinEMinus = 0.0;
  double MaxEMinus = -10.0;
  double MinEPlus = 10.0;
  double MaxEPlus = 0.0;
  Complex MixingTerm = mixingTermNorm * Phase(mixingTermArg);
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  for (int kz = 0; kz < nbrSitesZ; ++kz)
	    {
	      HermitianMatrix TmpOneBodyHamiltonian(4, true);
	      Complex B1 = 4.0 * nnHoping * Complex (cos (1.0 * M_PI * ((double) kx) / ((double) nbrSitesX)) * cos (1.0 * M_PI * ((double) ky) / ((double) nbrSitesY)) * cos(M_PI * 0.25), 
						     sin (1.0 * M_PI * ((double) kx) / ((double) nbrSitesX)) * sin (1.0 * M_PI * ((double) ky) / ((double) nbrSitesY)) * sin(M_PI * 0.25));
	      double d1 = 4.0 * nnnnHoping * cos (2.0 * M_PI * ((double) kx) / ((double) nbrSitesX)) * cos (2.0 * M_PI * ((double) ky) / ((double) nbrSitesY));
	      double d3 = mus + (2.0 * nnnHoping * (cos (2.0 * M_PI * ((double) kx) / ((double) nbrSitesX))
						    - cos (2.0 * M_PI * ((double) ky) / ((double) nbrSitesY))));
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	      B1 = 4.0 * nnHoping * Complex (cos (1.0 * M_PI * ((double) -kx) / ((double) nbrSitesX)) * cos (1.0 * M_PI * ((double) -ky) / ((double) nbrSitesY)) * cos(M_PI * 0.25), 
					     sin (1.0 * M_PI * ((double) -kx) / ((double) nbrSitesX)) * sin (1.0 * M_PI * ((double) -ky) / ((double) nbrSitesY)) * sin(M_PI * 0.25));
	      d1 = 4.0 * nnnnHoping * cos (2.0 * M_PI * ((double) -kx) / ((double) nbrSitesX)) * cos (2.0 * M_PI * ((double) -ky) / ((double) nbrSitesY));
	      d3 = mus + (2.0 * nnnHoping * (cos (2.0 * M_PI * ((double) -kx) / ((double) nbrSitesX))
					     - cos (2.0 * M_PI * ((double) -ky) / ((double) nbrSitesY))));
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 2, d1 + d3);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 3, Conj(B1));
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 3, d1 - d3);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 3, - I() * MixingTerm);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 2, I() * MixingTerm);
	      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	      TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
	      TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif   
	      if (MaxEMinus < TmpDiag(0, 0))
		{
		  MaxEMinus = TmpDiag(0, 0);
		}
	      if (MinEMinus > TmpDiag(0, 0))
		{
		  MinEMinus = TmpDiag(0, 0);
		}
	      if (MaxEPlus < TmpDiag(2, 2))
		{
		  MaxEPlus = TmpDiag(2, 2);
		}
	      if (MinEPlus > TmpDiag(2, 2))
		{
		  MinEPlus = TmpDiag(2, 2);
		}
	      File << (2.0 * M_PI * ((double) kx) / ((double) nbrSitesX)) << " " << (2.0 * M_PI * ((double) ky) / ((double) nbrSitesY)) << " " << TmpDiag(0, 0) << " " << TmpDiag(1, 1) <<  " " << TmpDiag(2, 2) << " " << TmpDiag(3, 3) << endl;
	    }
	  File << endl;
	}
    }
  cout << "Spread = " << (MaxEMinus - MinEMinus) << "  Gap = " <<  (MinEPlus - MaxEMinus) << "  Flatening = " << ((MaxEMinus - MinEMinus) / (MinEPlus - MaxEMinus)) << endl;
}

