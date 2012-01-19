#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

//#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian.h"
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
// nbrSitesY = number of sites in the x direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda2 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice staggered chemical potential 
// mixingTermNorm = norm of the mixing term coupling the two copies of the kagome lattice
// mixingTermArgv = argument of the mixing term coupling the two copies of the kagome lattice
void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, double t1, double t2, double lambda1, double lambda2, double mus, double mixingTermNorm, double mixingTermArg);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEQuantumSpinHallKagomeModelTwoBands" , "0.01");
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
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength between identical spins", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive on-site potential strength between opposite spins", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive nearest neighbor potential strength between opposite spins", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "real part of the nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "real part of the next nearest neighbor hopping amplitude", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l1", "imaginary part of the nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l2", "imaginary part of the next nearest neighbor hopping amplitude", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice chemical potential on A site", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mixing-norm", "norm of the mixing term coupling the two copies of the kagome lattice", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mixing-arg", "argument of the mixing term coupling the two copies of the kagome lattice (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new BooleanOption ('\n', "decoupled", "assume two decoupled copies of the kagome lattice");
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
      cout << "see man page for option syntax or type FQHEQuantumSpinHallKagomeModelTwoBands -h" << endl;
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
  sprintf (CommentLine, "eigenvalues\n# kx ky E");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetBoolean("flat-band") == true)
    {
      if (Manager.GetDouble("mu-s") == 0.0)
	sprintf (EigenvalueOutputFile, "fermions_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
      else
	sprintf (EigenvalueOutputFile, "fermions_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f_mus_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
    }
  else
    {
      if (Manager.GetDouble("mu-s") == 0.0)
	sprintf (EigenvalueOutputFile, "fermions_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_u_%f_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
      else
	sprintf (EigenvalueOutputFile, "fermions_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_u_%f_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f_mus_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
    }

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg") * 2.0 * M_PI);
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
  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  if (Manager.GetBoolean("decoupled") == false)
	    {
	      cout << "(kx=" << i << ",ky=" << j << ") " << endl;
	      FermionOnSquareLatticeWithSpinMomentumSpace Space(NbrParticles, NbrSitesX, NbrSitesY, i, j);
	      cout << "dim = " << Space.GetHilbertSpaceDimension()  << endl;
	      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
	      AbstractQHEHamiltonian* Hamiltonian = 0;
// new ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian(&Space, NbrParticles, NbrSitesX, NbrSitesY,
// 											       Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
// 											       Manager.GetDouble("tpp"), Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg") * 2.0 * M_PI, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
// 											       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", i, j);
	      char* EigenstateOutputFile = new char [512];
	      if (Manager.GetBoolean("flat-band") == true)
		{
		  sprintf (EigenstateOutputFile, "fermions_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f_kx_%d_ky_%d", NbrParticles, NbrSitesX, NbrSitesY, 
			   Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
		}
	      else
		{
		  sprintf (EigenstateOutputFile, "fermions_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_u_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f_kx_%d_ky_%d", NbrParticles, NbrSitesX, NbrSitesY, 
			   Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
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
	  else
	    {
	      for (int Sz = -NbrParticles; Sz <= NbrParticles; Sz += 2)
		{
		  cout << "(kx=" << i << ",ky=" << j << ") Sz=" << Sz << " : " << endl;
		  FermionOnSquareLatticeWithSpinMomentumSpace Space(NbrParticles, (Sz + NbrParticles) / 2, NbrSitesX, NbrSitesY, i, j);
		  cout << "dim = " << Space.GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
		  AbstractQHEHamiltonian* Hamiltonian = 0;
		  Hamiltonian = new ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian(&Space, NbrParticles, NbrSitesX, NbrSitesY,
												      Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
												      Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"),
												      Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
												      Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  char* ContentPrefix = new char[256];
		  sprintf (ContentPrefix, "%d %d %d", i, j, Sz);
		  char* EigenstateOutputFile = new char [512];
		  if (Manager.GetBoolean("flat-band") == true)
		    {
		      if (Manager.GetDouble("mu-s") == 0.0)
			sprintf (EigenstateOutputFile, "fermions_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f_kx_%d_ky_%d_sz_%d", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
				 Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j, Sz);
		      else
			sprintf (EigenstateOutputFile, "fermions_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d_sz_%d", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
				 Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"), i, j, Sz);
		    }
		  else
		    {
		      if (Manager.GetDouble("mu-s") == 0.0)
			sprintf (EigenstateOutputFile, "fermions_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_u_%f_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f_kx_%d_ky_%d_sz_%d", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), 
				 Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j, Sz);
		      else
			sprintf (EigenstateOutputFile, "fermions_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_u_%f_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d_sz_%d", NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), 
				 Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"), i, j, Sz);
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
    }
  return 0;
}

// compute the single particle spectrum 
//
// outputFileName = name of the output file
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the x direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda2 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice staggered chemical potential 
// mixingTermNorm = norm of the mixing term coupling the two copies of the kagome lattice
// mixingTermArgv = argument of the mixing term coupling the two copies of the kagome lattice

void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, double t1, double t2, double lambda1, double lambda2, double mus, double mixingTermNorm, double mixingTermArg)
{
  ofstream File;
  File.open(outputFileName);
  File << "# kx    ky     E_{-,1}    E_{-,2}    E_{+,1}    E_{+,2}" << endl;
  double MinEMinus = 0.0;
  double MaxEMinus = -10.0;
  double MinEPlus = 10.0;
  double MaxEPlus = 0.0;
  Complex MixingTerm = mixingTermNorm * Phase(mixingTermArg);
  double NNHopping = t1;
  double NNSpinOrbit = lambda1;
  double NextNNHopping = t2;
  double NextNNSpinOrbit = lambda2;
  double KxFactor =  2.0 * M_PI / ((double) nbrSitesX);
  double KyFactor = 2.0 * M_PI / ((double) nbrSitesY);

  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  HermitianMatrix TmpOneBodyHamiltonian(6, true);
	  double KX = 0.5 * KxFactor * ((double) kx);
	  double KY = 0.5 * KyFactor * ((double) ky);
	  Complex HAB (-2.0 * NNHopping, -2.0 * NNSpinOrbit);
	  HAB *= cos (KX);
	  Complex HAC(-2.0 * NNHopping, 2.0 * NNSpinOrbit);
	  HAC *= cos (KY);
	  Complex HBC(-2.0 * NNHopping, -2.0 * NNSpinOrbit);
	  HBC *= cos(KX - KY);
	  
	  Complex HAB2 (-2.0 * NextNNHopping, 2.0 * NextNNSpinOrbit);
	  HAB2 *= cos (KX - 2.0 * KY);
	  Complex HAC2 (-2.0 * NextNNHopping, -2.0 * NextNNSpinOrbit);
	  HAC2 *= cos (2.0 * KX - KY);
	  Complex HBC2 (-2.0 * NextNNHopping, 2.0  *  NextNNSpinOrbit);
	  HBC2 *= cos (KX + KY);
	  
	  HAB += HAB2;
	  HAC += HAC2;
	  HBC += HBC2;

	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HAB);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HAC);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HBC);

	  KX = 0.5 * KxFactor * ((double) -kx);
	  KY = 0.5 * KyFactor * ((double) -ky);
	  HAB = Complex(-2.0 * NNHopping, -2.0 * NNSpinOrbit);
	  HAB *= cos (KX);
	  HAC = Complex(-2.0 * NNHopping, 2.0 * NNSpinOrbit);
	  HAC *= cos (KY);
	  HBC = Complex(-2.0 * NNHopping, -2.0 * NNSpinOrbit);
	  HBC *= cos(KX - KY);
	  
	  HAB2 = Complex(-2.0 * NextNNHopping, 2.0 * NextNNSpinOrbit);
	  HAB2 *= cos (KX - 2.0 * KY);
	  HAC2 = Complex(-2.0 * NextNNHopping, -2.0 * NextNNSpinOrbit);
	  HAC2 *= cos (2.0 * KX - KY);
	  HBC2 = Complex(-2.0 * NextNNHopping, 2.0 * NextNNSpinOrbit);
	  HBC2 *= cos (KX + KY);
	  
	  HAB += HAB2;
	  HAC += HAC2;
	  HBC += HBC2;
	  
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 4, Conj(HAB));
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 5, Conj(HAC));
	  TmpOneBodyHamiltonian.SetMatrixElement(4, 5, Conj(HBC));
	  
// 	  TmpOneBodyHamiltonian.SetMatrixElement(0, 3, - I() * MixingTerm);
// 	  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, I() * MixingTerm);

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
  cout << "Spread = " << (MaxEMinus - MinEMinus) << "  Gap = " <<  (MinEPlus - MaxEMinus) << "  Flatening = " << ((MaxEMinus - MinEMinus) / (MinEPlus - MaxEMinus)) << endl;
}

