#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"

#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandCheckerboardHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian.h"
#include "Hamiltonian/ExplicitHamiltonian.h"
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
// nnHoping = nearest neighbor hoping amplitude
// nnnHoping =  next nearest neighbor hoping amplitude
// nnnnHoping =  second next nearest neighbor hoping amplitude
// mus = sublattice staggered chemical potential 
// mixingTermNorm = norm of the mixing term coupling the two copies of the checkerboard lattice
// mixingTermArgv = argument of the mixing term coupling the two copies of the checkerboard lattice
void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, double nnHoping, double nnnHoping, double nnnnHoping, double mus, double mixingTermNorm, double mixingTermArg);

// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nnHoping = nearest neighbor hoping amplitude
// nnnHoping =  next nearest neighbor hoping amplitude
// nnnnHoping =  second next nearest neighbor hoping amplitude
// mus = sublattice staggered chemical potential 
// mixingTermNorm = norm of the mixing term coupling the two copies of the checkerboard lattice
// mixingTermArgv = argument of the mixing term coupling the two copies of the checkerboard lattice
ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, double nnHoping, double nnnHoping, double nnnnHoping, double mus, double mixingTermNorm, double mixingTermArg);


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
  (*SystemGroup) += new BooleanOption ('\n', "decoupled", "assume two decoupled copies of the checkerboard lattice");
  (*SystemGroup) += new BooleanOption ('\n', "four-bands", "perform the calculations within the full four band model");
  (*SystemGroup) += new BooleanOption ('\n', "project-fourbands", "project the hamiltonian from the four band model to the two band model");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "show matrix representation of the hamiltonian");
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
  int TotalNbrSites = NbrSitesX * NbrSitesY;
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [16];
//   if (Manager.GetBoolean("boson") == false)
//     {
      sprintf (StatisticPrefix, "fermions");
//     }
//   else
//     {
//       sprintf (StatisticPrefix, "bosons");
//     }

  char* FilePrefix = new char [512];
  if (Manager.GetBoolean("four-bands") == false)
    {
      sprintf (FilePrefix, "%s_twoband_quantumspinhall_checkerboardlattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
    }
  else
    {
      sprintf (FilePrefix, "%s_fourband_quantumspinhall_checkerboardlattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
    }
  char* FileBandParameters =  new char [512];
  if (Manager.GetDouble("mu-s") == 0.0)
    {
      if (Manager.GetBoolean("decoupled") == true)
	{
	  if (Manager.GetBoolean("flat-band") == true)
	    {
	      sprintf (FileBandParameters, "v_%f_w_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f", 
		       Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
		       Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
		       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	    }
	  else
	    {
	      sprintf (FileBandParameters, "u_%f_v_%f_w_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f", 
		       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
		       Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
		       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	    }
	}
      else
	{
	  if (Manager.GetBoolean("flat-band") == true)
	    {
	      sprintf (FileBandParameters, "v_%f_w_%f_t1_%f_t2_%f_tpp_%f_delta_%f_phi_%f_gx_%f_gy_%f", 
		       Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
		       Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
		       Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg"), 
		       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	    }
	  else
	    {
	      sprintf (FileBandParameters, "u_%f_v_%f_w_%f_t1_%f_t2_%f_tpp_%f_delta_%f_phi_%f_gx_%f_gy_%f", 
		       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
		       Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
		       Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg"), 
		       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("flat-band") == true)
	{
	  if (Manager.GetBoolean("decoupled") == true)
	    {
	      sprintf (FileBandParameters, "v_%f_w_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f_mus_%f", 
		       Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
		       Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
		       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	    }
	  else
	    {
	      sprintf (FileBandParameters, "v_%f_w_t1_%f_t2_%f_tpp_%f_delta_%f_phi_%f_gx_%f_gy_%f_mus_%f", 
		       Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
		       Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
		       Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg"), 
		       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	    }
	}
      else
	{
	  if (Manager.GetBoolean("decoupled") == true)
	    {
	      sprintf (FileBandParameters, "u_%f_v_%f_w_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f_mus_%f", 
		       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
		       Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
		       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	    }
	  else
	    {
	      sprintf (FileBandParameters, "u_%f_v_%f_w_t1_%f_t2_%f_tpp_%f_delta_%f_phi_%f_gx_%f_gy_%f_mus_%f", 
		       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
		       Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
		       Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg"), 
		       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	    }
	}
    }

  char* CommentLine = new char [256];
  if (Manager.GetBoolean("decoupled") == true)
    {
      sprintf (CommentLine, "eigenvalues\n# kx ky Sz ");
    }
  else
    {
      sprintf (CommentLine, "eigenvalues\n# kx ky ");
    }
  char* EigenvalueOutputFile = new char [512 + strlen(FilePrefix) + strlen(FileBandParameters)];
  sprintf (EigenvalueOutputFile, "%s_%s.dat", FilePrefix, FileBandParameters);

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"), Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg") * 2.0 * M_PI);
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
	      if (Manager.GetBoolean("four-bands") == false)
		{
		  FermionOnSquareLatticeWithSpinMomentumSpace Space(NbrParticles, NbrSitesX, NbrSitesY, i, j);
		  cout << "dim = " << Space.GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
		  AbstractQHEHamiltonian* Hamiltonian = 0;
		  Hamiltonian = new ParticleOnLatticeQuantumSpinHallTwoBandCheckerboardHamiltonian(&Space, NbrParticles, NbrSitesX, NbrSitesY,
												   Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
												   Manager.GetDouble("tpp"), Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg") * 2.0 * M_PI, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
												   Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  char* ContentPrefix = new char[256];
		  sprintf (ContentPrefix, "%d %d", i, j);
		  char* EigenstateOutputFile = new char [512 + strlen(FilePrefix) + strlen(FileBandParameters)];
		  sprintf (EigenstateOutputFile, "%s_%s_kx_%d_ky_%d", FilePrefix, FileBandParameters, i, j);
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
		  ParticleOnSphere* Space = 0;
#ifdef __128_BIT_LONGLONG__
		  if (TotalNbrSites <= 15)
#else
		    if (TotalNbrSites <= 7)
#endif
		      {
			Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		      }
		    else
		      {
			Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		      }
		  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
		  AbstractHamiltonian* Hamiltonian = 0;
		  Hamiltonian = new ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian ((ParticleOnSphereWithSU4Spin*) Space, NbrParticles, NbrSitesX, NbrSitesY,
												     Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
												     Manager.GetDouble("t1"), Manager.GetDouble("t2"),
												     Manager.GetDouble("tpp"), Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg") * 2.0 * M_PI, 
												     Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),  		     
												     Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  
		  if (Manager.GetBoolean("project-fourbands") == true)
		    {
		      ComplexMatrix* OneBodyBasis = ComputeSingleParticleTransformationMatrices(NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"),
												Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"), Manager.GetDouble("mixing-norm"), Manager.GetDouble("mixing-arg") * 2.0 * M_PI);
		      ComplexMatrix NBodyTransformationMatrix = ((FermionOnSquareLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixOneBodyBasis(OneBodyBasis);
		      ComplexMatrix HRep (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension());
		      Hamiltonian->GetHamiltonian(HRep);
		      ComplexMatrix TransformedHRep = HRep.Conjugate(NBodyTransformationMatrix);
		      FermionOnSquareLatticeWithSpinMomentumSpace* TargetSpace = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		      ComplexMatrix SU4SU2TransformationMatrix = ((FermionOnSquareLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixSU4ToSU2(TargetSpace, 2, 3);
		      ComplexMatrix TransformedHRep2 = TransformedHRep.InvConjugate(SU4SU2TransformationMatrix);
		      if (Manager.GetDouble("u-potential") != 0.0)
			TransformedHRep2 /= Manager.GetDouble("u-potential");
		      
		      RealDiagonalMatrix TmpDiag;
		      HermitianMatrix HRep2(TransformedHRep2);
		      delete Hamiltonian;
		      delete Space;
		      delete[] OneBodyBasis;
		      Hamiltonian = new ExplicitHamiltonian(TargetSpace, &HRep2);
		      Space = TargetSpace;
		    }
		  
		  char* ContentPrefix = new char[256];
		  sprintf (ContentPrefix, "%d %d", i, j);
		  char* EigenstateOutputFile = new char [512 + strlen(FilePrefix) + strlen(FileBandParameters)];
		  sprintf (EigenstateOutputFile, "%s_%s_kx_%d_ky_%d", FilePrefix, FileBandParameters, i, j);
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
	  else
	    {
	      for (int Sz = -NbrParticles; Sz <= NbrParticles; Sz += 2)
		{
		  cout << "(kx=" << i << ",ky=" << j << ") Sz=" << Sz << " : " << endl;
		  FermionOnSquareLatticeWithSpinMomentumSpace Space(NbrParticles, (Sz + NbrParticles) / 2, NbrSitesX, NbrSitesY, i, j);
		  cout << "dim = " << Space.GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
		  AbstractQHEHamiltonian* Hamiltonian = 0;
		  Hamiltonian = new ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian(&Space, NbrParticles, NbrSitesX, NbrSitesY,
													    Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
													    Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
													    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  char* ContentPrefix = new char[256];
		  sprintf (ContentPrefix, "%d %d %d", i, j, Sz);
		  char* EigenstateOutputFile = new char [512 + strlen(FilePrefix) + strlen(FileBandParameters)];
		  sprintf (EigenstateOutputFile, "%s_%s_kx_%d_ky_%d_sz_%d", FilePrefix, FileBandParameters, i, j, Sz);
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
// nnHoping = nearest neighbor hoping amplitude
// nnnHoping =  next nearest neighbor hoping amplitude
// nnnnHoping =  second next nearest neighbor hoping amplitude
// mus = sublattice staggered chemical potential 
// mixingTermNorm = norm of the mixing term coupling the two copies of the checkerboard lattice
// mixingTermArgv = argument of the mixing term coupling the two copies of the checkerboard lattice

void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, double nnHoping, double nnnHoping, double nnnnHoping, double mus, double mixingTermNorm, double mixingTermArg)
{
  ofstream File;
  File.open(outputFileName);
  File << "# kx    ky     E_{-,1}    E_{-,2}    E_{+,1}    E_{+,2}" << endl;
  double MinEMinus = 0.0;
  double MaxEMinus = -10.0;
  double MinEPlus = 10.0;
  double MaxEPlus = 0.0;
  Complex MixingTerm = mixingTermNorm * Phase(mixingTermArg);
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
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
  cout << "Spread = " << (MaxEMinus - MinEMinus) << "  Gap = " <<  (MinEPlus - MaxEMinus) << "  Flatening = " << ((MaxEMinus - MinEMinus) / (MinEPlus - MaxEMinus)) << endl;
}


// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nnHoping = nearest neighbor hoping amplitude
// nnnHoping =  next nearest neighbor hoping amplitude
// nnnnHoping =  second next nearest neighbor hoping amplitude
// mus = sublattice staggered chemical potential 
// mixingTermNorm = norm of the mixing term coupling the two copies of the checkerboard lattice
// mixingTermArgv = argument of the mixing term coupling the two copies of the checkerboard lattice

ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, double nnHoping, double nnnHoping, double nnnnHoping, double mus, double mixingTermNorm, double mixingTermArg)
{
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[nbrSitesX * nbrSitesY];
  Complex MixingTerm = mixingTermNorm * Phase(mixingTermArg);
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  int Index = (kx * nbrSitesY) + ky;
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
	  OneBodyBasis[Index] = ComplexMatrix (4, 4);
	  OneBodyBasis[Index].SetToIdentity();
#ifdef __LAPACK__
	  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, OneBodyBasis[Index]);
#else
	  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, OneBodyBasis[Index]);
#endif   
	}
    }
  return OneBodyBasis;
}
