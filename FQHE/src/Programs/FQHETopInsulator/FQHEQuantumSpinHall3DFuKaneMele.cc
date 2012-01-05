#include "Options/Options.h"

#include "HilbertSpace/FermionOnCubicLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong.h"

#include "Hamiltonian/ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian.h"
#include "Hamiltonian/ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian.h"
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
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// nnHopingDistortion111 = distortion of nearest neighbor hoping amplitude in the (111) direction
// spinOrbit = amplitude of the spin orbit coupling
void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, int nbrSitesZ, double nnHopingDistortion111, double spinOrbit);

// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// nnHopingDistortion111 = distortion of nearest neighbor hoping amplitude in the (111) direction
// spinOrbit = amplitude of the spin orbit coupling
ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, int nbrSitesZ, double nnHopingDistortion111, double spinOrbit);


int main(int argc, char** argv)
{
  cout.precision(14);

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
  (*SystemGroup) += new SingleDoubleOption  ('d', "deltat-111", "distortion of the nearest neighbor coupling in the (111) direction", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('l', "lambda-so", "spin orbit coupling (in t unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-z", "boundary condition twisting angle along z (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new BooleanOption ('\n', "four-bands", "perform the calculations within the full four band model");
  (*SystemGroup) += new BooleanOption ('\n', "project-fourbands", "project the hamiltonian from the four band model to the two band model");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "show matrix representation of the hamiltonian");
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
  int TotalNbrSites = NbrSitesX * NbrSitesY * NbrSitesZ;
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky kz ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetBoolean("four-bands") == true)
    {
      if (Manager.GetBoolean("flat-band") == true)
	{
	  if (Manager.GetDouble("mu-s") == 0.0)
	    sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_fourbands_n_%d_x_%d_y_%d_z_%d_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));
	  else
	    sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_fourbands_n_%d_x_%d_y_%d_z_%d_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f_mus_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), Manager.GetDouble("mu-s"));
	}
      else
	{
	  if (Manager.GetDouble("mu-s") == 0.0)
	    sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_fourbands_n_%d_x_%d_y_%d_z_%d_u_%f_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));
	  else
	    sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_fourbands_n_%d_x_%d_y_%d_z_%d_u_%f_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f_mus_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), Manager.GetDouble("mu-s"));
	}
    }
  else
    {
      if (Manager.GetBoolean("flat-band") == true)
	{
	  if (Manager.GetDouble("mu-s") == 0.0)
	    sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));
	  else
	    sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f_mus_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), Manager.GetDouble("mu-s"));
	}
      else
	{
	  if (Manager.GetDouble("mu-s") == 0.0)
	    sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_u_%f_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));
	  else
	    sprintf (EigenvalueOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_u_%f_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f_mus_%f.dat", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), Manager.GetDouble("mu-s"));
	}
    }
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"));
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
	      if (Manager.GetBoolean("four-bands") == false)
		{
		  FermionOnCubicLatticeWithSpinMomentumSpace Space(NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
		  cout << "dim = " << Space.GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
		  AbstractQHEHamiltonian* Hamiltonian = 0;
		  Hamiltonian = new ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian(&Space, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ,
										       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"),
										       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), 		     
										       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  char* ContentPrefix = new char[256];
		  sprintf (ContentPrefix, "%d %d %d", i, j, k);
		  char* EigenstateOutputFile = new char [512];
		  if (Manager.GetBoolean("flat-band") == true)
		    {
		      sprintf (EigenstateOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_t1_%f_t2_%f_gx_%f_gy_%f_gz_%f_kx_%d_ky_%d_kz_%d", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, 
			       Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), i, j, k);
		    }
		  else
		    {
		      sprintf (EigenstateOutputFile, "fermions_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_u_%f_t1_%f_t2_%f_gx_%f_gy_%f_gz_%f_kx_%d_ky_%d_kz_%d", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, 
			       Manager.GetDouble("u-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), i, j, k);
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
		  ParticleOnSphere* Space = 0;
#ifdef __128_BIT_LONGLONG__
		  if (TotalNbrSites <= 15)
#else
		  if (TotalNbrSites <= 7)
#endif
		    {
		      Space = new FermionOnCubicLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
		    }
		  else
		    {
		      Space = new FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
		    }
		  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
		  AbstractHamiltonian* Hamiltonian = 0;
		  Hamiltonian = new ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian((ParticleOnSphereWithSU4Spin*) Space, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ,
											Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 1.0, Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"),
											Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), 		     
											Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);

		  if (Manager.GetBoolean("project-fourbands") == true)
		    {
		      ComplexMatrix* OneBodyBasis = ComputeSingleParticleTransformationMatrices(NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"));
		      ComplexMatrix NBodyTransformationMatrix = ((FermionOnCubicLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixOneBodyBasis(OneBodyBasis);
		      ComplexMatrix HRep (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension());
		      Hamiltonian->GetHamiltonian(HRep);
		      ComplexMatrix TransformedHRep = HRep.Conjugate(NBodyTransformationMatrix);
		      FermionOnCubicLatticeWithSpinMomentumSpace* TargetSpace = new FermionOnCubicLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
		      ComplexMatrix SU4SU2TransformationMatrix = ((FermionOnCubicLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixSU4ToSU2(TargetSpace, 0, 1);
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
		  sprintf (ContentPrefix, "%d %d %d", i, j, k);
		  char* EigenstateOutputFile = new char [512];
		  if (Manager.GetBoolean("flat-band") == true)
		    {
		      sprintf (EigenstateOutputFile, "fermions_quantumspinhall3d_fukanemele_fourbands_n_%d_x_%d_y_%d_z_%d_t1_%f_t2_%f_gx_%f_gy_%f_gz_%f_kx_%d_ky_%d_kz_%d", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, 
			       Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), i, j, k);
		    }
		  else
		    {
		      sprintf (EigenstateOutputFile, "fermions_quantumspinhall3d_fukanemele_fourbands_n_%d_x_%d_y_%d_z_%d_u_%f_t1_%f_t2_%f_gx_%f_gy_%f_gz_%f_kx_%d_ky_%d_kz_%d", NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, 
			       Manager.GetDouble("u-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), i, j, k);
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
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// nnHopingDistortion111 = distortion of nearest neighbor hoping amplitude in the (111) direction
// spinOrbit = amplitude of the spin orbit coupling

void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, int nbrSitesZ, double nnHopingDistortion111, double spinOrbit)
{
  ofstream File;
  File.open(outputFileName);
  File << "# kx    ky     kz     E_{-,1}    E_{-,2}    E_{+,1}    E_{+,2}" << endl;
  double MinEMinus = 0.0;
  double MaxEMinus = -10.0;
  double MinEPlus = 10.0;
  double MaxEPlus = 0.0;
  double KxFactor = 2.0 * M_PI / ((double) nbrSitesX);
  double KyFactor = 2.0 * M_PI / ((double) nbrSitesY);
  double KzFactor = 2.0 * M_PI / ((double) nbrSitesZ);
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  for (int kz = 0; kz < nbrSitesZ; ++kz)
	    {
	      HermitianMatrix TmpOneBodyHamiltonian(4, true);
	      Complex B1 = 1.0 + nnHopingDistortion111 + Phase(0.5 * (((double) ky) * KyFactor) + 0.5 * (((double) kz) * KzFactor))   + Phase(0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) kz) * KzFactor))  + Phase(0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) ky) * KyFactor)) ;
	      double d3 = spinOrbit * (sin (0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) kz) * KzFactor))
				       - sin (0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) ky) * KyFactor))
				       - sin (0.5 * (((double) kx) * KxFactor) - 0.5 * (((double) ky) * KyFactor))
				       + sin (0.5 * (((double) kx) * KxFactor) - 0.5 * (((double) kz) * KzFactor)));
	      double d4 = spinOrbit * (sin (0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) kz) * KzFactor))
				       - sin (0.5 * (((double) ky) * KyFactor) + 0.5 * (((double) kz) * KzFactor))
				       - sin (0.5 * (((double) ky) * KyFactor) - 0.5 * (((double) kz) * KzFactor))
				       + sin (0.5 * (((double) ky) * KyFactor) - 0.5 * (((double) kx) * KxFactor)));
	      double d5 = spinOrbit * (sin (0.5 * (((double) ky) * KyFactor) + 0.5 * (((double) kx) * KxFactor))
				       - sin (0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) kz) * KzFactor))
				       - sin (0.5 * (((double) kz) * KzFactor) - 0.5 * (((double) kx) * KxFactor))
				       + sin (0.5 * (((double) kz) * KzFactor) - 0.5 * (((double) ky) * KyFactor)));
	      Complex B2 = d3 + I() * d4;
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 2, -d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 3, d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 2, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 3, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B2);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 3, -B2);
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
	      File << (KxFactor * ((double) kx)) << " " << (KyFactor * ((double) ky)) << " " << (KzFactor * ((double) kz)) << " " << TmpDiag(0, 0) << " " << TmpDiag(1, 1) <<  " " << TmpDiag(2, 2) << " " << TmpDiag(3, 3) << endl;
	    }
	  File << endl;
	}
    }
  cout << "Spread = " << (MaxEMinus - MinEMinus) << "  Gap = " <<  (MinEPlus - MaxEMinus) << "  Flatening = " << ((MaxEMinus - MinEMinus) / (MinEPlus - MaxEMinus)) << endl;
}

// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// nnHopingDistortion111 = distortion of nearest neighbor hoping amplitude in the (111) direction
// spinOrbit = amplitude of the spin orbit coupling

ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, int nbrSitesZ, double nnHopingDistortion111, double spinOrbit)
{
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[nbrSitesX * nbrSitesY * nbrSitesZ];
  double KxFactor = 2.0 * M_PI / ((double) nbrSitesX);
  double KyFactor = 2.0 * M_PI / ((double) nbrSitesY);
  double KzFactor = 2.0 * M_PI / ((double) nbrSitesZ);
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  for (int kz = 0; kz < nbrSitesZ; ++kz)
	    {
	      int Index = ((kx * nbrSitesY) + ky) * nbrSitesZ + kz;
	      HermitianMatrix TmpOneBodyHamiltonian(4, true);
	      Complex B1 = 1.0 + nnHopingDistortion111 + Phase(0.5 * (((double) ky) * KyFactor) + 0.5 * (((double) kz) * KzFactor))   + Phase(0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) kz) * KzFactor))  + Phase(0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) ky) * KyFactor)) ;
	      double d3 = spinOrbit * (sin (0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) kz) * KzFactor))
				       - sin (0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) ky) * KyFactor))
				       - sin (0.5 * (((double) kx) * KxFactor) - 0.5 * (((double) ky) * KyFactor))
				       + sin (0.5 * (((double) kx) * KxFactor) - 0.5 * (((double) kz) * KzFactor)));
	      double d4 = spinOrbit * (sin (0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) kz) * KzFactor))
				       - sin (0.5 * (((double) ky) * KyFactor) + 0.5 * (((double) kz) * KzFactor))
				       - sin (0.5 * (((double) ky) * KyFactor) - 0.5 * (((double) kz) * KzFactor))
				       + sin (0.5 * (((double) ky) * KyFactor) - 0.5 * (((double) kx) * KxFactor)));
	      double d5 = spinOrbit * (sin (0.5 * (((double) ky) * KyFactor) + 0.5 * (((double) kx) * KxFactor))
				       - sin (0.5 * (((double) kx) * KxFactor) + 0.5 * (((double) kz) * KzFactor))
				       - sin (0.5 * (((double) kz) * KzFactor) - 0.5 * (((double) kx) * KxFactor))
				       + sin (0.5 * (((double) kz) * KzFactor) - 0.5 * (((double) ky) * KyFactor)));
	      Complex B2 = d3 + I() * d4;
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 2, -d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 3, d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 2, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 3, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B2);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 3, -B2);
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
    }
  return OneBodyBasis;
}

