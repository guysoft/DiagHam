#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian.h"

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
void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, double nnHoping, double nnnHoping, double nnnnHoping, double mus);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHECheckerboardLatticeModel" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three body interaction instead of a two body interaction");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECheckerboardLatticeModel -h" << endl;
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
      sprintf (FilePrefix, "fermions_checkerboardlattice_n_%d_x_%d_y",  NbrParticles, NbrSitesX, NbrSitesY);
    }
  else
    {
      if (Manager.GetBoolean("three-body") == false)
	sprintf (FilePrefix, "fermions_singleband_checkerboardlattice_n_%d_x_%d_y_%d",  NbrParticles, NbrSitesX, NbrSitesY);
      else
	sprintf (FilePrefix, "fermions_singleband_threebody_checkerboardlattice_n_%d_x_%d_y_%d",  NbrParticles, NbrSitesX, NbrSitesY);
    }
  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
  {
    if (Manager.GetBoolean("single-band") == false)
      {
        if (Manager.GetDouble("mu-s") == 0.0)
          sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_t1_%f_t2_%f_tpp_%f.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"));
        else
          sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_t1_%f_t2_%f_tpp_%f_mus_%f.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"));
      }
    else
      {
        if (Manager.GetBoolean("flat-band") == true)
          {
            if (Manager.GetDouble("mu-s") == 0.0)
              sprintf (EigenvalueOutputFile, "%s_v_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f.dat",FilePrefix, Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
            else
              sprintf (EigenvalueOutputFile, "%s_v_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f_mus_%f.dat",FilePrefix, Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
          }
        else
          {
            if (Manager.GetDouble("mu-s") == 0.0)
              sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
            else
              sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f_mus_%f.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
          }
      }
  }

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"));
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
	  cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
 	  if (Manager.GetBoolean("single-band") == false)
 	    {
	      ParticleOnSphereWithSpin* Space = 0;
	      if ((NbrSitesX * NbrSitesY) <= 31)
		{
		  Space = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      else
		{
		  Space = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
														Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
														Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x") * 2.0 * M_PI, Manager.GetDouble("gamma-y") * 2.0 * M_PI, 		     
														Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", i, j);
	      char* EigenstateOutputFile = new char [512];
	      sprintf (EigenstateOutputFile, "%s_u_%f_v_%f_t1_%f_t2_%f_kx_%d_ky_%d",FilePrefix, 
		       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), i, j);
	      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	      FirstRunFlag = false;
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      cout << "------------------------------------" << endl;
	      delete Hamiltonian;
	      delete Space;
	      delete[] EigenstateOutputFile;
	      delete[] ContentPrefix;
 	    }
 	  else
 	    {
	      ParticleOnSphere* Space = 0;
	      if ((NbrSitesX * NbrSitesY) <= 63)
		{
		  Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      else
		{
		  Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
 	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
// 	      for (long k = 0; k < Space->GetHilbertSpaceDimension(); ++k)
// 		{
// 		  Space->PrintState(cout, k) << endl;
// 		}
// 	      return 0;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
 	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
 	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      if (Manager.GetBoolean("three-body") == false)
		{ 
		  Hamiltonian = new ParticleOnLatticeCheckerboardLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
											      Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
											      Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
											      Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		}
	      else
		{ 
		  Hamiltonian = new ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
												       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
												       Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
												       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		}
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", i, j);
	      char* EigenstateOutputFile = new char [512];
	      if (Manager.GetBoolean("flat-band") == true)
		{
		  if (Manager.GetDouble("mu-s") == 0.0)
		    sprintf (EigenstateOutputFile, "%s_v_%f_t1_%f_t2_%f_gx_%f_gy_%f_kx_%d_ky_%d",FilePrefix, 
			     Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
		  else
		    sprintf (EigenstateOutputFile, "%s_v_%f_t1_%f_t2_%f_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d",FilePrefix, 
			     Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"), i, j);
		}
	      else
		{
		  if (Manager.GetDouble("mu-s") == 0.0)
		    sprintf (EigenstateOutputFile, "%s_u_%f_v_%f_t1_%f_t2_%f_gx_%f_gy_%f_kx_%d_ky_%d",FilePrefix, 
			     Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
		  else
		    sprintf (EigenstateOutputFile, "%s_u_%f_v_%f_t1_%f_t2_%f_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d",FilePrefix, 
			     Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"), i, j);
		}
	      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	      FirstRunFlag = false;
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      cout << "------------------------------------" << endl;
	      delete Hamiltonian;
	      delete Space;
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
// nbrSitesY = number of sites in the x direction
// nnHoping = nearest neighbor hoping amplitude
// nnnHoping =  next nearest neighbor hoping amplitude
// nnnnHoping =  second next nearest neighbor hoping amplitude
// mus = sublattice staggered chemical potential 

void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, double nnHoping, double nnnHoping, double nnnnHoping, double mus)
{
  ofstream File;
  File.open(outputFileName);
  File << "# kx    ky     E_-    E_-" << endl;
  double MinEMinus = 0.0;
  double MaxEMinus = -10.0;
  double MinEPlus = 10.0;
  double MaxEPlus = 0.0;
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  Complex B1 = 4.0 * nnHoping * Complex (cos (1.0 * M_PI * ((double) kx) / ((double) nbrSitesX)) * cos (1.0 * M_PI * ((double) ky) / ((double) nbrSitesY)) * cos(M_PI * 0.25), 
					   sin (1.0 * M_PI * ((double) kx) / ((double) nbrSitesX)) * sin (1.0 * M_PI * ((double) ky) / ((double) nbrSitesY)) * sin(M_PI * 0.25));
	  double d1 = 4.0 * nnnnHoping * cos (2.0 * M_PI * ((double) kx) / ((double) nbrSitesX)) * cos (2.0 * M_PI * ((double) ky) / ((double) nbrSitesY));
	  double d3 = mus + (2.0 * nnnHoping * (cos (2.0 * M_PI * ((double) kx) / ((double) nbrSitesX))
						- cos (2.0 * M_PI * ((double) ky) / ((double) nbrSitesY))));
	  HermitianMatrix TmpOneBobyHamiltonian(2, true);
	  TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	  TmpOneBobyHamiltonian.SetMatrixElement(0, 1, B1);
	  TmpOneBobyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag);
#else
	  TmpOneBobyHamiltonian.Diagonalize(TmpDiag);
#endif   
	  if (MaxEMinus < TmpDiag(0, 0))
	    {
	      MaxEMinus = TmpDiag(0, 0);
	    }
	  if (MinEMinus > TmpDiag(0, 0))
	    {
	      MinEMinus = TmpDiag(0, 0);
	    }
	  if (MaxEPlus < TmpDiag(1, 1))
	    {
	      MaxEPlus = TmpDiag(1, 1);
	    }
	  if (MinEPlus > TmpDiag(1, 1))
	    {
	      MinEPlus = TmpDiag(1, 1);
	    }
	  File << (2.0 * M_PI * ((double) kx) / ((double) nbrSitesX)) << " " << (2.0 * M_PI * ((double) ky) / ((double) nbrSitesY)) << " " << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << endl;
	}
      File << endl;
    }
  cout << "Spread = " << (MaxEMinus - MinEMinus) << "  Gap = " <<  (MinEPlus - MaxEMinus) << "  Flatening = " << ((MaxEMinus - MinEMinus) / (MinEPlus - MaxEMinus)) << endl;
}
