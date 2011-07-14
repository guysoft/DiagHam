#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

//#include "Hamiltonian/ParticleOnLatticeWithSpinHaldaneModelHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeHaldaneModelSingleBandHamiltonian.h"
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
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the x direction
// nnHoping = nearest neighbor hoping amplitude
// nnnHoping =  next nearest neighbor hoping amplitude
// phi =  Haldane phase on nnn hopping
// mus = sublattice staggered chemical potential 
void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSiteX, int nbrSiteY, double nnHoping, double nnnHoping, double phi, double mus);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEHaldaneModel" , "0.01");
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
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "phi", "Haldane phase on nnn hopping", M_PI/3);
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
      cout << "see man page for option syntax or type FQHEHaldaneModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSiteX = Manager.GetInteger("nbr-sitex"); 
  int NbrSiteY = Manager.GetInteger("nbr-sitey"); 
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n#");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
  {
    if (Manager.GetBoolean("single-band") == false)
      {
        if (Manager.GetDouble("mu-s") == 0.0)
          sprintf (EigenvalueOutputFile, "fermions_haldane_n_%d_x_%d_y_%d_u_%f_v_%f_t1_%f_t2_%f_phi_%f.dat", NbrParticles, NbrSiteX, NbrSiteY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("phi"));
        else
          sprintf (EigenvalueOutputFile, "fermions_haldane_n_%d_x_%d_y_%d_u_%f_v_%f_t1_%f_t2_%f_phi_%f_mus_%f.dat", NbrParticles, NbrSiteX, NbrSiteY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("phi"), Manager.GetDouble("mu-s"));
      }
    else
      {
        if (Manager.GetBoolean("flat-band") == true)
          {
            if (Manager.GetDouble("mu-s") == 0.0)
              sprintf (EigenvalueOutputFile, "fermions_singleband_haldane_n_%d_x_%d_y_%d_v_%f_t1_%f_t2_%f_phi_%f_gx_%f_gy_%f.dat", NbrParticles, NbrSiteX, NbrSiteY, Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("phi"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
            else
              sprintf (EigenvalueOutputFile, "fermions_singleband_haldane_n_%d_x_%d_y_%d_v_%f_t1_%f_t2_%f_phi_%f_gx_%f_gy_%f_mus_%f.dat", NbrParticles, NbrSiteX, NbrSiteY, Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("phi"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
          }
        else
          {
            if (Manager.GetDouble("mu-s") == 0.0)
              sprintf (EigenvalueOutputFile, "fermions_singleband_haldane_n_%d_x_%d_y_%d_u_%f_v_%f_t1_%f_t2_%f_phi_%f_gx_%f_gy_%f.dat", NbrParticles, NbrSiteX, NbrSiteY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("phi"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
            else
              sprintf (EigenvalueOutputFile, "fermions_singleband_haldane_n_%d_x_%d_y_%d_u_%f_v_%f_t1_%f_t2_%f_phi_%f_gx_%f_gy_%f_mus_%f.dat", NbrParticles, NbrSiteX, NbrSiteY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("phi"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
          }
      }
  }

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSiteX, NbrSiteY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("phi"), Manager.GetDouble("mu-s"));
      return 0;
    }

  int MinKx = 0;
  int MaxKx = NbrSiteX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrSiteY - 1;
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
                cout << "NotImplemented: only the single-band case is supported at the moment." <<endl;
                return 1;
//	      FermionOnSquareLatticeWithSpinMomentumSpace Space(NbrParticles, NbrSiteX, NbrSiteY, i, j);
//	      cout << "dim = " << Space.GetHilbertSpaceDimension()  << endl;
//	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
//		Memory = Architecture.GetArchitecture()->GetLocalMemory();
//	      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
//	      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnLatticeWithSpinHaldaneModelHamiltonian(&Space, NbrParticles, NbrSiteX, NbrSiteY,
//														Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
//														Manager.GetDouble("phi"), Manager.GetDouble("gamma-x") * 2.0 * M_PI, Manager.GetDouble("gamma-y") * 2.0 * M_PI, 		     
//														Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
//	      char* ContentPrefix = new char[256];
//	      sprintf (ContentPrefix, "%d %d", i, j);
//	      char* EigenstateOutputFile = new char [512];
//	      sprintf (EigenstateOutputFile, "fermions_haldane_n_%d_x_%d_y_%d_u_%f_v_%f_t1_%f_t2_%f_kx_%d_ky_%d", NbrParticles, NbrSiteX, NbrSiteY, 
//		       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), i, j);
//	      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
//	      FirstRunFlag = false;
//	      MainTaskOperation TaskOperation (&Task);
//	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
//	      cout << "------------------------------------" << endl;
//	      delete Hamiltonian;
//	      delete[] EigenstateOutputFile;
//	      delete[] ContentPrefix;
 	    }
 	  else
 	    {
 	      FermionOnSquareLatticeMomentumSpace Space(NbrParticles, NbrSiteX, NbrSiteY, i, j);
 	      cout << "dim = " << Space.GetHilbertSpaceDimension()  << endl;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
 	      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
 	      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnLatticeHaldaneModelSingleBandHamiltonian(&Space, NbrParticles, NbrSiteX, NbrSiteY,
														  Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"),
														  Manager.GetDouble("phi"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
														  Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", i, j);
	      char* EigenstateOutputFile = new char [512];
	      if (Manager.GetBoolean("flat-band") == true)
		{
		  if (Manager.GetDouble("mu-s") == 0.0)
		    sprintf (EigenstateOutputFile, "fermions_singleband_haldane_n_%d_x_%d_y_%d_v_%f_t1_%f_t2_%f_gx_%f_gy_%f_kx_%d_ky_%d", NbrParticles, NbrSiteX, NbrSiteY, 
			     Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
		  else
		    sprintf (EigenstateOutputFile, "fermions_singleband_haldane_n_%d_x_%d_y_%d_v_%f_t1_%f_t2_%f_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d", NbrParticles, NbrSiteX, NbrSiteY, 
			     Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"), i, j);
		}
	      else
		{
		  if (Manager.GetDouble("mu-s") == 0.0)
		    sprintf (EigenstateOutputFile, "fermions_singleband_haldane_n_%d_x_%d_y_%d_u_%f_v_%f_t1_%f_t2_%f_gx_%f_gy_%f_kx_%d_ky_%d", NbrParticles, NbrSiteX, NbrSiteY, 
			     Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
		  else
		    sprintf (EigenstateOutputFile, "fermions_singleband_haldane_n_%d_x_%d_y_%d_u_%f_v_%f_t1_%f_t2_%f_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d", NbrParticles, NbrSiteX, NbrSiteY, 
			     Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"), i, j);
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
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the x direction
// nnHoping = nearest neighbor hoping amplitude
// nnnHoping =  next nearest neighbor hoping amplitude
// phase =  Haldane phase on nnn hopping
// mus = sublattice staggered chemical potential 

void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSiteX, int nbrSiteY, double nnHoping, double nnnHoping, double phase, double mus)
{
  ofstream File;
  File.open(outputFileName);
  File << "# kx    ky     E_-    E_-" << endl;
  double MinEMinus = 0.0;
  double MaxEMinus = -10.0;
  double MinEPlus = 10.0;
  double MaxEPlus = 0.0;
  for (int kx = 0; kx < nbrSiteX; ++kx)
  {
      double x=2*M_PI*((double)kx)/nbrSiteX;
      for (int ky = 0; ky < nbrSiteY; ++ky)
      {
          double y=2*M_PI*((double)ky)/nbrSiteY;
          Complex B1 = nnHoping * Complex(1 + cos(x+y) + cos(y), + sin(x+y) + sin(y));
          double d0 = + 2.0 * nnnHoping * cos(phase) * (cos(x) + cos(y) + cos(x+y));
          double d3 = + 2.0 * nnnHoping * sin(phase) * (sin(x) + sin(y) - sin(x+y)) + mus;
	  HermitianMatrix TmpOneBobyHamiltonian(2, true);
	  TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d0 + d3);
	  TmpOneBobyHamiltonian.SetMatrixElement(0, 1, B1);
	  TmpOneBobyHamiltonian.SetMatrixElement(1, 1, d0 - d3);
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
	  File << (2.0 * M_PI * ((double) kx) / ((double) nbrSiteX)) << " " << (2.0 * M_PI * ((double) ky) / ((double) nbrSiteY)) << " " << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << endl;
	}
      File << endl;
    }
  cout << "Spread = " << (MaxEMinus - MinEMinus) << "  Gap = " <<  (MinEPlus - MaxEMinus) << "  Flatening = " << ((MaxEMinus - MinEMinus) / (MinEPlus - MaxEMinus)) << endl;
}
