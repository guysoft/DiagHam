#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian.h"
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
// l1 = Rashba coupling between nearest neighbor sites
// l2 = Rashba coupling between next nearest neighbor sites
void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSiteX, int nbrSiteY, double nnHoping, double nnnHoping, double l1, double l2);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEAlternativeKagomeLatticeModel" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive two-body nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive two-body nearest next neighbor potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive three-body nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "s-potential", "repulsive three-body next-to-nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three-body interaction in addition to a two-body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four-body interaction in addition to a two-body interaction");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", -0.3);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l1", "Rashba coupling between nearest neighbor sites", 0.28);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l2", "Rashba coupling between next nearest neighbor sites", 0.2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-index", "index of the band that has to be partially filled, should be 0 (lower band), 1 or 2 (upper band)", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model. The n-body interaction strength with largest n is set to unity");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
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
      cout << "see man page for option syntax or type FQHEAlternativeKagomeLatticeModel -h" << endl;
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
  int BandIndex = Manager.GetInteger("band-index");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }
  
  char* FilePrefix = new char [512];
  int lenFilePrefix=0;
  if (Manager.GetBoolean("single-band") == false)
    {
      cout << "NotImplemented: only the single-band case is supported at the moment." <<endl;
      return 1;
    }
  else
    {
      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
	lenFilePrefix += sprintf (FilePrefix, "%s_singleband_kagome_band_%d_n_%d_x_%d_y_%d", StatisticPrefix, BandIndex, NbrParticles, NbrSiteX, NbrSiteY);
      else
	{
	  if (Manager.GetBoolean("three-body") == true)
	    lenFilePrefix += sprintf (FilePrefix, "%s_singleband_threebody_kagome_band_%d_n_%d_x_%d_y_%d", StatisticPrefix, BandIndex, NbrParticles, NbrSiteX, NbrSiteY);
          else
	    lenFilePrefix += sprintf (FilePrefix, "%s_singleband_fourbody_kagome_band_%d_n_%d_x_%d_y_%d", StatisticPrefix, BandIndex, NbrParticles, NbrSiteX, NbrSiteY);
	}
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) && Manager.GetBoolean("flat-band") == false)
	lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_w_%f", Manager.GetDouble("w-potential"));
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) && Manager.GetDouble("s-potential") != 0.0)
	lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_s_%f", Manager.GetDouble("s-potential"));
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) || Manager.GetBoolean("flat-band") == false)
	lenFilePrefix += sprintf (FilePrefix + lenFilePrefix, "_u_%f", Manager.GetDouble("u-potential"));
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_v_%f", Manager.GetDouble("v-potential"));
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
    }
  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
      sprintf (EigenvalueOutputFile, "%s.dat", FilePrefix);

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSiteX, NbrSiteY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"));
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
 	    }
 	  else
 	    {
	      ParticleOnSphere* Space = 0;
	      if (Manager.GetBoolean("boson") == false)
		{
		  if ((NbrSiteX * NbrSiteY) <= 63)
		    {
		      Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, i, j);
		    }
		  else
		    {
		      Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, i, j);
		    }
		}
	      else
		{
		  Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, i, j);
		}
 	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
 	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
 	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
              {
                  Hamiltonian = new ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
                          Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
												   Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), BandIndex,
                          Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
              }
              else
              {
		  if (Manager.GetBoolean("three-body") == true)
                  {
                      Hamiltonian = new ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
                              Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
														Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), BandIndex,
                              Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
                  }
                  else
                  {
                      cout << "NotImplemented: four-body interaction is not supported at the moment." <<endl;
                      return 1;
                  }
              }

	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", i, j);
	      char* EigenstateOutputFile = new char [512];
              if (Manager.GetString("eigenstate-file")!=0)
                  sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), i, j);
              else
                  sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", FilePrefix, i, j);
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
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the x direction
// nnHoping = nearest neighbor hoping amplitude
// nnnHoping =  next nearest neighbor hoping amplitude
// nnRashba = Rashba coupling between nearest neighbor sites
// nnnRashba = Rashba coupling between next nearest neighbor sites

void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSiteX, int nbrSiteY, double nnHoping, double nnnHoping, double nnRashba, double nnnRashba)
{
  ofstream File;
  File.open(outputFileName);
  File << "# kx    ky     E_1    E_2    E_3" << endl;
  double MinE0 = 100.0;
  double MinE1 = 100.0;
  double MinE2 = 100.0;
  double MaxE0 = -100.0;
  double MaxE1 = -100.0;
  double MaxE2 = -100.0;
  for (int kx = 0; kx < nbrSiteX; ++kx)
  {
      double x=2*M_PI*((double)kx)/nbrSiteX;
      for (int ky = 0; ky < nbrSiteY; ++ky)
      {
          double y=2*M_PI*((double)ky)/nbrSiteY;
          int Index = (kx * nbrSiteY) + ky;
  
          Complex nnBA = Complex(-nnHoping, -nnRashba) * (1 + Phase(x));
          Complex nnCA = Complex(-nnHoping, +nnRashba) * (1 + Phase(y));
          Complex nnCB = Complex(-nnHoping, -nnRashba) * (1 + Phase(y-x));
          Complex nnnBA = Complex(-nnnHoping, +nnnRashba) * (Phase(y) + Phase(x-y));
          Complex nnnCA = Complex(-nnnHoping, -nnnRashba) * (Phase(x) + Phase(y-x));
          Complex nnnCB = Complex(-nnnHoping, +nnnRashba) * (Phase(-x) + Phase(y));
  
          HermitianMatrix TmpOneBobyHamiltonian(3, true);
          TmpOneBobyHamiltonian.SetMatrixElement(1, 0, nnBA + nnnBA);
          TmpOneBobyHamiltonian.SetMatrixElement(2, 0, nnCA + nnnCA);
          TmpOneBobyHamiltonian.SetMatrixElement(2, 1, nnCB + nnnCB);
          ComplexMatrix TmpMatrix(3, 3, true);
          TmpMatrix[0][0] = 1.0;
          TmpMatrix[1][1] = 1.0;
          TmpMatrix[2][2] = 1.0;
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag);
#else
	  TmpOneBobyHamiltonian.Diagonalize(TmpDiag);
#endif   
	  if (MaxE0 < TmpDiag(0, 0))
              MaxE0 = TmpDiag(0, 0);
	  if (MinE0 > TmpDiag(0, 0))
	      MinE0 = TmpDiag(0, 0);
	  if (MaxE1 < TmpDiag(1, 1))
	      MaxE1 = TmpDiag(1, 1);
	  if (MinE1 > TmpDiag(1, 1))
	      MinE1 = TmpDiag(1, 1);
	  if (MaxE2 < TmpDiag(2, 2))
	      MaxE2 = TmpDiag(2, 2);
	  if (MinE2 > TmpDiag(2, 2))
	      MinE2 = TmpDiag(2, 2);
	  File << (2.0 * M_PI * ((double) kx) / ((double) nbrSiteX)) << " " << (2.0 * M_PI * ((double) ky) / ((double) nbrSiteY)) << " " << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << endl;
	}
      File << endl;
    }
  cout << "Delta12 = " << (MinE1 - MaxE0) << "  Delta13 = " << (MinE2-MaxE0) << "  W = " << (MaxE0 - MinE0) << endl;
  cout << "Delta12 / W = " << ((MinE1 - MaxE0) / (MaxE0 - MinE0)) << "  Delta13 / W = " << ((MinE2 - MaxE0) / (MaxE0 - MinE0)) <<endl;
}
