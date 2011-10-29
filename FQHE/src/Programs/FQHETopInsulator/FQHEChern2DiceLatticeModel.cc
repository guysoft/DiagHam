#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeChern2DiceLatticeSingleBandFourBodyHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeChern2DiceLatticeSingleBandFiveBodyHamiltonian.h"

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
// t = nearest neighbor hopping amplitude
// epsilon = on site energy for site 3
// lambda = Rashba spin orbit coupling strength
// bfield1 = magnetic field strength on sites 1 and 2
// bfield3 = magnetic field strength on site 3
// gammaX = boundary condition twisting angle along e_a (measured in units of 2pi)
// gammaY = boundary condition twisting angle along e_b (measured in units of 2pi)
void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, double t, double epsilon, double lambda, double bfield1, double bfield3, double gammaX, double gammaY);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEChern2DiceLatticeModel" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "five-body", "use a five body interaction instead of a two body interaction");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t", "nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "epsilon", "on site energy for site 3", 0.6);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "lambda", "Rashba spin orbit coupling strength", 0.3);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "B1", "magnetic field strength on sites 1 and 2", 0.2440);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "B3", "magnetic field strength on site 3", -0.0162);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
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
      cout << "see man page for option syntax or type FQHEChern2DiceLatticeModel -h" << endl;
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

  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }

  char* FilePrefix = new char [256];

  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
    { 
      sprintf (FilePrefix, "%s_singleband_chern2dicelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
    }
  else
    {
      if (Manager.GetBoolean("three-body") == true)
	{
	  sprintf (FilePrefix, "%s_singleband_threebody_chern2dicelattice_n_%d_x_%d_y_%d",StatisticPrefix,  NbrParticles, NbrSitesX, NbrSitesY);
	}
      else
	{
	  if (Manager.GetBoolean("four-body") == true)
	    {
	      sprintf (FilePrefix, "%s_singleband_fourbody_chern2dicelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
	    }
	  else
	    {
	      sprintf (FilePrefix, "%s_singleband_fivebody_chern2dicelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
	    }
	}
    }

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    {
      if (Manager.GetBoolean("flat-band") == true)
	{
	  sprintf (EigenvalueOutputFile, "%s_t_%g_epsilon_%g_lambda_%g_B1_%g_B3_%g_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("t"), Manager.GetDouble("epsilon"), Manager.GetDouble("lambda"), Manager.GetDouble("B1"), Manager.GetDouble("B3"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	}
      else
	{
	  sprintf (EigenvalueOutputFile, "%s_u_%g_t_%g_epsilon_%g_lambda_%g_B1_%g_B3_%g_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("t"), Manager.GetDouble("epsilon"), Manager.GetDouble("lambda"), Manager.GetDouble("B1"), Manager.GetDouble("B3"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	}
    }

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSitesX, NbrSitesY, Manager.GetDouble("t"), Manager.GetDouble("epsilon"), Manager.GetDouble("lambda"), Manager.GetDouble("B1"), Manager.GetDouble("B3"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
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

	  ParticleOnSphere* Space = 0;
	  if (Manager.GetBoolean("boson") == false)
	    {
	      if ((NbrSitesX * NbrSitesY) <= 63)
		{
		  Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      else
		{
		  Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	    }
	  else
	    {
	      Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
	    }
	  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	  AbstractQHEHamiltonian* Hamiltonian = 0;
	  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
	    { 
	      Hamiltonian = new ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
										    Manager.GetDouble("u-potential"), Manager.GetDouble("t"), Manager.GetDouble("epsilon"), Manager.GetDouble("lambda"), Manager.GetDouble("B1"),
										    Manager.GetDouble("B3"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
										    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	    }
	  else
	    { 
	      if (Manager.GetBoolean("three-body") == true)
		{
		  Hamiltonian = 0;
// 		  Hamiltonian = new ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
// 												 Manager.GetDouble("u-potential"), Manager.GetDouble("t"), Manager.GetDouble("epsilon"), Manager.GetDouble("lambda"), Manager.GetDouble("B1"),
// 												 Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
// 												 Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		}
	      else
		{
		  if (Manager.GetBoolean("four-body") == true)
		    {
		      Hamiltonian = 0;
// 		      Hamiltonian = new ParticleOnLatticeChern2DiceLatticeSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
// 												    Manager.GetDouble("u-potential"), Manager.GetDouble("t"), Manager.GetDouble("epsilon"), Manager.GetDouble("lambda"), Manager.GetDouble("B1"),
// 												    Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
// 												    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		    }
		  else
		    {
		      Hamiltonian = 0;
// 		      Hamiltonian = new ParticleOnLatticeChern2DiceLatticeSingleBandFiveBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
// 												    Manager.GetDouble("u-potential"), Manager.GetDouble("t"), Manager.GetDouble("epsilon"), Manager.GetDouble("lambda"), Manager.GetDouble("B1"),
// 												    Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
// 												    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		    }
		  
		}
	    }

	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d %d", i, j);
	  char* EigenstateOutputFile = new char [512];
	  if (Manager.GetString("eigenstate-file")!=0)
	    sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), i, j);
	  else
	    {
	      if (Manager.GetBoolean("flat-band") == true)
		{
		  sprintf (EigenstateOutputFile, "%s_t_%g_epsilon_%g_lambda_%g_B1_%g_B3_%g_gx_%g_gy_%g_kx_%d_ky_%d",FilePrefix, 
			   Manager.GetDouble("t"), Manager.GetDouble("epsilon"), Manager.GetDouble("lambda"), Manager.GetDouble("B1"), Manager.GetDouble("B3"),
			   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
		}
	      else
		{
		  sprintf (EigenstateOutputFile, "%s_u_%g_t_%g_epsilon_%g_lambda_%g_B1_%g_B3_%g_gx_%g_gy_%g_kx_%d_ky_%d",FilePrefix, 
			   Manager.GetDouble("u-potential"), Manager.GetDouble("t"), Manager.GetDouble("epsilon"), Manager.GetDouble("lambda"),
			   Manager.GetDouble("B1"), Manager.GetDouble("B3"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
		}
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
  return 0;
}

// compute the single particle spectrum 
//
// outputFileName = name of the output file
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the x direction
// t = nearest neighbor hopping amplitude
// epsilon = on site energy for site 3
// lambda = Rashba spin orbit coupling strength
// bfield1 = magnetic field strength on sites 1 and 2
// bfield3 = magnetic field strength on site 3
// gammaX = boundary condition twisting angle along e_a (measured in units of 2pi)
// gammaY = boundary condition twisting angle along e_b (measured in units of 2pi)

void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, double t, double epsilon, double lambda, double bfield1, double bfield3, double gammaX, double gammaY)
{
  ofstream File;
  File.open(outputFileName);
  File << "# kx    ky     E_1    E_2    E_3    E_4    E_5    E_6" << endl;
  double MinEMinus = 0.0;
  double MaxEMinus = -10.0;
  double MinEPlus = 10.0;
  double MaxEPlus = 0.0;
  double DirectGap = 20.0;
  double KX, KY;
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  KX = 2.0 * M_PI / ((double) nbrSitesX) * ((double) kx);
	  KY = 2.0 * M_PI / ((double) nbrSitesY) * ((double) ky);
 	  Complex GammaK = t * (1.0 + Phase(KX) + Phase(KY));
 	  Complex GammaKPlus = I() * lambda * (1.0 + Phase(KX + (2.0 * M_PI / 3.0)) + Phase(KY + (4.0 * M_PI / 3.0)));
 	  Complex GammaKMinus = I() * lambda * (1.0 + Phase(KX - (2.0 * M_PI / 3.0)) + Phase(KY - (4.0 * M_PI / 3.0)));
	  
	  HermitianMatrix TmpOneBodyHamiltonian(6, true);
 	  TmpOneBodyHamiltonian.SetMatrixElement(0, 4, Conj(GammaK));
 	  TmpOneBodyHamiltonian.SetMatrixElement(1, 5, Conj(GammaK));
 	  TmpOneBodyHamiltonian.SetMatrixElement(0, 5, Conj(GammaKPlus));
 	  TmpOneBodyHamiltonian.SetMatrixElement(1, 4, Conj(GammaKMinus));

 	  TmpOneBodyHamiltonian.SetMatrixElement(2, 4, GammaK);
 	  TmpOneBodyHamiltonian.SetMatrixElement(3, 5, GammaK);
 	  TmpOneBodyHamiltonian.SetMatrixElement(2, 5, GammaKMinus);
 	  TmpOneBodyHamiltonian.SetMatrixElement(3, 4, GammaKPlus);

	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, bfield1);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -bfield1);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 2, bfield1);
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 3, -bfield1);
	  TmpOneBodyHamiltonian.SetMatrixElement(4, 4, epsilon + bfield3);
	  TmpOneBodyHamiltonian.SetMatrixElement(5, 5, epsilon - bfield3);

	  TmpOneBodyHamiltonian *= -1.0;

	  ComplexMatrix TmpMatrix(6, 6, true);
	  TmpMatrix.SetToIdentity();
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
	  if (MaxEMinus < TmpDiag(2, 2))
	    {
	      MaxEMinus = TmpDiag(2, 2);
	    }
	  if (MinEMinus > TmpDiag(2, 2))
	    {
	      MinEMinus = TmpDiag(2, 2);
	    }
	  if (MaxEPlus < TmpDiag(3, 3))
	    {
	      MaxEPlus = TmpDiag(3, 3);
	    }
	  if (MinEPlus > TmpDiag(3, 3))
	    {
	      MinEPlus = TmpDiag(3, 3);
	    }
	  if ((TmpDiag(3, 3) - TmpDiag(2, 2)) < DirectGap)
	    DirectGap = TmpDiag(3, 3) - TmpDiag(2, 2);
	  double Kx = KX;
	  if (Kx > M_PI) 
	    Kx -= M_PI;
	  double Ky = (KX + 2.0 * KY) /sqrt(3.0);
	  if (Ky  > M_PI) 
	    Ky -= M_PI;
	  if (Ky < -M_PI) 
	    Ky += M_PI;
	  File << Kx << " " << Ky << " " << KX << " " << KY << " " << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2)  
	       << " " << TmpDiag(3, 3) << " " << TmpDiag(4, 4) << " " << TmpDiag(5, 5) << endl;
	}
      File << endl;
    }
  cout << "Spread = " << (MaxEMinus - MinEMinus) << "  Gap = " <<  (MinEPlus - MaxEMinus) << "  Flattening = " << ((MaxEMinus - MinEMinus) / (MinEPlus - MaxEMinus)) << "  Direct Gap = " << DirectGap << endl;
}
