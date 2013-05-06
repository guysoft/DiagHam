#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

#include "Tools/FTITightBinding/TightBindingModelTimeReversalKagomeLattice.h"

#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "GeneralTools/FilenameTools.h"

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
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
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
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mixing-12", "mixing term coupling the two copies of the kagome lattice (sites 1 and 2)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mixing-13", "mixing term coupling the two copies of the kagome lattice (sites 1 and 3)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mixing-23", "mixing term coupling the two copies of the kagome lattice (sites 2 and 3)", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-z2invariant", "compute the z2 invariant of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytheta", "export the one-body topological information (phase of the eigenvalues of the D matrix) in an ASCII text file");
  (*SystemGroup) += new BooleanOption ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new BooleanOption ('\n', "decoupled", "assume two decoupled copies of the kagome lattice");
  (*SystemGroup) += new BooleanOption ('\n', "fixed-sz", "fix the Sz value when considering two decoupled copies of the checkerboard lattice");
  (*SystemGroup) += new SingleIntegerOption ('\n', "sz-value", "twice the fixed Sz value", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-hamiltonian", "export the hamiltonian in a column formatted ASCII file");
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
  if (Manager.GetBoolean("decoupled") == false)
    {
      sprintf (CommentLine, "eigenvalues\n# kx ky E");
    }
  else
    {
      sprintf (CommentLine, "eigenvalues\n# kx ky Sz E");
    }

  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetDouble("mu-s") == 0.0)
    {
      if (Manager.GetBoolean("decoupled") == true)
	{
	  sprintf (EigenvalueOutputFile, "%s_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_u_%f_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f.dat", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	}
      else
	{
	  sprintf (EigenvalueOutputFile, "%s_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_u_%f_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_mix12_%f_mix13_%f_mix23_%f_gx_%f_gy_%f.dat", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mixing-12"), Manager.GetDouble("mixing-13"), Manager.GetDouble("mixing-23"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	}
    }
  else
    {
      if (Manager.GetBoolean("decoupled") == true)
	{
	  sprintf (EigenvalueOutputFile, "%s_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_u_%f_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_gx_%f_gy_%f_mus_%f.dat", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	}
      else
	{
	  sprintf (EigenvalueOutputFile, "%s_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_u_%f_v_%f_w_%f_t1_%f_t2_%f_l1_%f_l2_%f_mix12_%f_mix13_%f_mix23_%f_gx_%f_gy_%f_mus_%f.dat", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mixing-12"), Manager.GetDouble("mixing-13"), Manager.GetDouble("mixing-23"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	}
    }
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || 
	  (Manager.GetBoolean("singleparticle-z2invariant") == true) || (Manager.GetBoolean("export-onebodytheta") == true))
	ExportOneBody = true;
      TightBindingModelTimeReversalKagomeLattice TightBindingModel(NbrSitesX, NbrSitesY,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), 
								   Manager.GetDouble("l1"), Manager.GetDouble("l2"), 
								   Manager.GetDouble("mixing-12"), Manager.GetDouble("mixing-13"), Manager.GetDouble("mixing-23"),
								   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);
      /*
      TightBindingModel.WriteAsciiSpectrum(EigenvalueOutputFile);*/
//       cout << "Chern number = " << TightBindingModel.ComputeChernNumber(0) << endl;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	{
	  char* BandStructureOutputFile = new char [512];
	  if (Manager.GetString("export-onebodyname") != 0)
	    strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	  else
	    sprintf (BandStructureOutputFile, "%s_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_t1_%f_t2_%f_l1_%f_l2_%f_mix12_%f_mix13_%f_mix23_%f_gx_%f_gy_%f_tightbinding.dat", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mixing-12"), Manager.GetDouble("mixing-13"), Manager.GetDouble("mixing-23"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	  if (Manager.GetBoolean("export-onebody") == true)
	    {
	      TightBindingModel.WriteBandStructure(BandStructureOutputFile);
	    }
	  else
	    {
	      TightBindingModel.WriteBandStructureASCII(BandStructureOutputFile);
	    }
	  delete[] BandStructureOutputFile;
	}
      if (Manager.GetBoolean("singleparticle-z2invariant") == true)
	cout << "Z2 invariant = " << TightBindingModel.ComputeZ2Invariant(2) << endl;
      
      if (Manager.GetBoolean("export-onebodytheta") == true)
      {
	cout << "Z2 invariant = " << TightBindingModel.ComputeZ2Invariant(2) << endl;
	char* ThetaOutputFile = new char [512];
	sprintf(ThetaOutputFile, "%s_twoband_quantumspinhall_kagome_n_%d_x_%d_y_%d_t1_%f_t2_%f_l1_%f_l2_%f_mix12_%f_mix13_%f_mix23_%f_gx_%f_gy_%f_theta.dat", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mixing-12"), Manager.GetDouble("mixing-13"), Manager.GetDouble("mixing-23"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	TightBindingModel.WriteAsciiDMatrixEigenValues(ThetaOutputFile, 2);
      }
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
  
  TightBindingModelTimeReversalKagomeLattice *TightBindingModel;
  TightBindingModel = new TightBindingModelTimeReversalKagomeLattice(NbrSitesX, NbrSitesY,  
								   Manager.GetDouble("t1"), Manager.GetDouble("t2"), 
								   Manager.GetDouble("l1"), Manager.GetDouble("l2"), 
								   Manager.GetDouble("mixing-12"), Manager.GetDouble("mixing-13"), Manager.GetDouble("mixing-23"),
								   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true);
  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  if (Manager.GetBoolean("decoupled") == false)
	    {
	      cout << "(kx=" << i << ",ky=" << j << ") " << endl;
	      ParticleOnSphereWithSpin* Space = 0;
	      if (Manager.GetBoolean("boson") == false)
		{
		  Space = new FermionOnSquareLatticeWithSpinMomentumSpace(NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      else
		{
		  Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace(NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
														 Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"),TightBindingModel, 		     
														 Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", i, j);
	      char* EigenstateOutputFile = new char [512];
	      char* TmpExtention = new char [512];
	      sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
	      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	      FirstRunFlag = false;
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      cout << "------------------------------------" << endl;
	      delete Hamiltonian;
	      delete[] EigenstateOutputFile;
	      delete[] ContentPrefix;
	      delete Space;
	    }
	  else
	    {
	      int MinSz = -NbrParticles;
	      int MaxSz = NbrParticles;
	      if (Manager.GetBoolean("fixed-sz") == true)
		{
		  MinSz = Manager.GetInteger("sz-value");
		  MaxSz = MinSz;
		}
	      for (int Sz = MinSz; Sz <= MaxSz; Sz += 2)
		{
		  cout << "(kx=" << i << ",ky=" << j << ") Sz=" << Sz << " : " << endl;
		  ParticleOnSphereWithSpin* Space = 0;
		  if (Manager.GetBoolean("boson") == false)
		    {
		      Space = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, (Sz + NbrParticles) / 2, NbrSitesX, NbrSitesY, i, j);
		    }
		  else
		    {
		      Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, (Sz + NbrParticles) / 2, NbrSitesX, NbrSitesY, i, j);
		    }
		  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
		  AbstractQHEHamiltonian* Hamiltonian = 0;
		  Hamiltonian = new ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
												      Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
												      Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"),
												      Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
												      Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  char* ContentPrefix = new char[256];
		  sprintf (ContentPrefix, "%d %d %d", i, j, Sz);
		  char* EigenstateOutputFile = new char [512];
		  char* TmpExtention = new char [512];
		  sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d", i, j, Sz);
		  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
		  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
		  FirstRunFlag = false;
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  cout << "------------------------------------" << endl;
		  delete Hamiltonian;
		  delete[] EigenstateOutputFile;
		  delete[] ContentPrefix;
		  delete Space;
		}
	    }
	}
    }
  return 0;
}

