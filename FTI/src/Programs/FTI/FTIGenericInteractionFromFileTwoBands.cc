#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU8SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian.h"
 
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("FTIGenericInteractionFromFileTwoBands" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new SingleStringOption  ('\n', "interaction-file", "name of the file containing the two-body interaction matrix elements");
  (*SystemGroup) += new BooleanOption  ('\n', "real-interaction", "assume that the two-body interaction matrix elements are real");
  (*SystemGroup) += new SingleStringOption  ('\n', "interaction-name", "name of the two-body interaction", "noname");
  (*SystemGroup) += new BooleanOption  ('\n', "add-valley", "add valley-like degree of freedom (i.e. U(1) symmetry) included in --interaction-file");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pz-value", "twice the valley Pz value", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ez-value", "twice the Ez =1/2(N_{1u}+N_{2d}-N_{1d}-N_{2u}) value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "add-spin", "add spin 1/2 degree of freedom while assuming an SU(2) invariant interaction");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-value", "twice the spin Sz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-valleyspin", "use the spin per valley instead of --sz-value and --ez-value");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz1-value", "twice the Sz value in valley 1", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz2-value", "twice the Sz value in valley 2", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "flatband-gap", "when using the flat band model with two bands, set the one-body gap between the two bands", 0.0);
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
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
  (*ToolsGroup) += new SingleDoubleOption  ('\n',"testhermitian-error", "precision of the hermeticity test",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIGenericInteractionFromFileTwoBands -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "no interaction file defined" << endl;
      cout << "see man page for option syntax or type FTIGenericInteractionFromFileTwoBands -h" << endl;
      return -1;
    }
  if (!(IsFile(Manager.GetString("interaction-file"))))
    {
      cout << "interaction file " << Manager.GetString("interaction-file")<< " does not exist" << endl;
      return -1;
    }
  
    
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSites = 2 * NbrSitesX * NbrSitesY;
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
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
  int MinSz = NbrParticles & 1;
  int MaxSz = MinSz;
  if (Manager.GetBoolean("add-spin") == true)
    {
      MinSz = Manager.GetInteger("sz-value") | (NbrParticles & 1);
      if (Manager.GetBoolean("use-valleyspin") == true)
	{
	  MinSz = (Manager.GetInteger("sz1-value") + Manager.GetInteger("sz2-value")) | (NbrParticles & 1);
	}
      MaxSz = MinSz;
    }  
  int MinPz = NbrParticles & 1;
  int MaxPz = MinPz;
  if (Manager.GetBoolean("add-valley") == true)
    {
      MinPz = Manager.GetInteger("pz-value") | (NbrParticles & 1);
      MaxPz = MinPz;
    }  
  int MinEz = NbrParticles & 1;
  int MaxEz = MinEz;
  if ((Manager.GetBoolean("add-valley") == true) && (Manager.GetBoolean("add-spin") == true))
    {
      MinEz = Manager.GetInteger("ez-value") | (NbrParticles & 1);
      if (Manager.GetBoolean("use-valleyspin") == true)
	{
	  MinEz = (Manager.GetInteger("sz1-value") - Manager.GetInteger("sz2-value")) | (NbrParticles & 1);
	}
      MaxEz = MinEz;
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


  char* FileSystemGeometry = new char [512];
  char* CommentLine = new char [256];
  if (Manager.GetBoolean("add-valley") == false)
    {
      if (Manager.GetBoolean("add-spin") == false)
	{
	  sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	  sprintf (CommentLine, "eigenvalues\n# kx ky");
	}
      else
	{
	  sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinSz);
	  sprintf (CommentLine, "eigenvalues\n# Sz kx ky");
	}
   }
  else
    {
      if (Manager.GetBoolean("add-spin") == false)
	{
	  sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinPz);
	  sprintf (CommentLine, "eigenvalues\n# Pz kx ky");
	}
      else
	{
	  sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d_ez_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinPz, MinEz, MinSz);
	  sprintf (CommentLine, "eigenvalues\n# Pz Sz Ez kx ky");
	}
    }
  char* FilePrefix = new char [512 + strlen(FileSystemGeometry)];
  sprintf (FilePrefix, "%s_%s_%s", StatisticPrefix, Manager.GetString("interaction-name"), FileSystemGeometry);
  
  char* EigenvalueOutputFile = new char [512 + strlen(FilePrefix)];
  
  if (Manager.GetString("eigenvalue-file") != 0)
    {
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
    }
  else
    {
      sprintf (EigenvalueOutputFile, "%s.dat", FilePrefix);
    }
  
  Abstract2DTightBindingModel* TightBindingModel;
  double* DummyChemicalPotentials = new double[2];
  DummyChemicalPotentials[0] = 0.0;
  DummyChemicalPotentials[1] = 0.0;
  
  if (Manager.GetString("import-onebody") == 0)
    {
      TightBindingModel = new TightBindingModel2DAtomicLimitLattice (NbrSitesX, NbrSitesY, 2, DummyChemicalPotentials,
								     0.0, 0.0, Architecture.GetArchitecture(), true);
      
      // char* BandStructureOutputFile = new char [64 + strlen(FilePrefix) + strlen(FileParameterString) + strlen(FileTwistedBoundaryConditions)];
      // sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString, FileTwistedBoundaryConditions);
      // TightBindingModel->WriteBandStructure(BandStructureOutputFile);
    }
  else
    {
      TightBindingModel = new Generic2DTightBindingModel(Manager.GetString("import-onebody")); 
    }

  bool FirstRunFlag = true;

  if (Manager.GetBoolean("real-interaction"))
    {
      Lanczos.SetRealAlgorithms();
    }
  
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  if (Manager.GetBoolean("add-valley") == false)
	    {
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
		}
	      else
		{
		  cout << "(kx=" << i << ",ky=" << j << ",2sz=" << MinSz << ") : " << endl;
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  cout << "(kx=" << i << ",ky=" << j << ",2pz" << MinPz << ") : " << endl;
		}
	      else
		{
		  cout << "(kx=" << i << ",ky=" << j << ",2pz" << MinPz << ",2sz=" << MinSz << ",2ez=" << MinEz<< ") : " << endl;
		}
	    }
	  ParticleOnSphereWithSpin* Space = 0;
	  AbstractQHEHamiltonian* Hamiltonian = 0;
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  if (Manager.GetBoolean("boson") == false)
	    {
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  if (Manager.GetBoolean("add-valley") == false)
		    {
		      if ((NbrSitesX * NbrSitesY) <= 32)
			{
			  Space = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
			}
		      else
			{
			  Space = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
			}
		    }
		  else
		    {
		      if ((NbrSitesX * NbrSitesY) <= 16)
			{
			  Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j,
										      MinPz, 10000000ul);
			}
		      else
			{
			  cout << "SU(4) not supported with more than 16 momenta" << endl;
			  Space = 0;
			  //		      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j, MinSz);
			}
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("add-valley") == false)
		    {
		      if ((NbrSitesX * NbrSitesY) <= 16)
			{
			  Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j,
										      MinSz, 10000000ul);
			}
		      else
			{
			  cout << "SU(4) not supported with more than 16 momenta" << endl;
			  Space = 0;
			  //		      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j, MinSz);
			}
		    }
		  else
		    {
		      int NbrParticlesUpPlus = (NbrParticles + MinSz + MinPz + MinEz);
		      int NbrParticlesUpMinus = (NbrParticles + MinSz - MinPz - MinEz);
		      int NbrParticlesDownPlus = (NbrParticles - MinSz + MinPz - MinEz);
		      int NbrParticlesDownMinus = (NbrParticles - MinSz - MinPz + MinEz);			  
		      if ((NbrParticlesUpPlus < 0) || (NbrParticlesUpMinus < 0) || (NbrParticlesDownPlus < 0) || (NbrParticlesDownMinus < 0)
			  || ((NbrParticlesUpPlus & 3) != 0) ||  ((NbrParticlesUpMinus & 3) != 0)
			  || ((NbrParticlesDownPlus & 3) != 0) ||  ((NbrParticlesDownMinus & 3) != 0))
			{
			  cout << "Incompatible values of N, 2Sz, 2Pz and 2Ez, lead to 4N_{up,+}=" << NbrParticlesUpPlus
				   << " 4N_{up,-}=" << NbrParticlesUpMinus << " 4N_{down,+}=" << NbrParticlesDownPlus
			       << " 4N_{down,-}=" << NbrParticlesDownMinus << endl;
			  return 0;
			    }
		      NbrParticlesUpPlus /= 4;
		      NbrParticlesUpMinus /= 4;
		      NbrParticlesDownPlus /= 4;
		      NbrParticlesDownMinus /= 4;
		      if ((NbrParticlesUpPlus > NbrParticles) || (NbrParticlesUpMinus > NbrParticles)
			  || (NbrParticlesDownPlus > NbrParticles) || (NbrParticlesDownMinus > NbrParticles))
			{
			  cout << "Incompatible values of N, 2Sz, 2Pz and 2Ez, lead to N_{up,+}=" << NbrParticlesUpPlus
			       << " N_{up,-}=" << NbrParticlesUpMinus << " N_{down,+}=" << NbrParticlesDownPlus
			       << " N_{down,-}=" << NbrParticlesDownMinus << endl;
			  return 0;
			}
		      cout << "N_{up,+}=" << NbrParticlesUpPlus << " N_{up,-}=" << NbrParticlesUpMinus
			   << " N_{down,+}=" << NbrParticlesDownPlus << " N_{down,-}=" << NbrParticlesDownMinus << endl;
		      if ((NbrSitesX * NbrSitesY) <= 8)
			{
			  Space = new FermionOnSquareLatticeWithSU8SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j,
										      NbrParticlesDownMinus, NbrParticlesDownPlus,
										      NbrParticlesUpMinus, NbrParticlesUpPlus, 10000000ul);
			}
		      else
			{
			  if ((NbrSitesX * NbrSitesY) <= 16)
			    {			  
			      Space = new FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong(NbrParticles, NbrSitesX, NbrSitesY, i, j,
										      NbrParticlesDownMinus, NbrParticlesDownPlus,
										      NbrParticlesUpMinus, NbrParticlesUpPlus, 10000000ul);
			    }
			  else
			    {
			      cout << "SU(8) not supported with more than 16 momenta" << endl;
			      Space = 0;			      
			    }
			}
		    }
		}
	    }
	  else
	    {
	      Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
	    }
	  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

	  if (Manager.GetBoolean("real-interaction"))
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
												Manager.GetString("interaction-file"),
												TightBindingModel, Manager.GetBoolean("flat-band"), 
												Manager.GetDouble("flatband-gap"),
												Manager.GetBoolean("add-spin"),
												Architecture.GetArchitecture(), Memory);
		}
	      else
		{
		  Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
													Manager.GetString("interaction-file"),
													TightBindingModel, Manager.GetBoolean("flat-band"), 
													Manager.GetDouble("flatband-gap"),
													Manager.GetBoolean("add-spin"),
													Architecture.GetArchitecture(), Memory);
		}		
	    }
	  else
	    {
	      Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
											Manager.GetString("interaction-file"),
											TightBindingModel, Manager.GetBoolean("flat-band"), 
											Manager.GetDouble("flatband-gap"),
											Architecture.GetArchitecture(), Memory);
	    }
	  
	  
	  char* ContentPrefix = new char[256];
	  char* EigenstateOutputFile = new char [512];
	  char* TmpExtention = new char[256];
	  if (Manager.GetBoolean("add-valley") == false)
	    {
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  sprintf (ContentPrefix, "%d %d", i, j);
		  sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
		}
	      else
		{
		  sprintf (ContentPrefix, "%d %d %d", MinSz, i, j);
		  sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d", i, j, MinSz);
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  sprintf (ContentPrefix, "%d %d %d", MinPz, i, j);
		  sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d", i, j, MinPz);
		}
	      else
		{
		  sprintf (ContentPrefix, "%d %d %d %d %d", MinPz, MinSz, MinEz, i, j);
		  sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d_ez_%d_sz_%d", i, j, MinPz, MinEz, MinSz);
		}
	    }
	  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	  delete[] TmpExtention;
	  if (Manager.GetBoolean("real-interaction"))
	    {
	      GenericRealMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix,
				       CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	    }
	  else
	    {
	      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix,
					  CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	    }
	  FirstRunFlag = false;
	  
	  cout << "------------------------------------" << endl;
	  delete Hamiltonian;
	  delete Space;
	  delete[] EigenstateOutputFile;
	  delete[] ContentPrefix;
	}
    }
  return 0;
}

