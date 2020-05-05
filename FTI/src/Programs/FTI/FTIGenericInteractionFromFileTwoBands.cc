#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian.h"
 
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
  sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY)  ;
  char* FilePrefix = new char [512 + strlen(FileSystemGeometry)];
  sprintf (FilePrefix, "%s_%s_%s", StatisticPrefix, Manager.GetString("interaction-name"), FileSystemGeometry);

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky");
  
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
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
	  ParticleOnSphereWithSpin* Space = 0;
	  AbstractQHEHamiltonian* Hamiltonian = 0;
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  if (Manager.GetBoolean("boson") == false)
	    {
	      if ((NbrSitesX * NbrSitesY) <= 31)
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
	      Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
	    }
	  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

	  if (Manager.GetBoolean("real-interaction"))
	    {
	      Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
											    Manager.GetString("interaction-file"),
											    TightBindingModel, Manager.GetBoolean("flat-band"), 
											    Manager.GetDouble("flatband-gap"),
											    Architecture.GetArchitecture(), Memory);
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
	  sprintf (ContentPrefix, "%d %d", i, j);
	  char* EigenstateOutputFile = new char [512];
	  char* TmpExtention = new char[256];
	  sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
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

