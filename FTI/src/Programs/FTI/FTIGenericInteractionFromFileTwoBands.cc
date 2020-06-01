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
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandWithSpinHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian.h"
 
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DExplicitBandStructure.h"
#include "Tools/FTITightBinding/TightBindingModel2DExplicitBlochHamiltonian.h"

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
  (*SystemGroup) += new SingleDoubleOption  ('u', "interaction-rescaling", "global rescaling of the two-body intearction", 1.0);
  (*SystemGroup) += new SingleStringOption  ('\n', "singleparticle-file", "optional name of the file containing the one-body matrix elements");
  (*SystemGroup) += new BooleanOption  ('\n', "full-singleparticle", "the one-body matrix element file contains off-diagonal inter-band contributions");
  (*SystemGroup) += new BooleanOption  ('\n', "add-valley", "add valley-like degree of freedom (i.e. U(1) symmetry) included in --interaction-file");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pz-value", "twice the valley Pz value", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ez-value", "twice the Ez =1/2(N_{1u}+N_{2d}-N_{1d}-N_{2u}) value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "add-spin", "add spin 1/2 degree of freedom while assuming an SU(2) invariant interaction");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-value", "twice the spin Sz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-valleyspin", "use the spin per valley instead of --sz-value and --ez-value");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz1-value", "twice the Sz value in valley 1", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz2-value", "twice the Sz value in valley 2", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "conserve-bandoccuption", "assume that the interaction conserves the number of particles per band");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
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
	  if (Manager.GetBoolean("conserve-bandoccuption") == false)
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	      sprintf (CommentLine, "eigenvalues\n# kx ky");
	    }
	  else
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	      sprintf (CommentLine, "eigenvalues\n# kx ky n1 n2");
	    }
	}
      else
	{
	  if (Manager.GetBoolean("conserve-bandoccuption") == false)
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinSz);
	      sprintf (CommentLine, "eigenvalues\n# Sz kx ky");
	    }
	  else
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinSz);
	      sprintf (CommentLine, "eigenvalues\n# Sz kx ky n1up n1down n2up n2down");
	    }
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
  if (Manager.GetBoolean("flat-band"))
    {
      sprintf (FilePrefix, "%s_twoband_flatband_%s_%s", StatisticPrefix, Manager.GetString("interaction-name"), FileSystemGeometry);
    }
  else
    {
      sprintf (FilePrefix, "%s_twoband_u_%.3f_%s_%s", StatisticPrefix, Manager.GetDouble("interaction-rescaling"),
	       Manager.GetString("interaction-name"), FileSystemGeometry);
    }
  
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
  
  if (Manager.GetString("import-onebody") == 0)
    {
      if (Manager.GetString("singleparticle-file") == 0)
	{
	  double* DummyChemicalPotentials = new double[2];
	  DummyChemicalPotentials[0] = 0.0;
	  DummyChemicalPotentials[1] = 0.0;
	  TightBindingModel = new TightBindingModel2DAtomicLimitLattice (NbrSitesX, NbrSitesY, 2, DummyChemicalPotentials,
									 0.0, 0.0, Architecture.GetArchitecture(), true);
	}
      else
	{
	  double* DummyChemicalPotentials = new double[2];
	  DummyChemicalPotentials[0] = 0.0;
	  DummyChemicalPotentials[1] = 0.0;
	  int* TmpKxValues = 0;
	  int* TmpKyValues = 0;
	  int* TmpValleyIndices = 0;
	  int* TmpBandIndices1 = 0;	  
	  int* TmpBandIndices2 = 0;	  
	  double* TmpOneBodyEnergies = 0;
	  MultiColumnASCIIFile OneBodyEnergyFile;
	  if (OneBodyEnergyFile.Parse(Manager.GetString("singleparticle-file")) == false)
	    {
	      OneBodyEnergyFile.DumpErrors(cout) << endl;
	      return 0;
	    }
	  if (OneBodyEnergyFile.GetNbrLines() == 0)
	    {
	      cout << Manager.GetString("singleparticle-file") << " is an empty file" << endl;
	      return 0;
	    }
	  int NbrEnergies = OneBodyEnergyFile.GetNbrLines();
	  if (Manager.GetBoolean("full-singleparticle") == false)
	    {
	      if (((Manager.GetBoolean("add-valley") == false) && (NbrEnergies != (2 * NbrSitesX * NbrSitesY)))
		  || ((Manager.GetBoolean("add-valley") == true) && (NbrEnergies != (4 * NbrSitesX * NbrSitesY))))
		{
		  cout << Manager.GetString("singleparticle-file") << " has a wrong number of lines (has "
		       << NbrEnergies << ", should be " << (2 * NbrSitesX * NbrSitesY)
		       << " without valley, " << (4 * NbrSitesX * NbrSitesY) << " with vallley)" << endl;
		  return 0;
		}
	      if (((Manager.GetBoolean("add-valley") == false) && (OneBodyEnergyFile.GetNbrColumns() < 4))
		  || ((Manager.GetBoolean("add-valley") == true) && (OneBodyEnergyFile.GetNbrColumns() < 5)))
		{
		  cout << Manager.GetString("singleparticle-file") << " has a wrong number of columns (has "
		       << OneBodyEnergyFile.GetNbrColumns() << ", should be at least 4 without valley, 5 with valley)" << endl;
		  return 0;
		}
	    }
	  else
	    {
	      if (((Manager.GetBoolean("add-valley") == false) && (NbrEnergies != (4 * NbrSitesX * NbrSitesY)))
		  || ((Manager.GetBoolean("add-valley") == true) && (NbrEnergies != (8 * NbrSitesX * NbrSitesY))))
		{
		  cout << Manager.GetString("singleparticle-file") << " has a wrong number of lines (has "
		       << NbrEnergies << ", should be " << (4 * NbrSitesX * NbrSitesY)
		       << " without valley and --full-singleparticle, " << (8 * NbrSitesX * NbrSitesY) << " with vallley and --full-singleparticle)" << endl;
		  return 0;
		}
	      if (((Manager.GetBoolean("add-valley") == false) && (OneBodyEnergyFile.GetNbrColumns() < 5))
		  || ((Manager.GetBoolean("add-valley") == true) && (OneBodyEnergyFile.GetNbrColumns() < 6)))
		{
		  cout << Manager.GetString("singleparticle-file") << " has a wrong number of columns (has "
		       << OneBodyEnergyFile.GetNbrColumns() << ", should be at least 5 without valley, 6 with valley when using --full-singleparticle)" << endl;
		  return 0;
		}
	    }

	  int TmpNbrBands = 2;
	  if (Manager.GetBoolean("add-valley") == true)
	    {
	      TmpNbrBands = 4;
	    }
	  if (Manager.GetBoolean("add-spin") == true)
	    {
	      TmpNbrBands *= 2;
	    }
	  int* KxValues = new int[NbrSitesX * NbrSitesY];
	  int* KyValues = new int[NbrSitesX * NbrSitesY];
	  int TmpIndex = 0;
	  int* IndexToKIndex = new int[NbrSitesX * NbrSitesY];
	  int* NbrBandsPerKSector = new int[NbrSitesX * NbrSitesY];
	  for (int i = 0; i < NbrSitesX; ++i)
	    {
	      for (int j = 0; j < NbrSitesY; ++j)
		{
		  KxValues[TmpIndex] = i;
		  KyValues[TmpIndex] = j;
		  IndexToKIndex[(i * NbrSitesY) + j] = TmpIndex;
		  ++TmpIndex;
		}
	    }
	  TmpKxValues = OneBodyEnergyFile.GetAsIntegerArray(0);
	  TmpKyValues = OneBodyEnergyFile.GetAsIntegerArray(1);

	  if (Manager.GetBoolean("full-singleparticle") == false)
	    {
	      double** Energies = new double*[NbrSitesX * NbrSitesY];
	      TmpIndex = 0;
	      for (int i = 0; i < NbrSitesX; ++i)
		{
		  for (int j = 0; j < NbrSitesY; ++j)
		    {
		      Energies[TmpIndex] = new double[TmpNbrBands];
		      ++TmpIndex;
		    }
		}
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  TmpBandIndices1 = OneBodyEnergyFile.GetAsIntegerArray(2);
		  TmpOneBodyEnergies = OneBodyEnergyFile.GetAsDoubleArray(3);
		  
		  if (Manager.GetBoolean("add-spin") == false)
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][TmpBandIndices1[i]] = TmpOneBodyEnergies[i];
			}
		    }
		  else
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][TmpBandIndices1[i]] = TmpOneBodyEnergies[i];
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][2 + TmpBandIndices1[i]] = TmpOneBodyEnergies[i];
			}
		    }			      
		}
	      else
		{
		  TmpValleyIndices = OneBodyEnergyFile.GetAsIntegerArray(2);
		  TmpBandIndices1 = OneBodyEnergyFile.GetAsIntegerArray(3);
		  TmpOneBodyEnergies = OneBodyEnergyFile.GetAsDoubleArray(4);
		  if (Manager.GetBoolean("add-spin") == false)
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][TmpBandIndices1[i] + (TmpValleyIndices[i] + 1)] = TmpOneBodyEnergies[i];
			}
		    }
		  else
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][TmpBandIndices1[i] + (TmpValleyIndices[i] + 1)] = TmpOneBodyEnergies[i];
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][4 + TmpBandIndices1[i] + (TmpValleyIndices[i] + 1)] = TmpOneBodyEnergies[i];
			}
		    }			      
		}
	      TightBindingModel = new TightBindingModel2DExplicitBandStructure (NbrSitesX, NbrSitesY, TmpNbrBands, 0.0, 0.0,
										KxValues, KyValues, Energies,
										Architecture.GetArchitecture(), true);
	      TmpIndex = 0;
	      for (int i = 0; i < NbrSitesX; ++i)
		{
		  for (int j = 0; j < NbrSitesY; ++j)
		    {
		      delete[] Energies[TmpIndex];
		      ++TmpIndex;
		    }
		}
	      delete[] Energies;
	    }
	  else
	    {
	      HermitianMatrix* BlochHamiltonian = new HermitianMatrix[NbrSitesX * NbrSitesY];
	      TmpIndex = 0;
	      for (int i = 0; i < NbrSitesX; ++i)
		{
		  for (int j = 0; j < NbrSitesY; ++j)
		    {
		      BlochHamiltonian[TmpIndex] = HermitianMatrix(TmpNbrBands, true);
		      ++TmpIndex;
		    }
		}
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  TmpBandIndices1 = OneBodyEnergyFile.GetAsIntegerArray(2);
		  TmpBandIndices2 = OneBodyEnergyFile.GetAsIntegerArray(3);
		  TmpOneBodyEnergies = OneBodyEnergyFile.GetAsDoubleArray(4);
		  
		  if (Manager.GetBoolean("add-spin") == false)
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i], TmpBandIndices2[i], TmpOneBodyEnergies[i]);
			}
		    }
		  else
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i], TmpBandIndices2[i], TmpOneBodyEnergies[i]);
			  BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(2 + TmpBandIndices1[i], 2 + TmpBandIndices2[i], TmpOneBodyEnergies[i]);
			}
		    }			      
		}
	      else
		{
		  TmpBandIndices1 = OneBodyEnergyFile.GetAsIntegerArray(2);
		  TmpBandIndices2 = OneBodyEnergyFile.GetAsIntegerArray(3);
		  TmpValleyIndices = OneBodyEnergyFile.GetAsIntegerArray(4);
		  TmpOneBodyEnergies = OneBodyEnergyFile.GetAsDoubleArray(5);
		  if (Manager.GetBoolean("add-spin") == false)
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i] + (TmpValleyIndices[i] + 1), TmpBandIndices2[i] + (TmpValleyIndices[i] + 1), TmpOneBodyEnergies[i]);
			}
		    }
		  else
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i] + (TmpValleyIndices[i] + 1), TmpBandIndices2[i] + (TmpValleyIndices[i] + 1), TmpOneBodyEnergies[i]);
			  BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(4 + TmpBandIndices1[i] + (TmpValleyIndices[i] + 1), 4 + TmpBandIndices2[i] + (TmpValleyIndices[i] + 1), TmpOneBodyEnergies[i]);
			}
		    }			      
		}
	      TightBindingModel = new TightBindingModel2DExplicitBlochHamiltonian (NbrSitesX, NbrSitesY, TmpNbrBands, 0.0, 0.0,
										   KxValues, KyValues, BlochHamiltonian,
										   Architecture.GetArchitecture(), true);
	      delete[] BlochHamiltonian;
	    }
	  char* BandStructureOutputFile = new char [64 + strlen(FilePrefix)];
	  sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
	  TightBindingModel->WriteBandStructure(BandStructureOutputFile);

	  delete[] IndexToKIndex;
	  delete[] KxValues;
	  delete[] KyValues;
	}      
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

  int* NbrParticlesBand1UpPlus = 0;
  int* NbrParticlesBand2UpPlus = 0;
  int* NbrParticlesBand1UpMinus = 0;
  int* NbrParticlesBand2UpMinus = 0;
  int* NbrParticlesBand1DownPlus = 0;
  int* NbrParticlesBand2DownPlus = 0;
  int* NbrParticlesBand1DownMinus = 0;
  int* NbrParticlesBand2DownMinus = 0;
  int NbrSymmetrySectors = 1;
  if (Manager.GetBoolean("conserve-bandoccuption") == false)
    {
      if (Manager.GetBoolean("add-spin") == false)
	{
	  if (Manager.GetBoolean("add-valley") == false)
	    {
	      NbrSymmetrySectors = 1;
	      NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1UpPlus[0] = NbrParticles;
	      NbrParticlesBand2UpPlus[0] = 0;
	      NbrParticlesBand1UpMinus[0] = 0;
	      NbrParticlesBand2UpMinus[0] = 0;
	      NbrParticlesBand1DownPlus[0] = 0;
	      NbrParticlesBand2DownPlus[0] = 0;
	      NbrParticlesBand1DownMinus[0] = 0;
	      NbrParticlesBand2DownMinus[0] = 0;	      
	    }
	  else
	    {
	      NbrSymmetrySectors = 1;
	      NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1UpPlus[0] = (NbrParticles + MinPz) / 2;
	      NbrParticlesBand2UpPlus[0] = 0;
	      NbrParticlesBand1UpMinus[0] = (NbrParticles - MinPz) / 2;
	      NbrParticlesBand2UpMinus[0] = 0;
	      NbrParticlesBand1DownPlus[0] = 0;
	      NbrParticlesBand2DownPlus[0] = 0;
	      NbrParticlesBand1DownMinus[0] = 0;
	      NbrParticlesBand2DownMinus[0] = 0;	      
	    }
	}
      else
	{
	  if (Manager.GetBoolean("add-valley") == false)
	    {
	      NbrSymmetrySectors = 1;
	      NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1UpPlus[0] = (NbrParticles + MinSz) / 2;
	      NbrParticlesBand2UpPlus[0] = 0;
	      NbrParticlesBand1UpMinus[0] = 0;
	      NbrParticlesBand2UpMinus[0] = 0;
	      NbrParticlesBand1DownPlus[0] = (NbrParticles - MinSz) / 2;
	      NbrParticlesBand2DownPlus[0] = 0;
	      NbrParticlesBand1DownMinus[0] = 0;
	      NbrParticlesBand2DownMinus[0] = 0;	      
	    }
	  else
	    {
	      NbrSymmetrySectors = 1;
	      NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1UpPlus[0] = (NbrParticles + MinSz + MinPz + MinEz);
	      NbrParticlesBand1UpMinus[0] = (NbrParticles + MinSz - MinPz - MinEz);
	      NbrParticlesBand1DownPlus[0] = (NbrParticles - MinSz + MinPz - MinEz);
	      NbrParticlesBand1DownMinus[0] = (NbrParticles - MinSz - MinPz + MinEz);			  
	      NbrParticlesBand2UpPlus[0] = 0;
	      NbrParticlesBand2UpMinus[0] = 0;
	      NbrParticlesBand2DownPlus[0] = 0;
	      NbrParticlesBand2DownMinus[0] = 0;			  
	      if ((NbrParticlesBand1UpPlus[0] < 0) || (NbrParticlesBand1UpMinus[0] < 0) || (NbrParticlesBand1DownPlus[0] < 0) || (NbrParticlesBand1DownMinus[0] < 0)
		  || ((NbrParticlesBand1UpPlus[0] & 3) != 0) ||  ((NbrParticlesBand1UpMinus[0] & 3) != 0)
		  || ((NbrParticlesBand1DownPlus[0] & 3) != 0) ||  ((NbrParticlesBand1DownMinus[0] & 3) != 0))
		{
		  cout << "Incompatible values of N, 2Sz, 2Pz and 2Ez, lead to 4N_{up,+}=" << NbrParticlesBand1UpPlus[0]
		       << " 4N_{up,-}=" << NbrParticlesBand1UpMinus[0] << " 4N_{down,+}=" << NbrParticlesBand1DownPlus[0]
		       << " 4N_{down,-}=" << NbrParticlesBand1DownMinus[0] << endl;
		  return 0;
		}
	      NbrParticlesBand1UpPlus[0] /= 4;
	      NbrParticlesBand1UpMinus[0] /= 4;
	      NbrParticlesBand1DownPlus[0] /= 4;
	      NbrParticlesBand1DownMinus[0] /= 4;
	      if ((NbrParticlesBand1UpPlus[0] > NbrParticles) || (NbrParticlesBand1UpMinus[0] > NbrParticles)
		  || (NbrParticlesBand1DownPlus[0] > NbrParticles) || (NbrParticlesBand1DownMinus[0] > NbrParticles))
		{
		  cout << "Incompatible values of N, 2Sz, 2Pz and 2Ez, lead to N_{up,+}=" << NbrParticlesBand1UpPlus[0]
		       << " N_{up,-}=" << NbrParticlesBand1UpMinus[0] << " N_{down,+}=" << NbrParticlesBand1DownPlus[0]
		       << " N_{down,-}=" << NbrParticlesBand1DownMinus[0] << endl;
		  return 0;
		}
	      cout << "N_{up,+}=" << NbrParticlesBand1UpPlus[0] << " N_{up,-}=" << NbrParticlesBand1UpMinus[0]
		   << " N_{down,+}=" << NbrParticlesBand1DownPlus[0] << " N_{down,-}=" << NbrParticlesBand1DownMinus[0] << endl;
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("add-spin") == false)
	{
	  if (Manager.GetBoolean("add-valley") == false)
	    {
	      NbrSymmetrySectors = NbrParticles + 1;
	      NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
	      for (int i = 0; i < NbrSymmetrySectors; ++i)
		{
		  NbrParticlesBand1UpPlus[i] = i;
		  NbrParticlesBand2UpPlus[i] = NbrParticles - i;
		  NbrParticlesBand1UpMinus[i] = 0;
		  NbrParticlesBand2UpMinus[i] = 0;
		  NbrParticlesBand1DownPlus[i] = 0;
		  NbrParticlesBand2DownPlus[i] = 0;
		  NbrParticlesBand1DownMinus[i] = 0;
		  NbrParticlesBand2DownMinus[i] = 0;
		}
	    }
	  else
	    {
	    }
	}
      else
	{
	  if (Manager.GetBoolean("add-valley") == false)
	    {
	      NbrSymmetrySectors = 0;
	      for (int TmpUp1 = 0; TmpUp1 <= NbrParticles; ++TmpUp1)
		{
		  int TmpUp2 = ((MinSz + NbrParticles) / 2) - TmpUp1;
		  if ((TmpUp2 >= 0) && (TmpUp2 <= NbrParticles))
		    {
		      for (int TmpDown1 = 0; TmpDown1 <= NbrParticles; ++TmpDown1)
			{
			  int TmpDown2 = ((NbrParticles - MinSz) / 2) - TmpDown1;
			  if ((TmpDown2 >= 0) && (TmpDown2 <= NbrParticles))
			    {
			      ++NbrSymmetrySectors;
			    }
			}
		    }
		}
	      NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
	      NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
	      NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
	      NbrSymmetrySectors = 0;
	      for (int TmpUp1 = 0; TmpUp1 <= NbrParticles; ++TmpUp1)
		{
		  int TmpUp2 = ((MinSz + NbrParticles) / 2) - TmpUp1;
		  if ((TmpUp2 >= 0) && (TmpUp2 <= NbrParticles))
		    {
		      for (int TmpDown1 = 0; TmpDown1 <= NbrParticles; ++TmpDown1)
			{
			  int TmpDown2 = ((NbrParticles - MinSz) / 2) - TmpDown1;
			  if ((TmpDown2 >= 0) && (TmpDown2 <= NbrParticles))
			    {
			      NbrParticlesBand1UpPlus[NbrSymmetrySectors] = TmpUp1;
			      NbrParticlesBand2UpPlus[NbrSymmetrySectors] = TmpUp2;
			      NbrParticlesBand1UpMinus[NbrSymmetrySectors] = 0;
			      NbrParticlesBand2UpMinus[NbrSymmetrySectors] = 0;
			      NbrParticlesBand1DownPlus[NbrSymmetrySectors] = TmpDown1;
			      NbrParticlesBand2DownPlus[NbrSymmetrySectors] = TmpDown2;
			      NbrParticlesBand1DownMinus[NbrSymmetrySectors] = 0;
			      NbrParticlesBand2DownMinus[NbrSymmetrySectors] = 0;
			      ++NbrSymmetrySectors;
			    }
			}
		    }
		}
	    }
	  else
	    {
	    }
	}
    }

  int TotalDim = 0;
  for (int SymmetrySectorIndex = 0; SymmetrySectorIndex < NbrSymmetrySectors; ++SymmetrySectorIndex)
    {
      for (int i = MinKx; i <= MaxKx; ++i)
	{
	  for (int j = MinKy; j <= MaxKy; ++j)
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  if (Manager.GetBoolean("add-spin") == false)
		    {
		      if (Manager.GetBoolean("conserve-bandoccuption") == false)
			{
			  cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
			}
		      else
			{
			  cout << "(kx=" << i << ",ky=" << j 
			       << ",n1=" << NbrParticlesBand1UpPlus[SymmetrySectorIndex]
			       << ",n2=" << NbrParticlesBand2UpPlus[SymmetrySectorIndex] << ") : " << endl;
			}
		    }
		  else
		    {
		      if (Manager.GetBoolean("conserve-bandoccuption") == false)
			{
			  cout << "(kx=" << i << ",ky=" << j << ",2sz=" << MinSz << ") : " << endl;
			}
		      else
			{
			  cout << "(kx=" << i << ",ky=" << j 
			       << ",n1up=" << NbrParticlesBand1UpPlus[SymmetrySectorIndex]
			       << ",n2up=" << NbrParticlesBand2UpPlus[SymmetrySectorIndex]
			       << ",n1down=" << NbrParticlesBand1DownPlus[SymmetrySectorIndex]
			       << ",n2down=" << NbrParticlesBand2DownPlus[SymmetrySectorIndex] << ") : " << endl;
			}
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("add-spin") == false)
		    {
		      if (Manager.GetBoolean("conserve-bandoccuption") == false)
			{
			  cout << "(kx=" << i << ",ky=" << j << ",2pz=" << MinPz << ") : " << endl;
			}
		      else
			{
			  cout << "(kx=" << i << ",ky=" << j << ",2pz=" << MinPz
			       << ",n1plus=" << NbrParticlesBand1UpPlus[SymmetrySectorIndex]
			       << ",n2plus=" << NbrParticlesBand2UpPlus[SymmetrySectorIndex]
			       << ",n1minus=" << NbrParticlesBand1UpMinus[SymmetrySectorIndex]
			       << ",n2minus=" << NbrParticlesBand2UpMinus[SymmetrySectorIndex] << ") : " << endl;
			}
		    }
		  else
		    {
		      if (Manager.GetBoolean("conserve-bandoccuption") == false)
			{
			  cout << "(kx=" << i << ",ky=" << j << ",2pz=" << MinPz << ",2sz=" << MinSz << ",2ez=" << MinEz<< ") : " << endl;
			}
		      else
			{
			  cout << "(kx=" << i << ",ky=" << j << ",2pz=" << MinPz
			       << ",n1upplus=" << NbrParticlesBand1UpPlus[SymmetrySectorIndex]
			       << ",n2upplus=" << NbrParticlesBand2UpPlus[SymmetrySectorIndex]
			       << ",n1upminus=" << NbrParticlesBand1UpMinus[SymmetrySectorIndex]
			       << ",n2upminus=" << NbrParticlesBand2UpMinus[SymmetrySectorIndex]
			       << ",n1downplus=" << NbrParticlesBand1DownPlus[SymmetrySectorIndex]
			       << ",n2downplus=" << NbrParticlesBand2DownPlus[SymmetrySectorIndex]
			       << ",n1downminus=" << NbrParticlesBand1DownMinus[SymmetrySectorIndex]
			       << ",n2downminus=" << NbrParticlesBand2DownMinus[SymmetrySectorIndex] << ") : " << endl;
			}
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
			  if (Manager.GetBoolean("conserve-bandoccuption") == false)
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
			      if ((NbrSitesX * NbrSitesY) <= 32)
				{
				  Space = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrParticlesBand1UpPlus[SymmetrySectorIndex], NbrSitesX, NbrSitesY, i, j);
				}
			      else
				{
				  Space = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrParticlesBand1UpPlus[SymmetrySectorIndex], NbrSitesX, NbrSitesY, i, j);
				}
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
			  if (Manager.GetBoolean("conserve-bandoccuption") == false)
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
			      int FakePz = (NbrParticlesBand1UpPlus[SymmetrySectorIndex] - NbrParticlesBand2UpPlus[SymmetrySectorIndex]
					    + NbrParticlesBand1DownPlus[SymmetrySectorIndex] - NbrParticlesBand2DownPlus[SymmetrySectorIndex]);
			      if ((NbrSitesX * NbrSitesY) <= 16)
				{
				  Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j,
											      MinSz, FakePz, 10000000ul);
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
			  if (Manager.GetBoolean("conserve-bandoccuption") == false)
			    {
			      if ((NbrSitesX * NbrSitesY) <= 8)
				{
				  Space = new FermionOnSquareLatticeWithSU8SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j,
											      NbrParticlesBand1DownMinus[SymmetrySectorIndex], NbrParticlesBand1DownPlus[SymmetrySectorIndex],
											      NbrParticlesBand1UpMinus[SymmetrySectorIndex], NbrParticlesBand1UpPlus[SymmetrySectorIndex], 10000000ul);
				}
			      else
				{
				  if ((NbrSitesX * NbrSitesY) <= 16)
				    {			  
				      Space = new FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong(NbrParticles, NbrSitesX, NbrSitesY, i, j,
												     NbrParticlesBand1DownMinus[SymmetrySectorIndex], NbrParticlesBand1DownPlus[SymmetrySectorIndex],
												     NbrParticlesBand1UpMinus[SymmetrySectorIndex], NbrParticlesBand1UpPlus[SymmetrySectorIndex], 10000000ul);
				    }
				  else
				    {
				      cout << "SU(8) not supported with more than 16 momenta" << endl;
				      Space = 0;			      
				    }
				}
			    }
			  else
			    {
			      int TmpNbrParticles[8];
			      TmpNbrParticles[0] = NbrParticlesBand1DownMinus[SymmetrySectorIndex];
			      TmpNbrParticles[1] = NbrParticlesBand2DownMinus[SymmetrySectorIndex];
			      TmpNbrParticles[2] = NbrParticlesBand1DownPlus[SymmetrySectorIndex];
			      TmpNbrParticles[3] = NbrParticlesBand2DownPlus[SymmetrySectorIndex];
			      TmpNbrParticles[4] = NbrParticlesBand1UpMinus[SymmetrySectorIndex];
			      TmpNbrParticles[5] = NbrParticlesBand2UpMinus[SymmetrySectorIndex];
			      TmpNbrParticles[6] = NbrParticlesBand1UpPlus[SymmetrySectorIndex];
			      TmpNbrParticles[7] = NbrParticlesBand2UpPlus[SymmetrySectorIndex];
			      if ((NbrSitesX * NbrSitesY) <= 8)
				{
				  Space = new FermionOnSquareLatticeWithSU8SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j, TmpNbrParticles, 10000000ul);
				}
			      else
				{
				  if ((NbrSitesX * NbrSitesY) <= 16)
				    {			  
				      Space = new FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong(NbrParticles, NbrSitesX, NbrSitesY, i, j, TmpNbrParticles, 10000000ul);
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
		}
	      else
		{
		  Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      TotalDim += Space->GetHilbertSpaceDimension();
	      if (Space->GetHilbertSpaceDimension() > 0)
		{
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		  
		  if (Manager.GetBoolean("real-interaction"))
		    {
		      if (Manager.GetBoolean("add-valley") == false)
			{
			  Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
													Manager.GetString("interaction-file"),
													TightBindingModel, Manager.GetBoolean("flat-band"), 
													Manager.GetDouble("interaction-rescaling"),
													Manager.GetBoolean("add-spin"),
													Architecture.GetArchitecture(), Memory);
			}
		      else
			{
			  Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
														Manager.GetString("interaction-file"),
														TightBindingModel, Manager.GetBoolean("flat-band"), 
														Manager.GetDouble("interaction-rescaling"),
														Manager.GetBoolean("add-spin"),
														Architecture.GetArchitecture(), Memory);
			}		
		    }
		  else
		    {
		      if (Manager.GetBoolean("add-valley") == false)
			{
			  Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
												    Manager.GetString("interaction-file"),
												    TightBindingModel, Manager.GetBoolean("flat-band"), 
												    Manager.GetDouble("interaction-rescaling"),
												    Manager.GetBoolean("add-spin"),
												    Architecture.GetArchitecture(), Memory);
			}
		      else
			{
			  Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandWithSpinHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
													    Manager.GetString("interaction-file"),
													    TightBindingModel, Manager.GetBoolean("flat-band"), 
													    Manager.GetDouble("interaction-rescaling"),
													    Manager.GetBoolean("add-spin"),
													    Architecture.GetArchitecture(), Memory);
			}				
		    }
		  
		  
		  char* ContentPrefix = new char[256];
		  char* EigenstateOutputFile = new char [512];
		  char* TmpExtention = new char[256];
		  if (Manager.GetBoolean("add-valley") == false)
		    {
		      if (Manager.GetBoolean("add-spin") == false)
			{
			  if (Manager.GetBoolean("conserve-bandoccuption") == false)
			    {
			      sprintf (ContentPrefix, "%d %d", i, j);
			      sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
			    }
			  else
			    {
			      sprintf (ContentPrefix, "%d %d %d %d", i, j, NbrParticlesBand1UpPlus[SymmetrySectorIndex],
				       NbrParticlesBand2UpPlus[SymmetrySectorIndex]);
			      sprintf (TmpExtention, "_kx_%d_ky_%d_bz_%d", i, j, (NbrParticlesBand1UpPlus[SymmetrySectorIndex] - NbrParticlesBand2UpPlus[SymmetrySectorIndex]));
			    }
			}
		      else
			{
			  if (Manager.GetBoolean("conserve-bandoccuption") == false)
			    {
			      sprintf (ContentPrefix, "%d %d %d", MinSz, i, j);
			      sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d", i, j, MinSz);
			    }
			  else
			    {
			      sprintf (ContentPrefix, "%d %d %d %d %d %d %d", MinSz, i, j, NbrParticlesBand1UpPlus[SymmetrySectorIndex],
				       NbrParticlesBand2UpPlus[SymmetrySectorIndex], NbrParticlesBand1DownPlus[SymmetrySectorIndex],
				       NbrParticlesBand2DownPlus[SymmetrySectorIndex]);
			      sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d", i, j, MinSz);
			    }
			}
		    }
		  else
		    {
		      if (Manager.GetBoolean("add-spin") == false)
			{
			  if (Manager.GetBoolean("conserve-bandoccuption") == false)
			    {
			      sprintf (ContentPrefix, "%d %d %d", MinPz, i, j);
			      sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d", i, j, MinPz);
			    }
			  else
			    {
			    }
			}
		      else
			{
			  if (Manager.GetBoolean("conserve-bandoccuption") == false)
			    {
			      sprintf (ContentPrefix, "%d %d %d %d %d", MinPz, MinSz, MinEz, i, j);
			      sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d_ez_%d_sz_%d", i, j, MinPz, MinEz, MinSz);
			    }
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
		  delete[] EigenstateOutputFile;
		  delete[] ContentPrefix;
		}
	      delete Space;
	    }
	}
    }
  cout << "Total dim=" << TotalDim << endl;
  return 0;
}

