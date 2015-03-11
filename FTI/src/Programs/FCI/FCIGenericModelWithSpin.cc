#include "Options/Options.h"

#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Tools/FTITightBinding/TightBindingModelCheckerboardLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"

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


// compute the description of the density-density interaction for the unit cell at the origin for a Checkerboard lattice model
//
// nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsPotentials = intensity of each density-density term 
// bosonFlag = true if we are dealing with bosons
// uPotential = nearest neighbor (for fermions) or on-site (for bosons) interaction amplitude
// vPotential = next nearest neighbor (for fermions) or nearest neighbor (for bosons) interaction amplitude
void FCICheckerboardLatticeModelComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
							   int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
							   bool bosonFlag, double uPotential, double vPotential);


int main(int argc, char** argv)
{
  OptionManager Manager ("FCIGenericModelWithSpin" , "0.01");
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
  (*SystemGroup) += new BooleanOption ('\n', "fixed-sz", "fix the Sz value");
  (*SystemGroup) += new SingleIntegerOption ('\n', "sz-value", "twice the fixed Sz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new SingleStringOption  ('\n', "model-name", "name of the one-body model", "unknown");
  (*SystemGroup) += new BooleanOption  ('\n', "checkerboard", "use checkerboard tightbiding model for both spins");
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the gutzwiller projected Hilbert space");
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "vupup-potential", "repulsive nearest nearest neighbor potential strength between up spins", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "vdowndown-potential", "repulsive nearest neighbor potential strength between down spins", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "vupdown-potential", "repulsive nearest neighbor potential strength between up and down spins", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "spin-flux", "adiabatic flux insertion corresponds to opposite phases for up and down spins");
  
  (*SystemGroup) += new BooleanOption  ('\n', "atomic-limit", "compute the spectrum in the atomic limit");
  
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebodyup", "import information on the spin up tight binding model from a file");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebodydown", "import information on the spin down tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "flatband-gap", "when using the flat band model with two bands, set the one-body gap between the two bands", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "real-space", "use the real space representation when considering the system with all bands");
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "use the real space representation when considering the system with all bands without the translations");
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
      cout << "see man page for option syntax or type FCIGenericModelWithSpin -h" << endl;
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
  int NbrSites = 2*NbrSitesX * NbrSitesY;
  int SzValue = Manager.GetInteger("sz-value");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  bool CheckerboardFlag = Manager.GetBoolean("checkerboard");

  char* StatisticPrefix = new char [64];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }
  
  char* FilePrefix = new char [512];
  if (Manager.GetBoolean("single-band") == false)
    {
      if (Manager.GetBoolean("real-space") == false)
	{
	  sprintf (FilePrefix, "%s_checkerboardlatticewithspin_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	}
      else
	{
         if ( Manager.GetBoolean("no-translation") == false)
	   if ( Manager.GetBoolean ("gutzwiller") == false)
	    sprintf (FilePrefix, "%s_realspace_checkerboardlatticewithspin_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	   else
	     sprintf (FilePrefix, "%s_realspace_gutzwiller_checkerboardlatticewithspin_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	else
	  if ( Manager.GetBoolean ("gutzwiller") == false)
	    sprintf (FilePrefix, "%s_realspace_notranslation_checkerboardlatticewithspin_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	  else
	    sprintf (FilePrefix, "%s_realspace_gutzwiller_notranslation_checkerboardlatticewithspin_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	}
    }
  else
  {
    cout << "Single band not implemented" << endl;
    return -1;
  }
    
  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* FileParameterString = new char [256];
//   if ((Manager.GetDouble("gamma-x") == 0.0) && (Manager.GetDouble("gamma-y") == 0.0))
//     sprintf (FileParameterString, "t1_%f_t2_%f_tpp_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"));
//   else
  if (Manager.GetBoolean("spin-flux") == false)
    sprintf (FileParameterString, "t1_%f_t2_%f_tpp_%f_gx_%f_gy_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
  else
    sprintf (FileParameterString, "t1_%f_t2_%f_tpp_%f_spingx_%f_spingy_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));

  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    {
      if (Manager.GetBoolean("single-band") == false)
	{
	  if (Manager.GetDouble("mu-s") == 0.0)
	  {
	    if (Manager.GetBoolean("fixed-sz") == false)
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%f_vuu_%f_vdd_%f_vud_%f_%s.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("vupup-potential"), Manager.GetDouble("vdowndown-potential"),Manager.GetDouble("vupdown-potential"),FileParameterString);
	    }
	    else
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%f_vuu_%f_vdd_%f_vud_%f_%s_sz_%d.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("vupup-potential"), Manager.GetDouble("vdowndown-potential"),Manager.GetDouble("vupdown-potential"),FileParameterString, SzValue);
	    }
	  }
	  else
	  {
	    if (Manager.GetBoolean("fixed-sz") == false)
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%f_vuu_%f_vdd_%f_vud_%f_%s_mus_%f.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("vupup-potential"), Manager.GetDouble("vdowndown-potential"), Manager.GetDouble("vupdown-potential"), FileParameterString, Manager.GetDouble("mu-s"));
	    }
	    else
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%f_vuu_%f_vdd_%f_vud_%f_%s_sz_%d_mus_%f.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("vupup-potential"), Manager.GetDouble("vdowndown-potential"), Manager.GetDouble("vupdown-potential"), FileParameterString, SzValue, Manager.GetDouble("mu-s"));
	    }
	  }
	}
      else
	{
	  if (Manager.GetBoolean("flat-band") == true)
	    {
	      if (Manager.GetDouble("mu-s") == 0.0)
	      {
		if (Manager.GetBoolean("fixed-sz") == false)
		{
		  if (Manager.GetBoolean("spin-flux") == false)
		    sprintf (EigenvalueOutputFile, "%s_v_%f_%s_gx_%f_gy_%f.dat",FilePrefix, Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		  else
		    sprintf (EigenvalueOutputFile, "%s_v_%f_%s_spingx_%f_spingy_%f.dat",FilePrefix, Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		}
		else
		{
		  if (Manager.GetBoolean("spin-flux") == false)
		    sprintf (EigenvalueOutputFile, "%s_v_%f_%s_sz_%d_gx_%f_gy_%f.dat",FilePrefix, Manager.GetDouble("v-potential"), FileParameterString, SzValue, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		  else
		    sprintf (EigenvalueOutputFile, "%s_v_%f_%s_sz_%d_spingx_%f_spingy_%f.dat",FilePrefix, Manager.GetDouble("v-potential"), FileParameterString, SzValue, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		}
	      }
	      else
	      {
		if (Manager.GetBoolean("fixed-sz") == false)
		{
		  if (Manager.GetBoolean("spin-flux") == false)
		    sprintf (EigenvalueOutputFile, "%s_v_%f_%s_gx_%f_gy_%f_mus_%f.dat",FilePrefix, Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
		  else
		    sprintf (EigenvalueOutputFile, "%s_v_%f_%s_spingx_%f_spingy_%f_mus_%f.dat",FilePrefix, Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
		}
		else
		{
		  sprintf (EigenvalueOutputFile, "%s_v_%f_%s_sz_%d_gx_%f_gy_%f_mus_%f.dat",FilePrefix, Manager.GetDouble("v-potential"), FileParameterString, SzValue, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
		}
	      }
	    }
	  else
	    {
	      if (Manager.GetDouble("mu-s") == 0.0)
	      {
		if (Manager.GetBoolean("spin-flux") == false)
		  sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s_gx_%f_gy_%f.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		else
		  sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s_spingx_%f_spingy_%f.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	      }
	      else
	      {
		if (Manager.GetBoolean("spin-flux") == false)
		  sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s_gx_%f_gy_%f_mus_%f.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
		else
		  sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s_spingx_%f_spingy_%f_mus_%f.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	      }
	    }
	}
    }
   
   
  Abstract2DTightBindingModel* TightBindingModelUp;
  Abstract2DTightBindingModel* TightBindingModelDown;

  if (Manager.GetString("import-onebodyup") == 0)
    {
      TightBindingModelUp = new TightBindingModelCheckerboardLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true, !(Manager.GetBoolean("single-band")));
      char* BandStructureOutputFile = new char [1024];
      sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
      TightBindingModelUp->WriteBandStructure(BandStructureOutputFile);
    }
  else
    {
      TightBindingModelUp = new Generic2DTightBindingModel(Manager.GetString("import-onebodyup")); 
    }
    
   
   int Sign = 1;
   if (Manager.GetBoolean("spin-flux") == true)
    Sign = -1; 
   if (Manager.GetString("import-onebodydown") == 0)
    {
      TightBindingModelDown = new TightBindingModelCheckerboardLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"), Sign * Manager.GetDouble("gamma-x"), Sign * Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true, !(Manager.GetBoolean("single-band")));
      char* BandStructureOutputFile = new char [1024];
      sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
      TightBindingModelDown->WriteBandStructure(BandStructureOutputFile);
    }
  else
    {
      TightBindingModelDown = new Generic2DTightBindingModel(Manager.GetString("import-onebodydown")); 
    }
    
  if (Manager.GetBoolean("atomic-limit"))
  {
    TightBindingModelUp = new TightBindingModel2DAtomicLimitLattice(NbrSitesX, NbrSitesY, 2, 1, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());  
    TightBindingModelDown = new TightBindingModel2DAtomicLimitLattice(NbrSitesX, NbrSitesY, 2, 1, Sign * Manager.GetDouble("gamma-x"), Sign * Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());  
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
  if (Manager.GetBoolean("no-translation") == true)
    {  
      MaxKx = 0;
      MaxKy = 0;
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
	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      int* NbrInteractingOrbitalsupup;
	      int** InteractingOrbitalsOrbitalIndicesupup;
	      int** InteractingOrbitalsSpatialIndicesupup;
	      double** InteractingOrbitalsPotentialsupup;
	      
	      int* NbrInteractingOrbitalsdowndown;
	      int** InteractingOrbitalsOrbitalIndicesdowndown;
	      int** InteractingOrbitalsSpatialIndicesdowndown;
	      double** InteractingOrbitalsPotentialsdowndown;
	      
	      int* NbrInteractingOrbitalsupdown;
	      int** InteractingOrbitalsOrbitalIndicesupdown;
	      int** InteractingOrbitalsSpatialIndicesupdown;
	      double** InteractingOrbitalsPotentialsupdown;
	      
	      if (Manager.GetBoolean("boson") == false)
		{
		  FCICheckerboardLatticeModelComputeInteractingOrbitals(NbrInteractingOrbitalsupup, InteractingOrbitalsOrbitalIndicesupup, 
								    InteractingOrbitalsSpatialIndicesupup, InteractingOrbitalsPotentialsupup,
								    false, Manager.GetDouble("vupup-potential"),0);
		  FCICheckerboardLatticeModelComputeInteractingOrbitals(NbrInteractingOrbitalsdowndown, InteractingOrbitalsOrbitalIndicesdowndown, 
								    InteractingOrbitalsSpatialIndicesdowndown, InteractingOrbitalsPotentialsdowndown,
								    false, Manager.GetDouble("vdowndown-potential"),0);
		  FCICheckerboardLatticeModelComputeInteractingOrbitals(NbrInteractingOrbitalsupdown, InteractingOrbitalsOrbitalIndicesupdown, 
								    InteractingOrbitalsSpatialIndicesupdown, InteractingOrbitalsPotentialsupdown,
								    true, Manager.GetDouble("u-potential"),Manager.GetDouble("vupdown-potential"));
		}
	      
	      
	      
	      if (Manager.GetBoolean("real-space") == false)
	      {
	  
	      }
	      
	      else
	      {
		RealSymmetricMatrix DensityDensityInteractionupup(TightBindingModelUp->GetNbrBands() * TightBindingModelUp->GetNbrStatePerBand(), true);
		RealSymmetricMatrix DensityDensityInteractiondowndown(TightBindingModelDown->GetNbrBands() * TightBindingModelDown->GetNbrStatePerBand(), true);
		RealSymmetricMatrix DensityDensityInteractionupdown(TightBindingModelUp->GetNbrBands() * TightBindingModelUp->GetNbrStatePerBand(), true);
		 for (int x = 0; x < NbrSitesX; ++x)
		    {
		      for (int y = 0; y < NbrSitesY; ++y)
			{
			  for (int OrbitalIndex = 0; OrbitalIndex < TightBindingModelUp->GetNbrBands(); ++OrbitalIndex)
			    {
			      for (int k = 0; k < NbrInteractingOrbitalsupup[OrbitalIndex]; ++k)
				{
				  DensityDensityInteractionupup.AddToMatrixElement(TightBindingModelUp->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, OrbitalIndex), 
									     TightBindingModelUp->GetRealSpaceTightBindingLinearizedIndexSafe(x + InteractingOrbitalsSpatialIndicesupup[OrbitalIndex][2 * k], 
																	    y + InteractingOrbitalsSpatialIndicesupup[OrbitalIndex][(2 * k) + 1], 
																	    InteractingOrbitalsOrbitalIndicesupup[OrbitalIndex][k]), 
									     InteractingOrbitalsPotentialsupup[OrbitalIndex][k]);
				  
				}
			    }
			    
			    for (int OrbitalIndex = 0; OrbitalIndex < TightBindingModelDown->GetNbrBands(); ++OrbitalIndex)
			    {
			      for (int k = 0; k < NbrInteractingOrbitalsdowndown[OrbitalIndex]; ++k)
				{
				  DensityDensityInteractiondowndown.AddToMatrixElement(TightBindingModelDown->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, OrbitalIndex), 
									     TightBindingModelDown->GetRealSpaceTightBindingLinearizedIndexSafe(x + InteractingOrbitalsSpatialIndicesdowndown[OrbitalIndex][2 * k], 
																	    y + InteractingOrbitalsSpatialIndicesdowndown[OrbitalIndex][(2 * k) + 1], 
																	    InteractingOrbitalsOrbitalIndicesdowndown[OrbitalIndex][k]), 
									     InteractingOrbitalsPotentialsdowndown[OrbitalIndex][k]);
				  
				}
			    }
			    
			    for (int OrbitalIndex = 0; OrbitalIndex < TightBindingModelUp->GetNbrBands() ; ++OrbitalIndex)
			    {
			      for (int k = 0; k < NbrInteractingOrbitalsupdown[OrbitalIndex]; ++k)
				{
				  DensityDensityInteractionupdown.AddToMatrixElement(TightBindingModelUp->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, OrbitalIndex), 
									     TightBindingModelUp->GetRealSpaceTightBindingLinearizedIndexSafe(x + InteractingOrbitalsSpatialIndicesupdown[OrbitalIndex][2 * k], y + InteractingOrbitalsSpatialIndicesupdown[OrbitalIndex][(2 * k) + 1], 						    InteractingOrbitalsOrbitalIndicesupdown[OrbitalIndex][k]), 
									     InteractingOrbitalsPotentialsupdown[OrbitalIndex][k]);
				  
				}
			    }
			}
		    }
		    
		    
		    
		if (Manager.GetBoolean("boson") == false)
		{		
		  if (Manager.GetBoolean("no-translation") == true)
		  {
		    if (Manager.GetBoolean("fixed-sz"))
		    {
		      if (Manager.GetBoolean("gutzwiller") == false)
			Space = new FermionOnLatticeWithSpinRealSpace(NbrParticles, SzValue, NbrSites, 10000000);
		      else
			Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, SzValue, NbrSites, 10000000);
		    }
		    else
		    {
		      if (Manager.GetBoolean("gutzwiller") == false)
			Space = new FermionOnLatticeWithSpinRealSpace(NbrParticles, NbrSites);
		      else
			Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
		    }
		  }
		  
		  else
		  {
		    if (Manager.GetBoolean("fixed-sz"))
		    {
		      if (Manager.GetBoolean("gutzwiller") == false)
			Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, i, NbrSitesX, j, NbrSitesY, 10000000);
		      else
			Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, SzValue, NbrSites, i, NbrSitesX, j, NbrSitesY, 10000000);
		    }
		    else
		    {
		      if (Manager.GetBoolean("gutzwiller") == false)
			Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, i, NbrSitesX, j, NbrSitesY); 
		      else
			Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, i, NbrSitesX, j, NbrSitesY);
		    }
		  }
		  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		  HermitianMatrix TightBindingMatrixUp = TightBindingModelUp->GetRealSpaceTightBindingHamiltonian();
		  HermitianMatrix TightBindingMatrixDown = TightBindingModelDown->GetRealSpaceTightBindingHamiltonian();
		  
// 		  HermitianMatrix TightBindingMatrixUp (NbrSites, true);
// 		  HermitianMatrix TightBindingMatrixDown (NbrSites, true);
		  
		  if (Manager.GetBoolean("no-translation") == true)
		    {
		      Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceHamiltonian (Space, NbrParticles, NbrSites, 
									       TightBindingMatrixUp, TightBindingMatrixDown, DensityDensityInteractionupup,DensityDensityInteractiondowndown, DensityDensityInteractionupdown,
									       Architecture.GetArchitecture(), Memory);
		    }
		  else
		    {
		      Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian (Space, NbrParticles, NbrSites, 
											       i, NbrSitesX, j, NbrSitesY,
											       TightBindingMatrixUp, TightBindingMatrixDown, DensityDensityInteractionupup, DensityDensityInteractiondowndown, DensityDensityInteractionupdown,
											       Architecture.GetArchitecture(), Memory);
		    }
		
		}
	      }
	      
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", i, j);
	      char* EigenstateOutputFile = new char [512];
	      if (Manager.GetBoolean("no-translation") == false)
	      {
		if (Manager.GetBoolean("fixed-sz") == false)
		{
		  sprintf (EigenstateOutputFile, "%s_u_%f_vuu_%f_vdd_%f_vud_%f_%s_kx_%d_ky_%d", FilePrefix, 
		       Manager.GetDouble("u-potential"), Manager.GetDouble("vupup-potential"), Manager.GetDouble("vdowndown-potential"), Manager.GetDouble("vupdown-potential"), FileParameterString, i, j);
		}
		else
		 {
		  sprintf (EigenstateOutputFile, "%s_u_%f_vuu_%f_vdd_%f_vud_%f_%s_kx_%d_ky_%d_sz_%d", FilePrefix, 
		       Manager.GetDouble("u-potential"), Manager.GetDouble("vupup-potential"), Manager.GetDouble("vdowndown-potential"), Manager.GetDouble("vupdown-potential"), FileParameterString, i, j, SzValue);
		}
	      }
	      else
	      {
		if (Manager.GetBoolean("fixed-sz") == false)
		{
		  sprintf (EigenstateOutputFile, "%s_u_%f_vuu_%f_vdd_%f_vud_%f_%s", FilePrefix, 
		       Manager.GetDouble("u-potential"), Manager.GetDouble("vupup-potential"), Manager.GetDouble("vdowndown-potential"), Manager.GetDouble("vupdown-potential"), FileParameterString);
		}
		else
		 {
		  sprintf (EigenstateOutputFile, "%s_u_%f_vuu_%f_vdd_%f_vud_%f_%s_sz_%d", FilePrefix, 
		       Manager.GetDouble("u-potential"), Manager.GetDouble("vupup-potential"), Manager.GetDouble("vdowndown-potential"), Manager.GetDouble("vupdown-potential"), FileParameterString, SzValue);
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
    }
  
 return 0;   
}

// compute the description of the density-density interaction for the unit cell at the origin
//
// nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsPotentials = intensity of each density-density term 
// bosonFlag = true if we are dealing with bosons
// uPotential = nearest neighbor (for fermions) or on-site (for bosons) interaction amplitude
// vPotential = next nearest neighbor (for fermions) or nearest neighbor (for bosons) interaction amplitude

void FCICheckerboardLatticeModelComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
							   int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
							   bool bosonFlag, double uPotential, double vPotential)
{
  nbrInteractingOrbitals = new int[2];
  interactingOrbitalsOrbitalIndices = new int*[2];
  interactingOrbitalsSpatialIndices = new int*[2];
  interactingOrbitalsPotentials = new double*[2];
  if (bosonFlag == false)
    {
      nbrInteractingOrbitals[0] = 1;
      nbrInteractingOrbitals[1] = 3;
      if (vPotential != 0.0)
	{
	  nbrInteractingOrbitals[0] += 2;
	  nbrInteractingOrbitals[1] += 2;
	}
      interactingOrbitalsOrbitalIndices[0] = new int[nbrInteractingOrbitals[0]];
      interactingOrbitalsSpatialIndices[0] = new int[nbrInteractingOrbitals[0] * 2];
      interactingOrbitalsPotentials[0] = new double[nbrInteractingOrbitals[0]];
      interactingOrbitalsOrbitalIndices[1] = new int[nbrInteractingOrbitals[1]];
      interactingOrbitalsSpatialIndices[1] = new int[nbrInteractingOrbitals[1] * 2];
      interactingOrbitalsPotentials[1] = new double[nbrInteractingOrbitals[1]];

      int Index = 0;
      interactingOrbitalsOrbitalIndices[0][Index] = 1;
      interactingOrbitalsSpatialIndices[0][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[0][Index] = uPotential;
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[0][Index] = 0;
	  interactingOrbitalsSpatialIndices[0][2 * Index] = 1;
	  interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 0;
	  interactingOrbitalsPotentials[0][Index] = vPotential;
	  ++Index;	  
	  interactingOrbitalsOrbitalIndices[0][Index] = 0;
	  interactingOrbitalsSpatialIndices[0][2 * Index] = 0;
	  interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 1;
	  interactingOrbitalsPotentials[0][Index] = vPotential;
	  ++Index;	  
	}
      Index = 0;
      interactingOrbitalsOrbitalIndices[1][Index] = 0;
      interactingOrbitalsSpatialIndices[1][2 * Index] = 1;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;
      interactingOrbitalsOrbitalIndices[1][Index] = 0;
      interactingOrbitalsSpatialIndices[1][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 1;
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;
      interactingOrbitalsOrbitalIndices[1][Index] = 0;
      interactingOrbitalsSpatialIndices[1][2 * Index] = 1;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 1;
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[1][Index] = 1;
	  interactingOrbitalsSpatialIndices[1][2 * Index] = 1;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 0;
	  interactingOrbitalsPotentials[1][Index] = vPotential;
	  ++Index;	  
	  interactingOrbitalsOrbitalIndices[1][Index] = 1;
	  interactingOrbitalsSpatialIndices[1][2 * Index] = 0;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 1;
	  interactingOrbitalsPotentials[1][Index] = vPotential;
	  ++Index;	  
	}
    }
  else
    {
      nbrInteractingOrbitals[0] = 1;
      nbrInteractingOrbitals[1] = 1;
      if (vPotential != 0.0)
	{
	  nbrInteractingOrbitals[0] += 1;
	  nbrInteractingOrbitals[1] += 3;
	}
      interactingOrbitalsOrbitalIndices[0] = new int[nbrInteractingOrbitals[0]];
      interactingOrbitalsSpatialIndices[0] = new int[nbrInteractingOrbitals[0] * 2];
      interactingOrbitalsPotentials[0] = new double[nbrInteractingOrbitals[0]];
      interactingOrbitalsOrbitalIndices[1] = new int[nbrInteractingOrbitals[1]];
      interactingOrbitalsSpatialIndices[1] = new int[nbrInteractingOrbitals[1] * 2];
      interactingOrbitalsPotentials[1] = new double[nbrInteractingOrbitals[1]];
  
      int Index = 0;
      interactingOrbitalsOrbitalIndices[0][Index] = 0;
      interactingOrbitalsSpatialIndices[0][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[0][Index] =  uPotential;
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[0][Index] = 1;
	  interactingOrbitalsSpatialIndices[0][2 * Index] = 0;
	  interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 0;
	  interactingOrbitalsPotentials[0][Index] = vPotential;	  
	}
      Index = 0;
      interactingOrbitalsOrbitalIndices[1][Index] = 1;
      interactingOrbitalsSpatialIndices[1][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[1][Index] =  uPotential;
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
	  interactingOrbitalsSpatialIndices[1][2 * Index] = 1;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 0;
	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
	  ++Index;
	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
	  interactingOrbitalsSpatialIndices[1][2 * Index] = 0;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 1;
	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
	  ++Index;
	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
	  interactingOrbitalsSpatialIndices[1][2 * Index] = 1;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 1;
	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
	  ++Index;
	}
    }
}
