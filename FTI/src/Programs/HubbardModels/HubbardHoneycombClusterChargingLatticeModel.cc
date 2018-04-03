#include "Options/Options.h"

#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion.h"
#include "HilbertSpace/FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation.h"
#include "HilbertSpace/FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian.h"
#include "Tools/FTITightBinding/TightBindingModelHaldaneHoneycombLattice.h"



#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"


#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
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


// find linearized index associated with lattice coordinates
//
// siteAlpha = A or B site of honeycomb lattice
// siteX = coordinate along X axis
// siteY = coordinate along Y axis
// return value = linearized index
//
int FindSiteIndex(int siteAlpha, int siteX, int siteY, int NbrSitesX, int NbrSitesY);

int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("HubbardHoneycombClusterChargingLatticeModel" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 2);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "cluster-exclusion", "restrict Hilbert space to plaquettes with at most 6 particles");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-sz", "only evalute a given spin sector (negative if all sz sectors have to be computed)", -1); 
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use the Sz <-> -Sz symmetry");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-parity", "select the  Sz <-> -Sz parity (can be 1 or -1, 0 if both sectors have to be computed", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site (Hubbard) potential strength", 0.0)
;
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hopping amplitude", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);

  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);

  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1); 
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "do not use the momentum as a good quantum number");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption ('\n', "import-onebody", "import information on the tight binding model from a file");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption ('\n', "testhermitian-error", "error threshold when testing hermiticy (0 for machine accuracy)", 0.0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardHaldaneLatticeModel -h" << endl;
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
  int TotalSz = Manager.GetInteger("only-sz");
  int NbrSites = 2*NbrSitesX * NbrSitesY;
  bool ClusterExclusionFlag = Manager.GetBoolean("cluster-exclusion");
  bool SzSymmetryFlag = Manager.GetBoolean("szsymmetrized-basis");
  bool NoTranslationFlag = Manager.GetBoolean("no-translation");  
  
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [64];
  if (Manager.GetBoolean("boson") == false)
    {
      if (SzSymmetryFlag == false)
	{
	  if (ClusterExclusionFlag == false)
	    sprintf (StatisticPrefix, "fermions_hubbard_honeycomb_clustercharging");
	  else
	    sprintf (StatisticPrefix, "fermions_hubbard_honeycomb_clustercharging_exclusion");
	}
      else
	{
	  if (ClusterExclusionFlag == false)
	    sprintf (StatisticPrefix, "fermions_hubbard_honeycomb_clustercharging_szsym");
	  else
	    sprintf (StatisticPrefix, "fermions_hubbard_honeycomb_clustercharging_exclusion_szsym");
	}
    }
  else
    {
      if (SzSymmetryFlag == false)
	{
	  if (ClusterExclusionFlag == false)
	    sprintf (StatisticPrefix, "bosons_hubbard_honeycomb_clustercharging");
	  else
	    sprintf (StatisticPrefix, "bosons_hubbard_honeycomb_clustercharging_exclusion");
	}
      else
	{
	  if (ClusterExclusionFlag == false)
	    sprintf (StatisticPrefix, "bosons_hubbard_honeycomb_clustercharging_szsym");
	  else
	    sprintf (StatisticPrefix, "bosons_hubbard_honeycomb_clustercharging_exclusion_szsym");
	}
    }
    
  
  char* FilePrefix = new char [256];
  sprintf (FilePrefix, "%s_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
  
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t_%f", Manager.GetDouble("t1"));
  
  char* CommentLine = new char [256];
  if (NoTranslationFlag)
    {
      if (SzSymmetryFlag == false)
	{
	  sprintf (CommentLine, "sz");
	}
      else
	{
	  sprintf (CommentLine, "sz szp");
	}
    }
  else
    {
      if (SzSymmetryFlag == false)
	{
	  sprintf (CommentLine, "kx ky sz");
	}
      else
	{
	  sprintf (CommentLine, "kx ky sz szp");
	}
    }

  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetDouble("u-potential") == 0.0)
    sprintf(EigenvalueOutputFile, "%s_%s_sz_%d.dat", FilePrefix, FileParameterString, TotalSz);
  else
    sprintf(EigenvalueOutputFile, "%s_%s_u_%f_sz_%d.dat", FilePrefix, FileParameterString, Manager.GetDouble("u-potential"), TotalSz);

  Abstract2DTightBindingModel* TightBindingModel;
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	ExportOneBody = true;
      TightBindingModel = new TightBindingModelHaldaneHoneycombLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), 
									    0.0, Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);
      
      TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);
      if (ExportOneBody == true)
	{
	  char* BandStructureOutputFile = new char [512];
	  if (Manager.GetString("export-onebodyname") != 0)
	    strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	  else
	    sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
	  if (Manager.GetBoolean("export-onebody") == true)
	    {
	      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
	    }
	  else
	    {
	      TightBindingModel->WriteBandStructureASCII(BandStructureOutputFile);
	    }
	  delete[] BandStructureOutputFile;
	}	  
      return 0;
    }

  if (Manager.GetString("import-onebody") == 0)
    {
      TightBindingModel = new TightBindingModelHaldaneHoneycombLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), 
									    0.0, Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true);
      
      char* BandStructureOutputFile = new char [1024];
      sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
    }
  else
    {
      TightBindingModel = new Generic2DTightBindingModel(Manager.GetString("import-onebody")); 
    }

  RealSymmetricMatrix DensityDensityInteractionupup(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractiondowndown(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractionupdown(NbrSites, true);
  
  int TmpIndex1;
  int TmpIndex2;
  if (Manager.GetDouble("u-potential") != 0.0)
    {
      double UPotential = Manager.GetDouble("u-potential");
      for (int i = 0; i < NbrSites; ++i)
	{
	  DensityDensityInteractionupdown.SetMatrixElement(i, i, UPotential / 9.0);
	}
      for (int siteX = 0; siteX < NbrSitesX; ++siteX)
      {
	for (int siteY = 0; siteY < NbrSitesY; ++siteY)
	{
	  TmpIndex1 = FindSiteIndex(0, siteX, siteY, NbrSitesX, NbrSitesY);
	  TmpIndex2 = FindSiteIndex(1, siteX, siteY, NbrSitesX, NbrSitesY);
// 	  NN interaction
	  DensityDensityInteractionupup.SetMatrixElement(TmpIndex1, TmpIndex2, 2.0*UPotential/9.0);
	  DensityDensityInteractionupup.SetMatrixElement(TmpIndex1, FindSiteIndex(1, siteX + 1, siteY - 1, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  DensityDensityInteractionupup.SetMatrixElement(FindSiteIndex(0, siteX, siteY + 1, NbrSitesX, NbrSitesY), TmpIndex2, 2.0*UPotential/9.0);
	  
	  DensityDensityInteractiondowndown.SetMatrixElement(TmpIndex1, TmpIndex2, 2.0*UPotential/9.0);
	  DensityDensityInteractiondowndown.SetMatrixElement(TmpIndex1, FindSiteIndex(1, siteX + 1, siteY - 1, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  DensityDensityInteractiondowndown.SetMatrixElement(FindSiteIndex(0, siteX, siteY + 1, NbrSitesX, NbrSitesY), TmpIndex2, 2.0*UPotential/9.0);
	  
	  DensityDensityInteractionupdown.SetMatrixElement(TmpIndex1, TmpIndex2, 2.0*UPotential/9.0);
	  DensityDensityInteractionupdown.SetMatrixElement(TmpIndex1, FindSiteIndex(1, siteX + 1, siteY - 1, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  DensityDensityInteractionupdown.SetMatrixElement(FindSiteIndex(0, siteX, siteY + 1, NbrSitesX, NbrSitesY), TmpIndex2, 2.0*UPotential/9.0);
	  
// 	  NNN interaction
	  for (int alpha = 0; alpha < 2; ++alpha)
	  {
	    DensityDensityInteractionupup.SetMatrixElement(TmpIndex1, FindSiteIndex(alpha, siteX, siteY + 1, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	    DensityDensityInteractionupup.SetMatrixElement(TmpIndex1, FindSiteIndex(alpha, siteX + 1, siteY, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	    DensityDensityInteractionupup.SetMatrixElement(TmpIndex1, FindSiteIndex(alpha, siteX + 1, siteY - 1, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	    
	    DensityDensityInteractiondowndown.SetMatrixElement(TmpIndex1, FindSiteIndex(alpha, siteX, siteY + 1, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	    DensityDensityInteractiondowndown.SetMatrixElement(TmpIndex1, FindSiteIndex(alpha, siteX + 1, siteY, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	    DensityDensityInteractiondowndown.SetMatrixElement(TmpIndex1, FindSiteIndex(alpha, siteX + 1, siteY - 1, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	    
	    DensityDensityInteractionupdown.SetMatrixElement(TmpIndex1, FindSiteIndex(alpha, siteX, siteY + 1, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	    DensityDensityInteractionupdown.SetMatrixElement(TmpIndex1, FindSiteIndex(alpha, siteX + 1, siteY, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	    DensityDensityInteractionupdown.SetMatrixElement(TmpIndex1, FindSiteIndex(alpha, siteX + 1, siteY - 1, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  }
	  
// 	  NNNN interaction
	  DensityDensityInteractionupup.SetMatrixElement(TmpIndex1, FindSiteIndex(1, siteX + 1, siteY, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  DensityDensityInteractionupup.SetMatrixElement(TmpIndex1, FindSiteIndex(1, siteX + 1, siteY - 2, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  DensityDensityInteractionupup.SetMatrixElement(TmpIndex2, FindSiteIndex(0, siteX + 1, siteY, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  
	  DensityDensityInteractiondowndown.SetMatrixElement(TmpIndex1, FindSiteIndex(1, siteX + 1, siteY, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  DensityDensityInteractiondowndown.SetMatrixElement(TmpIndex1, FindSiteIndex(1, siteX + 1, siteY - 2, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  DensityDensityInteractiondowndown.SetMatrixElement(TmpIndex2, FindSiteIndex(0, siteX + 1, siteY, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  
	  DensityDensityInteractionupdown.SetMatrixElement(TmpIndex1, FindSiteIndex(1, siteX + 1, siteY, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  DensityDensityInteractionupdown.SetMatrixElement(TmpIndex1, FindSiteIndex(1, siteX + 1, siteY - 2, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  DensityDensityInteractionupdown.SetMatrixElement(TmpIndex2, FindSiteIndex(0, siteX + 1, siteY, NbrSitesX, NbrSitesY), 2.0*UPotential/9.0);
	  
	}
      }
    }

  bool FirstRunFlag = true;
  ParticleOnSphereWithSpin* Space = 0;
  AbstractHamiltonian* Hamiltonian = 0;

  if (NoTranslationFlag)
    {
      if (SzSymmetryFlag)
      {
	cout << "Spin flip symmetry not implemented in the absence of 2D translations" << endl;
	return 0;
      }      
      
      if (ClusterExclusionFlag == false)
	Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, TotalSz, NbrSites);
      else
	Space = new FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion (NbrParticles, TotalSz, NbrSitesX, NbrSitesY);
      
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  
      HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
      cout << TightBindingMatrix << endl;
      
      Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceHamiltonian(Space, NbrParticles, NbrSites,
									  TightBindingMatrix, TightBindingMatrix,
									  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
									  DensityDensityInteractionupdown, 
									  Architecture.GetArchitecture(), Memory);
	  
      char* ContentPrefix = new char[256];
      sprintf (ContentPrefix, "%d", TotalSz);
      char* EigenstateOutputFile;
      char* TmpExtention = new char [512];
      sprintf (TmpExtention, "_sz_%d", TotalSz);
      char* TmpExtentionSpectrum = new char [32];
      sprintf (TmpExtentionSpectrum, "_sz_%d.dat", TotalSz);
      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, TmpExtentionSpectrum, TmpExtention);
      
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
      int MinXMomentum = 0;
      int MaxXMomentum = NbrSitesX - 1;
      if (Manager.GetInteger("only-kx") >= 0)
	{
	  MaxXMomentum = Manager.GetInteger("only-kx");
	  MinXMomentum = MaxXMomentum;
	}
      for (int XMomentum = MinXMomentum; XMomentum <= MaxXMomentum; ++XMomentum)
	{
	  int MinYMomentum = 0;
	  int MaxYMomentum = NbrSitesY - 1;
	  if (Manager.GetInteger("only-ky") >= 0)
	    {
	      MaxYMomentum = Manager.GetInteger("only-ky");
	      MinYMomentum = MaxYMomentum;
	    }
	  for (int YMomentum = MinYMomentum; YMomentum <= MaxYMomentum; ++YMomentum)
	    {
	      int SzParitySector = -1;
	      int MaxSzParitySector = 1;
	      if (SzSymmetryFlag == false)
		{
		  SzParitySector = 1;
		}
	      else
		{
		  if ((Manager.GetInteger("sz-parity") != 0) && (TotalSz == 0))
		    {
		      SzParitySector = Manager.GetInteger("sz-parity");
		      MaxSzParitySector = SzParitySector;
		    }
		}
	      for (; SzParitySector <= MaxSzParitySector; SzParitySector += 2)
		{
// 		  ParticleOnSphereWithSpin* Space = 0;
// 		  AbstractHamiltonian* Hamiltonian = 0;
		  if (SzSymmetryFlag == false)
		    {
		      cout << "Kx = " << XMomentum << "  Ky = " << YMomentum << "  Sz = " << TotalSz <<  endl;
		      if (ClusterExclusionFlag == false)
			Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, TotalSz, NbrSites, XMomentum, NbrSitesX,YMomentum, NbrSitesY);
		      else
			Space = new FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation (NbrParticles, TotalSz, NbrSites, XMomentum, NbrSitesX, YMomentum, NbrSitesY);
		    }
		  else
		    {
		      bool MinusParitySector = true;
		      if (SzParitySector == 1)
			MinusParitySector = false;
		      cout << "Kx = " << XMomentum << "  Ky = " << YMomentum << "  Sz = " << TotalSz << "  SzParity = " << SzParitySector<< endl;
		      if (ClusterExclusionFlag == false)
			Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, TotalSz, NbrSites, MinusParitySector, 
												 XMomentum, NbrSitesX,
												 YMomentum, NbrSitesY);
		      else
			Space = new FermionOnHoneycombLatticeWithSpinSzSymmetryRealSpacePlaquetteExclusionAnd2DTranslation (NbrParticles, TotalSz, NbrSites, MinusParitySector, XMomentum, NbrSitesX, YMomentum, NbrSitesY);
		    }
		  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		    Memory = Architecture.GetArchitecture()->GetLocalMemory();
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		  
		  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
		  Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian(Space, NbrParticles, NbrSites,XMomentum, NbrSitesX, YMomentum, NbrSitesY, TightBindingMatrix, TightBindingMatrix, DensityDensityInteractionupup,  DensityDensityInteractiondowndown,DensityDensityInteractionupdown, Architecture.GetArchitecture(), Memory);
		  
		  char* ContentPrefix = new char[256];
		  if (SzSymmetryFlag == false)
		    {
		      sprintf (ContentPrefix, "%d %d %d", XMomentum, YMomentum, TotalSz);
		    }
		  else
		    {
		      sprintf (ContentPrefix, "%d %d %d %d", XMomentum, YMomentum, TotalSz, SzParitySector);
		    }
		  char* EigenstateOutputFile;
		  char* TmpExtention = new char [512];
		  if (SzSymmetryFlag == false)
		    {
		      sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d", XMomentum, YMomentum, TotalSz);
		    }
		  else
		    {
		      sprintf (TmpExtention, "_szp_%d_kx_%d_ky_%d_sz_%d", SzParitySector, XMomentum, YMomentum, TotalSz);
		    }
		  char* TmpExtentionSpectrum = new char [32];
		  sprintf (TmpExtentionSpectrum, "_sz_%d.dat", TotalSz);
		  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, TmpExtentionSpectrum, TmpExtention);
		  
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
    }  
  return 0;
}


// find linearized index associated with lattice coordinates
//
// siteAlpha = A or B site of honeycomb lattice
// siteX = coordinate along X axis
// siteY = coordinate along Y axis
// return value = linearized index
//
int FindSiteIndex(int siteAlpha, int siteX, int siteY, int NbrSitesX, int NbrSitesY)
{
  if (siteX < 0)
    siteX += NbrSitesX;
  if (siteX >= NbrSitesX)
    siteX -= NbrSitesX;
  
  if (siteY < 0)
    siteY += NbrSitesY;
  if (siteY >= NbrSitesY)
    siteY -= NbrSitesY;
  
  int TmpIndex = siteAlpha + 2 * (NbrSitesY * siteX + siteY);
  
  return TmpIndex;
}