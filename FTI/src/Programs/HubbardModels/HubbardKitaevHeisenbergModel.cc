#include "Options/Options.h"

#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergAnd2DTranslationHamiltonian.h"

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


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("HubbardKitaevHeisenbergModel" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sites", "total number of sites (if negative, guess it from the geometry file)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the Gutzwiller projection");
  (*SystemGroup) += new BooleanOption  ('\n', "stripe", "model geometry is a stripe");
  (*SystemGroup) += new BooleanOption  ('\n', "torus", "model geometry is a torus");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrsites-x", "number of sites in the x direction for the torus geometry", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrsites-y", "number of sites in the y direction for the torus geometry", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use the Sz <-> -Sz symmetry");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-parity", "select the  Sz <-> -Sz parity (can be 1 or -1, 0 if both sectors have to be computed", 0);

  (*SystemGroup) += new SingleStringOption  ('\n', "geometry-file", "name of an optional file that gives the position of the different bonds and their nature");
  (*SystemGroup) += new SingleStringOption  ('\n', "geometry-name", "name of the geometry used", "unknown");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site (Hubbard) potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "isotropic-t", "isotropic spin nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "anisotropic-t", "anisotropic nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "j1", "strength of the neareast neighbor Heisenberg interaction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "j2", "strength of the neareast neighbor anisotropic interaction", 1.0);

  (*SystemGroup) += new BooleanOption  ('\n', "xperiodic-boundary", "use periodic boundary conditions in the x direction");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "x-momentum", "set the momentum along the x direction (negative if all momentum sectors have to be evaluated)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "x-periodicity", "periodicity in the number of site index that implements the periodic boundary condition in the x direction", 4);
  (*SystemGroup) += new BooleanOption  ('\n', "2dperiodic-boundaries", "use periodic boundary conditions in the x and y directions");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "y-momentum", "set the momentum along the y direction (negative if all momentum sectors have to be evaluated)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "y-periodicity", "periodicity in the number of site index that implements the periodic boundary condition in the y direction", 2);

  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by .#.vec");
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
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardKitaevHeisenbergModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSites = Manager.GetInteger("nbr-sites"); 
  bool GutzwillerFlag = Manager.GetBoolean("gutzwiller");
  bool SzSymmetryFlag = Manager.GetBoolean("szsymmetrized-basis");
  bool StripeFlag = Manager.GetBoolean("stripe"); 
  bool TorusFlag = Manager.GetBoolean("torus");
  int NbrSitesX = Manager.GetInteger("nbrsites-x"); 
  int NbrSitesY = Manager.GetInteger("nbrsites-y"); 

 
  if ((StripeFlag == false) && (TorusFlag == false) && (Manager.GetString("geometry-file") == 0))
    {
      cout << "Error. A lattice geometry has to be specified" << endl; 
      return -1;
    }

  if ((Manager.GetString("geometry-file") != 0) && (NbrSites < 0))
    {
      MultiColumnASCIIFile GeometryFile;
      if (GeometryFile.Parse(Manager.GetString("geometry-file")) == false)
	{
	  GeometryFile.DumpErrors(cout);
	  return -1;
	}
      if (GeometryFile.GetNbrColumns() < 3)
	{
	  cout << "Error. " << GeometryFile.GetNbrColumns() << " has a wrong number of columns" << endl;
	  return -1;
	}
      int* TmpColumn = GeometryFile.GetAsIntegerArray(0);
      int LargestSiteIndex = -1; 
      for (int i = 0; i < GeometryFile.GetNbrLines(); ++i)
	{
	  if (TmpColumn[i] > LargestSiteIndex)
	    LargestSiteIndex = TmpColumn[i];
	}
      TmpColumn = GeometryFile.GetAsIntegerArray(1);
      for (int i = 0; i < GeometryFile.GetNbrLines(); ++i)
	{
	  if (TmpColumn[i] > LargestSiteIndex)
	    LargestSiteIndex = TmpColumn[i];
	}
      NbrSites = LargestSiteIndex + 1;
    }

  if ((Manager.GetBoolean("xperiodic-boundary") == true)  && 
      (((TorusFlag == false) && ((NbrSites % Manager.GetInteger("x-periodicity")) != 0)) ||
       ((TorusFlag == true) && ((NbrSitesX % Manager.GetInteger("x-periodicity")) != 0))))
    {
      cout << "Error. The number of sites is not compatible with the periodicity in the x direction" << endl; 
      return -1;
    }
    
//   if ((StripeFlag) && ((NbrSites % 4) != 2))
//   {
//    cout << "Error: number of sites should be of the form 4n + 2 for stripe geometry " << endl; 
//    return -1;
//   }
  

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [64];
  if (Manager.GetBoolean("boson") == false)
    {
      if (Manager.GetBoolean("xperiodic-boundary") == false)
	{
	  if (Manager.GetBoolean("2dperiodic-boundaries") == false)
	    {
	      if (SzSymmetryFlag == false)
		{
		  if (GutzwillerFlag == false)
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg");
		  else
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller");
		}
	      else
		{
		  if (GutzwillerFlag == false)
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_szsym");
		  else
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_szsym");
		}
	    }
	  else
	    {
	      if (GutzwillerFlag == false)
		sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_xymomentum_%d_%d", (int) Manager.GetInteger("x-periodicity"), 
			 (int) Manager.GetInteger("y-periodicity"));
	      else
		sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_xymomentum_%d_%d", (int) Manager.GetInteger("x-periodicity"),
			 (int) Manager.GetInteger("y-periodicity"));
	    }
	}
      else
	{
	  if (SzSymmetryFlag == false)
	    {
	      if (GutzwillerFlag == false)
		sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	      else
		sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	    }
	  else
	    {
	      if (GutzwillerFlag == false)
		sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_szsym_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	      else
		sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_szsym_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("xperiodic-boundary") == false)
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "bosons_kitaev_heisenberg");
	  else
	    sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_gutzwiller");
	}
      else
	{
	  if (SzSymmetryFlag == false)
	    {
	      if (GutzwillerFlag == false)
		sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	      else
		sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_gutzwiller_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	    }
	  else
	    {
	      if (GutzwillerFlag == false)
		sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_szsym_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	      else
		sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_gutzwiller_szsym_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	    }
	}
    }
    
  

  char* FilePrefix = new char [256];
  if (StripeFlag == true)
    {
      sprintf (FilePrefix, "%s_stripe_n_%d_x_%d", StatisticPrefix, NbrParticles, NbrSites);
    }
  else
    {
      if (TorusFlag == true)
	{
	  sprintf (FilePrefix, "%s_torus_nx_%d_ny_%d_n_%d_x_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
	}
      else
	{
	  sprintf (FilePrefix, "%s_%s_n_%d_x_%d", StatisticPrefix, Manager.GetString("geometry-name"), NbrParticles, NbrSites);
	}
    }
  
  
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t_%g_tK_%g_j1_%g_j2_%g", Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), Manager.GetDouble("j1"), Manager.GetDouble("j2"));

  char* CommentLine = new char [256];
  if (Manager.GetBoolean("xperiodic-boundary") == false)
    {
      if (Manager.GetBoolean("2dperiodic-boundaries") == false)
	{
	  sprintf (CommentLine, "");
	}
      else
	{
	  sprintf (CommentLine, "kx ky");
	}
    }
  else
    {
      if (SzSymmetryFlag == false)
	{
	  sprintf (CommentLine, "kx");
	}
      else
	{
	  sprintf (CommentLine, "kx szp");
	}
    }
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    {
      if (Manager.GetDouble("u-potential") == 0.0)
	sprintf(EigenvalueOutputFile, "%s_%s.dat", FilePrefix, FileParameterString);
      else
	sprintf(EigenvalueOutputFile, "%s_%s_u_%f.dat", FilePrefix, FileParameterString, Manager.GetDouble("u-potential"));
    }

  bool FirstRunFlag = true;
  if ((Manager.GetBoolean("xperiodic-boundary") == false) && (Manager.GetBoolean("2dperiodic-boundaries") == false))
    {
      int SzParitySector = -1;
      int MaxSzParitySector = 1;
      if (SzSymmetryFlag == false)
	{
	  SzParitySector = 1;
	}
      else
	{
	  if (Manager.GetInteger("sz-parity") != 0)
	    {
	      SzParitySector = Manager.GetInteger("sz-parity");
	      MaxSzParitySector = SzParitySector;
	    }
	}
      for (; SzParitySector <= MaxSzParitySector; SzParitySector += 2)
	{
	  ParticleOnSphereWithSpin* Space = 0;
	  AbstractHamiltonian* Hamiltonian = 0;
	  if (SzSymmetryFlag == false)
	    {
	      if (GutzwillerFlag == false)
		Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
	      else
		Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
	    }
	  else
	    {
	      bool MinusParitySector = true;
	      if (SzParitySector == 1)
		MinusParitySector = false;
	      if (GutzwillerFlag == false)
		Space = new FermionOnLatticeWithSpinSzSymmetryRealSpace (NbrParticles, NbrSites, MinusParitySector);
	      else
		Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites, MinusParitySector);
	    }
	  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  
	  Hamiltonian = new ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian(Space, NbrParticles, NbrSites, Manager.GetString("geometry-file"), Manager.GetDouble("isotropic-t"), 
										 Manager.GetDouble("anisotropic-t"), Manager.GetDouble("u-potential"), Manager.GetDouble("j1"), 
										 Manager.GetDouble("j2"), Architecture.GetArchitecture(), Memory);
	  
	  char* ContentPrefix = new char[256];
	  if (SzSymmetryFlag == false)
	    {
	      sprintf (ContentPrefix, "");
	    }
	  else
	    {
	      sprintf (ContentPrefix, "%d", SzParitySector);
	    }

	  char* EigenstateOutputFile;
	  if (Manager.GetString("eigenstate-file") != 0)
	    {
	      EigenstateOutputFile = new char [512];
	      sprintf (EigenstateOutputFile, "%s", Manager.GetString("eigenstate-file"));
	    }
	  else
	    {
	      char* TmpExtention = new char [512];
	      if (SzSymmetryFlag == false)
		{
		  sprintf (TmpExtention, "");
		}
	      else
		{
		  sprintf (TmpExtention, "_szp_%d", SzParitySector);
		}
	      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
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
  if (Manager.GetBoolean("xperiodic-boundary") == false)
    {
      int XPeriodicity = Manager.GetInteger("x-periodicity");
      int MinXMomentum = 0;
      int MaxXMomentum = (NbrSites / XPeriodicity) - 1;
      if (Manager.GetInteger("x-momentum") >= 0)
	{
	  MaxXMomentum = Manager.GetInteger("x-momentum");
	  MinXMomentum = MaxXMomentum;
	}
      for (int XMomentum = MinXMomentum; XMomentum <= MaxXMomentum; ++XMomentum)
	{
	  int SzParitySector = -1;
	  int MaxSzParitySector = 1;
	  if (SzSymmetryFlag == false)
	    {
	      SzParitySector = 1;
	    }
	  else
	    {
	      if (Manager.GetInteger("sz-parity") != 0)
		{
		  SzParitySector = Manager.GetInteger("sz-parity");
		  MaxSzParitySector = SzParitySector;
		}
	    }
	  for (; SzParitySector <= MaxSzParitySector; SzParitySector += 2)
	    {
	      ParticleOnSphereWithSpin* Space = 0;
	      AbstractHamiltonian* Hamiltonian = 0;
	      if (SzSymmetryFlag == false)
		{
		  if (GutzwillerFlag == false)
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, Manager.GetInteger("x-periodicity"));
		  else
		    Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, Manager.GetInteger("x-periodicity"));
		}
	      else
		{
		  bool MinusParitySector = true;
		  if (SzParitySector == 1)
		     MinusParitySector = false;
		  if (GutzwillerFlag == false)
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, Manager.GetInteger("x-periodicity"), MinusParitySector);
		  else
		    Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, 
														    Manager.GetInteger("x-periodicity"),MinusParitySector);
		}
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      
	      Hamiltonian = new ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian(Space, NbrParticles, NbrSites, XMomentum, Manager.GetInteger("x-periodicity"), 
												     Manager.GetString("geometry-file"), Manager.GetDouble("isotropic-t"), 
												     Manager.GetDouble("anisotropic-t"), Manager.GetDouble("u-potential"), 
												     Manager.GetDouble("j1"), Manager.GetDouble("j2"), 
												     Architecture.GetArchitecture(), Memory);
	      
	      char* ContentPrefix = new char[256];
	      if (SzSymmetryFlag == false)
		{
		  sprintf (ContentPrefix, "%d", XMomentum);
		}
	      else
		{
		  sprintf (ContentPrefix, "%d %d", XMomentum, SzParitySector);
		}
	      char* EigenstateOutputFile;
	      if (Manager.GetString("eigenstate-file") != 0)
		{
		  EigenstateOutputFile = new char [512];
		  sprintf (EigenstateOutputFile, "%s", Manager.GetString("eigenstate-file"));
		}
	      else
		{
		  char* TmpExtention = new char [512];
		  if (SzSymmetryFlag == false)
		    {
		      sprintf (TmpExtention, "_kx_%d", XMomentum);
		    }
		  else
		    {
		      sprintf (TmpExtention, "_szp_%d_kx_%d", SzParitySector, XMomentum);
		    }
		  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
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
  else
    {
      int XPeriodicity = Manager.GetInteger("x-periodicity");
      int MinXMomentum = 0;
      int MaxXMomentum = (NbrSites / XPeriodicity) - 1;
      if (Manager.GetInteger("x-momentum") >= 0)
	{
	  MaxXMomentum = Manager.GetInteger("x-momentum");
	  MinXMomentum = MaxXMomentum;
	}
      for (int XMomentum = MinXMomentum; XMomentum <= MaxXMomentum; ++XMomentum)
	{
	  int MinYMomentum = 0;
	  int MaxYMomentum = Manager.GetInteger("y-periodicity") - 1;
	  if (Manager.GetInteger("y-momentum") >= 0)
	    {
	      MaxYMomentum = Manager.GetInteger("y-momentum");
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
		  if (Manager.GetInteger("sz-parity") != 0)
		    {
		      SzParitySector = Manager.GetInteger("sz-parity");
		      MaxSzParitySector = SzParitySector;
		    }
		}
	      for (; SzParitySector <= MaxSzParitySector; SzParitySector += 2)
		{
		  ParticleOnSphereWithSpin* Space = 0;
		  AbstractHamiltonian* Hamiltonian = 0;
		  if (SzSymmetryFlag == false)
		    {
		      if (GutzwillerFlag == false)
			Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, Manager.GetInteger("x-periodicity"),
										       YMomentum, Manager.GetInteger("y-momentum"));
		      else
			Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, Manager.GetInteger("x-periodicity"),
													      YMomentum, Manager.GetInteger("y-momentum"));
		    }
		  else
		    {
		      bool MinusParitySector = true;
		      if (SzParitySector == 1)
			MinusParitySector = false;
// 		      if (GutzwillerFlag == false)
// 			Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, Manager.GetInteger("x-periodicity"), MinusParitySector);
// 		      else
// 			Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, 
// 															Manager.GetInteger("x-periodicity"),MinusParitySector);
		      cout << "error, 2d translations and Sz<->-Sz is not an implemented symmetry" << endl;
		      return -1;
		    }
		  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		    Memory = Architecture.GetArchitecture()->GetLocalMemory();
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		  
		  Hamiltonian = new ParticleOnLatticeWithSpinKitaevHeisenbergAnd2DTranslationHamiltonian(Space, NbrParticles, NbrSites, XMomentum, (NbrSites / XPeriodicity),
													 YMomentum, Manager.GetInteger("y-momentum"),
													 Manager.GetString("geometry-file"), Manager.GetDouble("isotropic-t"), 
													 Manager.GetDouble("anisotropic-t"), Manager.GetDouble("u-potential"), 
													 Manager.GetDouble("j1"), Manager.GetDouble("j2"), 
													 Architecture.GetArchitecture(), Memory);
		  
		  char* ContentPrefix = new char[256];
		  if (SzSymmetryFlag == false)
		    {
		      sprintf (ContentPrefix, "%d %d", XMomentum, YMomentum);
		    }
		  else
		    {
		      sprintf (ContentPrefix, "%d %d %d", XMomentum, YMomentum, SzParitySector);
		    }
		  char* EigenstateOutputFile;
		  if (Manager.GetString("eigenstate-file") != 0)
		    {
		      EigenstateOutputFile = new char [512];
		      sprintf (EigenstateOutputFile, "%s", Manager.GetString("eigenstate-file"));
		    }
		  else
		    {
		      char* TmpExtention = new char [512];
		      if (SzSymmetryFlag == false)
			{
			  sprintf (TmpExtention, "_kx_%d_ky%d", XMomentum, YMomentum);
			}
		      else
			{
			  sprintf (TmpExtention, "_szp_%d_kx_%d_ky_%d", SzParitySector, XMomentum, YMomentum);
			}
		      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
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
    }  
  return 0;
}
