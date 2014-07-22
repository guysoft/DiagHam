#include "Options/Options.h"

#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian.h"

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
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sites", "total number of sites (if negartive, guess it from the geometry file)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the Gutzwiller projection");
  (*SystemGroup) += new BooleanOption  ('\n', "stripe", "model geometry is a stripe");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");

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
  bool StripeFlag = Manager.GetBoolean("stripe");
  
  if ((StripeFlag == false) && (Manager.GetString("geometry-file") == 0))
    {
      cout << "Error. A lattice geometry has to be specified" << endl; 
      return -1;
    }

  if ((Manager.GetString("geometry-file") != 0) && (NbrSites< 0))
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

  if ((Manager.GetBoolean("xperiodic-boundary") == true)  && ((NbrSites % Manager.GetInteger("x-periodicity")) != 0))
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
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg");
	  else
	    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller");
	}
      else
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	  else
	    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
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
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	  else
	    sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_gutzwiller_xmomentum_%d", (int) Manager.GetInteger("x-periodicity"));
	}
    }
    
  

  char* FilePrefix = new char [256];
  if (StripeFlag)
    sprintf (FilePrefix, "%s_stripe_n_%d_x_%d", StatisticPrefix, NbrParticles, NbrSites);
  else
    sprintf (FilePrefix, "%s_%s_n_%d_x_%d", StatisticPrefix, Manager.GetString("geometry-name"), NbrParticles, NbrSites);
  
  
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t_%g_tK_%g_j1_%g_j2_%g", Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), Manager.GetDouble("j1"), Manager.GetDouble("j2"));

  char* CommentLine = new char [256];
  if (Manager.GetBoolean("xperiodic-boundary") == false)
    sprintf (CommentLine, "");
  else
    sprintf (CommentLine, "kx");
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
  if (Manager.GetBoolean("xperiodic-boundary") == false)
    {
      ParticleOnSphereWithSpin* Space = 0;
      AbstractHamiltonian* Hamiltonian = 0;
      if (GutzwillerFlag == false)
	Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
      else
	Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
      
      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      
      Hamiltonian = new ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian(Space, NbrParticles, NbrSites, Manager.GetString("geometry-file"), Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), Manager.GetDouble("u-potential"), Manager.GetDouble("j1"), Manager.GetDouble("j2"), Architecture.GetArchitecture(), Memory);
      
      char* ContentPrefix = new char[256];
      sprintf (ContentPrefix, "0");
      char* EigenstateOutputFile;
      if (Manager.GetString("eigenstate-file") != 0)
	{
	  EigenstateOutputFile = new char [512];
	  sprintf (EigenstateOutputFile, "%s", Manager.GetString("eigenstate-file"));
	}
      else
	{
	  char* TmpExtention = new char [512];
	  sprintf (TmpExtention, "");
	  //       sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
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
	  ParticleOnSphereWithSpin* Space = 0;
	  AbstractHamiltonian* Hamiltonian = 0;
	  if (GutzwillerFlag == false)
	    Space = new FermionOnLatticeWithSpinRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, Manager.GetInteger("x-periodicity"));
	  else
	    Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, Manager.GetInteger("x-periodicity"));

	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  
	  Hamiltonian = new ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian(Space, NbrParticles, NbrSites, XMomentum, Manager.GetInteger("x-periodicity"), Manager.GetString("geometry-file"), Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), Manager.GetDouble("u-potential"), Manager.GetDouble("j1"), Manager.GetDouble("j2"), Architecture.GetArchitecture(), Memory);

	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d", XMomentum);
	  char* EigenstateOutputFile;
	  if (Manager.GetString("eigenstate-file") != 0)
	    {
	      EigenstateOutputFile = new char [512];
	      sprintf (EigenstateOutputFile, "%s", Manager.GetString("eigenstate-file"));
	    }
	  else
	    {
	      char* TmpExtention = new char [512];
	      sprintf (TmpExtention, "_kx_%d", XMomentum);
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
  
  return 0;
}
