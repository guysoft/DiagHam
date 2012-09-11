#include "Options/Options.h"

#include "HilbertSpace/FermionOnCubicLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU4SpinMomentumSpace.h"


#include "Hamiltonian/ParticleOnCubicLatticeFourBandPyrochloreHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelPyrochloreLattice.h"
#include "Tools/FTITightBinding/TightBindingModel3DAtomicLimitLattice.h"

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
  cout.precision(14);

  OptionManager Manager ("FQHEQuantumSpinHall3DPyrochlore" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('z', "nbr-sitez", "number of sites along the z direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kz", "only evalute a given y momentum sector (negative if all kz sectors have to be computed)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site potential strength between identical spins for bosons or repulsive nearest neighbor site potential strength between identical spins for fermions", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive on-site potential strength between opposite spins", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive nearest neighbor site potential strength between opposite spins (fermions only)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "lambda-nn", "spin orbit coupling to neareast neighbor sites", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "lambda-nextnn", "spin orbit coupling to next neareast neighbor sites", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-z", "boundary condition twisting angle along z (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new BooleanOption ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new BooleanOption ('\n', "four-bands", "perform the calculations within the full four band model");
  (*SystemGroup) += new BooleanOption ('\n', "project-fourbands", "project the hamiltonian from the four band model to the two band model");
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
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-hamiltonian", "export the hamiltonian in a column formatted ASCII file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEQuantumSpinHall3DPyrochlore -h" << endl;
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
  int NbrSitesZ = Manager.GetInteger("nbr-sitez"); 
  int TotalNbrSites = NbrSitesX * NbrSitesY * NbrSitesZ;
  long Memory = ((long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky kz ");
  char* FilePrefix = new char [512];
  sprintf (FilePrefix, "%s_quantumspinhall3d_pyrochlore_n_%d_x_%d_y_%d_z_%d_sonn_%f_sonnn_%f",  StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("lambda-nn"), Manager.GetDouble("lambda-nextnn"));
  char* EigenvalueOutputFile = new char [1024];
  sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_w_%f_gx_%f_gy_%f_gz_%f.dat", FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
	   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	ExportOneBody = true;
      TightBindingModelPyrochloreLattice TightBindingModel(NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("lambda-nn"), Manager.GetDouble("lambda-nextnn"), 
							   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), ExportOneBody);
      TightBindingModel.WriteAsciiSpectrum(EigenvalueOutputFile);
      double BandSpread = TightBindingModel.ComputeBandSpread(0);
      double DirectBandGap = TightBindingModel.ComputeDirectBandGap(0);
      cout << "Spread = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	{
	  char* BandStructureOutputFile = new char [512];
	  if (Manager.GetString("export-onebodyname") != 0)
	    strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	  else
	    sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
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
  int MinKz = 0;
  int MaxKz = NbrSitesZ - 1;
  if (Manager.GetInteger("only-kz") >= 0)
    {						
      MinKz = Manager.GetInteger("only-kz");
      MaxKz = MinKz;
    }
  TightBindingModelPyrochloreLattice TightBindingModel(NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("lambda-nn"), Manager.GetDouble("lambda-nextnn"), 
  						       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));
  
  TightBindingModelPyrochloreLattice TightBindingModel2(NbrSitesY, NbrSitesX, NbrSitesZ, Manager.GetDouble("lambda-nn"), Manager.GetDouble("lambda-nextnn"), 
  						       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));
  
  
  ComplexMatrix OneBodyBasis1 =  TightBindingModel.GetOneBodyMatrix(1);
  ComplexMatrix OneBodyBasis2 =  TightBindingModel2.GetOneBodyMatrix(1);

  for (int i = 0; i < 8; ++i)
    {
      double Sum = 0.0;
      for (int j = 0; j < 8; ++j)
	{
	  Complex Tmp = 0.0;
	  for (int k = 0; k < 8; ++k)
	    Tmp += Conj(OneBodyBasis2[i][k]) * OneBodyBasis1[j][k];
	  cout << i << " " << j << " : " << SqrNorm(Tmp) << endl;
	  Sum += SqrNorm(Tmp);
	}
      cout << "sum = " << Sum << endl;
    }
//     double* TmpChem = new double[8];
//     TmpChem[0] = -1.0;
//     TmpChem[1] = -1.0;
//     TmpChem[2] = 1.0;
//     TmpChem[3] = 1.0;
//     TmpChem[4] = -1.0;
//     TmpChem[5] = -1.0;
//     TmpChem[6] = 1.0;
//     TmpChem[7] = 1.0;
//     TightBindingModel3DAtomicLimitLattice TightBindingModel(NbrSitesX, NbrSitesY, NbrSitesZ, 8, TmpChem, 
// 							    Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));

   bool FirstRunFlag = true;
   for (int i = MinKx; i <= MaxKx; ++i)
     {
       for (int j = MinKy; j <= MaxKy; ++j)
 	{
 	  for (int k = MinKz; k <= MaxKz; ++k)
 	    {
 	      cout << "(kx=" << i << ",ky=" << j << ",kz=" << k << ") " << endl;
	      ParticleOnSphere* Space = 0;
	      if (Manager.GetBoolean("boson") == false)
		{
#ifdef __128_BIT_LONGLONG__
		  if (TotalNbrSites <= 15)
#else
		    if (TotalNbrSites <= 7)
#endif
		      {
			Space = new FermionOnCubicLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
		      }
		    else
		      {
			Space = new FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
		      }
		}
	      else
		{
		  Space = new BosonOnCubicLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
		}
	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	      AbstractHamiltonian* Hamiltonian = 0;
	      Hamiltonian = new ParticleOnCubicLatticeFourBandPyrochloreHamiltonian((ParticleOnSphereWithSU4Spin*) Space, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ,
										    Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), &TightBindingModel,
										    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d %d", i, j, k);
	      char* EigenstateOutputFile = new char [512];
	      char* TmpExtention = new char [512];
	      sprintf (TmpExtention, "_kx_%d_ky_%d_kz_%d", i, j, k);
	      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	      
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

