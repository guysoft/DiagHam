#include "Options/Options.h"

#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"

#include "Tools/FTITightBinding/TightBindingModelBilayerTwistBilayer.h"
#include "Tools/FTITightBinding/TightBindingModelTrilayerGrapheneTwistBN.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include "Hamiltonian/ParticleOnSuperlatticeFourBandCoulombHamiltonian.h"
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
  OptionManager Manager ("FQSHFourBandProjectedFromContinuumModel" , "0.01");
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
    
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "kappa", "fraction of dielectric constant", 5.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-voltage", "strength of the electric field", -50.0);
  (*SystemGroup) += new SingleDoubleOption  ('d', "screening-distance", "contribution of screening", 5.0);
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx", "number of points to consider in the x direction for the sampling of the Brillouin zone ", 11);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny", "number of points to consider in the y direction for the sampling of the Brillouin zone ", 11);
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "truncate-qx", "number of points to consider in the x direction for the interaction in momentum space ", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "truncate-qy", "number of points to consider in the y direction for the interaction in momentum space ", 1);
  
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  
  (*SystemGroup) += new SingleStringOption ('\n', "tightbinding-name", "name for the tight-binding model");
  
  (*SystemGroup) += new BooleanOption ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-highsymmetryspectrum", "only compute the one body spectrum, restricting to lines connecting the high symmetry points");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-z2invariant", "compute the z2 invariant of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytheta", "export the one-body topological information (phase of the eigenvalues of the D matrix) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  
  (*SystemGroup) += new BooleanOption ('\n', "disable-flatband", "use dispersion of model instead of flat band");
  
  (*SystemGroup) += new BooleanOption ('\n', "fixed-sz", "fix the Sz value");
  (*SystemGroup) += new BooleanOption ('\n', "fixed-vz", "fix the many-body valley polarization");
  (*SystemGroup) += new SingleIntegerOption ('\n', "sz-value", "twice the fixed Sz value", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "vz-value", "twice the fixed value of valley quantum number", 0);
  
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
      cout << "see man page for option syntax or type FQSHFourBandProjectedFromContinuumModel -h" << endl;
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
  
  bool TiltedFlag = true;
  bool FlatBandFlag = !Manager.GetBoolean("disable-flatband");
  
  Complex** FormFactorPlus = 0;
  Complex** FormFactorMinus = 0;
  
  int NbrSamplesX = Manager.GetInteger("nx");
  int NbrSamplesY = Manager.GetInteger("ny");
  int NbrBZX = Manager.GetInteger("truncate-qx");
  int NbrBZY = Manager.GetInteger("truncate-qy");
  double UVoltage = Manager.GetDouble("u-voltage");
  double UPotential = 36780.0 / Manager.GetDouble("kappa");
  double ScreeningDistance = Manager.GetDouble("screening-distance");
      

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky Sz Vz");  

  char* StatisticPrefix = new char [16];
  sprintf (StatisticPrefix, "fermions");
   
  char* FilePrefix = new char [512 + strlen(StatisticPrefix)];
  if (Manager.GetString("tightbinding-name") == 0)
    sprintf(FilePrefix, "%s_TLGhBN", StatisticPrefix);
  else
    sprintf(FilePrefix, "%s_%s", StatisticPrefix, Manager.GetString("tightbinding-name"));  
  
  char* InteractionPrefix = new char [256];
  if (FlatBandFlag)
  {
      
      sprintf(InteractionPrefix, "flatband_coulomb_n_%d_x_%d_y_%d_nx_%d_ny_%d_U_%f_d_%f_gx_%f_gy_%f",
		       NbrParticles, NbrSitesX, NbrSitesY, NbrSamplesX, NbrSamplesY, UVoltage, ScreeningDistance, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
  }
  else
  {
      sprintf(InteractionPrefix, "coulomb_n_%d_x_%d_y_%d_nx_%d_ny_%d_U_%f_d_%f_kappa_%f_gx_%f_gy_%f", NbrParticles, NbrSitesX, NbrSitesY, NbrSamplesX, NbrSamplesY, UVoltage, ScreeningDistance, Manager.GetDouble("kappa"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
  }
 
  
  
  char* EigenvalueOutputFile = new char [512 + strlen(FilePrefix) + strlen(InteractionPrefix)];
  sprintf(EigenvalueOutputFile, "%s_%s.dat", FilePrefix, InteractionPrefix);
    
 char* TmpExtention = new char [512];
 sprintf (TmpExtention, "_sz_%ld_vz_%ld.dat", Manager.GetInteger("sz-value"), Manager.GetInteger("vz-value"));
 char* TmpEigenvalueOutputFile = EigenvalueOutputFile;
 EigenvalueOutputFile = ReplaceExtensionToFileName(TmpEigenvalueOutputFile, ".dat", TmpExtention);
 
 
 Abstract2DTightBindingModel* TightBindingModel;
 cout.precision(14);
  
  if ((Manager.GetBoolean("singleparticle-spectrum") == true) || (Manager.GetBoolean("singleparticle-highsymmetryspectrum") == true))
    {
      if (Manager.GetBoolean("singleparticle-highsymmetryspectrum") == false)
	{
	  bool ExportOneBody = false;
	  if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("export-onebodytheta") == true) || (Manager.GetBoolean("singleparticle-z2invariant") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
	    ExportOneBody = true;
	  
//       TightBindingModel = new TightBindingModelBilayerTwistBilayer (NbrSitesX, NbrSitesY, NbrSamplesX, NbrSamplesY, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);
      
      TightBindingModel = new TightBindingModelTrilayerGrapheneTwistBN (NbrSitesX, NbrSitesY, NbrSamplesX, NbrSamplesY, UVoltage, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);
	  
      TightBindingModel->WriteAsciiSpectrumColumn(EigenvalueOutputFile);
      int BandIndex = TightBindingModel->GetNbrBands() / 2 - 1;
      double BandSpread = TightBindingModel->ComputeBandSpread(BandIndex);
      double DirectBandGapAbove = TightBindingModel->ComputeDirectBandGap(BandIndex);
      double DirectBandGapBelow = TightBindingModel->ComputeDirectBandGap(BandIndex - 1);
      cout << "Spread = " << BandSpread << "  Direct Gap = " << DirectBandGapBelow << " " << DirectBandGapAbove  << "  Flattening = " << (BandSpread / DirectBandGapAbove) << endl;
      if (Manager.GetBoolean("singleparticle-chernnumber") == true)
        {
            cout << "Chern number = " << TightBindingModel->ComputeChernNumber(BandIndex) << endl;
        }
        
	  if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	    {
	      char* BandStructureOutputFile = new char [512 + strlen(FilePrefix) + strlen(InteractionPrefix)];
	      if (Manager.GetString("export-onebodyname") != 0)
		strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	      else
            {
            if (TiltedFlag == false)
                {
                    sprintf (BandStructureOutputFile, "%s_%s.tightbinding.dat", FilePrefix, InteractionPrefix);
                }
            else
                {
		      
                }
            }
            
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
	  if (Manager.GetBoolean("singleparticle-z2invariant") == true)
	    cout << "Z2 invariant = " << TightBindingModel->ComputeZ2Invariant(2) << endl;
	  
// 	  if (Manager.GetBoolean("export-onebodytheta") == true)
// 	    {
// 	      cout << "Z2 invariant = " << TightBindingModel->ComputeZ2Invariant(2) << endl;
// 	      char* ThetaOutputFile = new char [512];
// 	      sprintf(ThetaOutputFile, "%s_%s_n_%d_x_%d_y_%d_t1_%f_t2_%f_l1_%f_l2_%f_mix12_%f_mix13_%f_mix23_%f_gx_%f_gy_%f_theta.dat", StatisticPrefix, InteractionPrefix, 
// 		      NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), 
// 		      Manager.GetDouble("mixing-12"), Manager.GetDouble("mixing-13"), Manager.GetDouble("mixing-23"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
// 	      TightBindingModel->WriteAsciiDMatrixEigenValues(ThetaOutputFile, 2);
// 	    }
	  return 0;
	}
      
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
  
  if (Manager.GetString("import-onebody") == 0) 
    {
//         TightBindingModel = new TightBindingModelBilayerTwistBilayer (NbrSitesX, NbrSitesY, NbrSamplesX, NbrSamplesY,  Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true);
        
        TightBindingModel = new TightBindingModelTrilayerGrapheneTwistBN (NbrSitesX, NbrSitesY, NbrSamplesX, NbrSamplesY, UVoltage, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true);
        
      char* BandStructureOutputFile = new char [512 + strlen(FilePrefix) + strlen(InteractionPrefix)];
      sprintf (BandStructureOutputFile, "%s_%s.tightbinding.dat", FilePrefix, InteractionPrefix);
//       cout << BandStructureOutputFile << endl;
// 	  TightBindingModel->WriteBandStructure(BandStructureOutputFile);
      delete[] BandStructureOutputFile;
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
	      int MinSz = -NbrParticles;
	      int MaxSz = 0;
	      if (Manager.GetBoolean("fixed-sz") == true)
            {
                MinSz = Manager.GetInteger("sz-value");
                MaxSz = MinSz;
            }
	      
	      
	      int MinVz = -NbrParticles;
	      int MaxVz = 0;
	      if (Manager.GetBoolean("fixed-vz") == true)
            {
                MinVz = Manager.GetInteger("vz-value");
                MaxVz = MinVz;
            }
            
	      for (int Sz = MinSz; Sz <= MaxSz; Sz += 2)
            {
                for (int Vz = MinVz; Vz <= MaxVz; Vz += 2)
                {
                    ParticleOnSphereWithSU4Spin* Space = 0;
                    cout << "(kx=" << i << ",ky=" << j << ") Sz=" << Sz << " -- Vz=" << Vz <<  endl;
                    
                    if (NbrSitesX * NbrSitesY <= 16)
                        Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j, Sz, Vz);
                    else
                        Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j, Sz, Vz);
                    
                    cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
                    if (Space->GetHilbertSpaceDimension() > 0)
                    {
                        Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
                        AbstractQHEHamiltonian* Hamiltonian = 0;
                        
                        Hamiltonian = new ParticleOnSuperlatticeFourBandCoulombHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY, NbrBZX, NbrBZY, UPotential, ScreeningDistance, 
										  Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), TightBindingModel,
															FlatBandFlag, Architecture.GetArchitecture(), Memory);
                           
                        char* ContentPrefix = new char[256];
                        sprintf (ContentPrefix, "%d %d %d %d", i, j, Sz, Vz);
                        char* EigenstateOutputFile = new char [512];
                        char* TmpExtention = new char [512];
                        sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d_vz_%d", i, j, Sz, Vz);
                        EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
                        GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
                        FirstRunFlag = false;
                        MainTaskOperation TaskOperation (&Task);
                        TaskOperation.ApplyOperation(Architecture.GetArchitecture());
                        cout << "------------------------------------" << endl;
                        delete Hamiltonian;
                        delete[] EigenstateOutputFile;
                        delete[] ContentPrefix;
                    }
                    delete Space;
		}
	    }
	}
    }
  return 0;
}

