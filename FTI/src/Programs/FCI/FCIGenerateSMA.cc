#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU3SpinMomentumSpace.h"

#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModelKagomeLattice.h"
#include "Tools/FTITightBinding/TightBindingModelRubyLattice.h"
#include "Tools/FTITightBinding/TightBindingModelAlternativeKagomeLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include "Operator/AbstractOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Operator/ParticleOnLatticeProjectedDensityOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "GeneralTools/FilenameTools.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
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
  OptionManager Manager ("FCIGenerateSMA" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "name of the vector file to which the SMA should be applied");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerate-groundstate", "name of the file that gives the vector files to which the SMA should be applied (in all-bilinear mode)");	
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);  
  (*SystemGroup) += new BooleanOption  ('\n', "all-bilinear", "apply all bilinear operators to the ground state, without summing on the sector");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIGenerateSMA -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if ((Manager.GetString("eigenstate-file") == 0) && (Manager.GetString("degenerate-groundstate") == 0))
    {
      cout << "error, an eigenstate state file should be provided. See man page for option syntax or type FCIGenerateSMA -h" << endl;
      return -1;
    }
  if (Manager.GetString("import-onebody") == 0)
    {
      cout << "error, a file giving the one-body information should be provided. See man page for option syntax or type FCIGenerateSMA -h" << endl;
      return -1;
    }
 
  bool BilinearFlag = Manager.GetBoolean("all-bilinear");
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  int* TotalKx = 0;
  int* TotalKy = 0;
  int NbrSpaces = 1;
  char** GroundStateFiles = 0;
  bool Statistics = true;
  ComplexVector* GroundStates;
  
  if (Manager.GetString("degenerate-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalKx = new int[1];
      TotalKy = new int[1];
      GroundStates = new ComplexVector[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("eigenstate-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("eigenstate-file"));
      if (GroundStates[0].ReadVector (GroundStateFiles[0]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[0] << endl;
	return -1;      
      }	
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, NbrSiteX, NbrSiteY, TotalKx[0], TotalKy[0], Statistics) == false)
      {
	cout << "error while retrieving system parameters from file name " << GroundStateFiles[0] << endl;
	return -1;
      }
    }
  else
    {
	MultiColumnASCIIFile DegeneratedFile;
	if (DegeneratedFile.Parse(Manager.GetString("degenerate-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
	NbrSpaces = DegeneratedFile.GetNbrLines();
	GroundStateFiles = new char* [NbrSpaces];
	TotalKx = new int[NbrSpaces];
	TotalKy = new int[NbrSpaces];
	GroundStates = new ComplexVector[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));	
	   if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }	
	    if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	 }
    }
 
  
  int MinKx = 0;
  int MaxKx = NbrSiteX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrSiteY - 1;
  if (Manager.GetInteger("only-ky") >= 0)
    {						
      MinKy = Manager.GetInteger("only-ky");
      MaxKy = MinKy;
    }
    
  //Define tightbinding model paramaters. This part is for testing only, will be ignored to make the code more generic
  double t1 = 1.0;
  double l1 = 1.0;
  double t2 = 0.0;
  double l2 = 0.0;
  double mu = 0.0;
  double gammaX = 0.0;
  double gammaY = 0.0;
  bool ExportOneBody = true;
  
  
  ParticleOnSphere* SpaceSource = 0;
  
  
  Generic2DTightBindingModel TightBindingModel(Manager.GetString("import-onebody"));
//   TightBindingModelAlternativeKagomeLattice TightBindingModel(NbrSiteX, NbrSiteY, t1, t2, l1, l2, mu, gammaX, gammaY, Architecture.GetArchitecture(), ExportOneBody); 
//   TightBindingModelKagomeLattice TightBindingModel(NbrSiteX, NbrSiteY, t1, t2, l1, l2, mu, gammaX, gammaY, Architecture.GetArchitecture(), ExportOneBody); 
    
//   int nx1 = 3;
//   int ny1 = 0;
//   int nx2 = -3;
//   int ny2 = 4;
//   int offset = 3;
  
//   double tr = 1;
//   double ti = 1;
//   double t1r = -1.4;
//   double t1i = 2.4;
//   double t4 = -1.46;
//   int nx1 = 2;
//   int ny1 = 5;
//   int nx2 = -4;
//   int ny2 = -1;
//   int offset = 13;
//   TightBindingModelRubyLattice TightBindingModel(NbrSiteX, NbrSiteY, nx1, ny1, nx2, ny2, offset, tr, ti, t1r, t1i, t4, mu, gammaX, gammaY, Architecture.GetArchitecture(), ExportOneBody); 
//   cout << MinKx << " " << MinKy << " " << MaxKx << " " << MaxKy << endl;
  AbstractOperator* Projector;
  ParticleOnSphere* SpaceDestination = 0;
  if (BilinearFlag == false)
  {
    cout << "N = " << NbrParticles << " Nx = " << NbrSiteX << " Ny = " << NbrSiteY << " Kx = "<< TotalKx[0] << " Ky = "<< TotalKy[0] <<  endl;
    if (Statistics == true)
      SpaceSource = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx[0], TotalKy[0]);
    else
      SpaceSource = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx[0], TotalKy[0]);
  for (int kx = MinKx; kx <= MaxKx; ++kx)
  {
   for (int ky = MinKy; ky <= MaxKy; ++ky)
   {
    
    int TotalQx = (TotalKx[0] - kx + NbrSiteX) % NbrSiteX;
    int TotalQy = (TotalKy[0] - ky + NbrSiteY) % NbrSiteY;
    cout << "Total Kx = " << TotalQx << " Total Ky = " << TotalQy << endl;
//     cout << kx << " " << ky << endl;
//     cout << TotalKx << " " << TotalKy << " " << TotalQx << " " << TotalQy << endl;
    if (Statistics == true)
      SpaceDestination = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalQx, TotalQy);
    else
      SpaceDestination = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalQx, TotalQy);
    
    Projector = new ParticleOnLatticeProjectedDensityOperator(SpaceSource, SpaceDestination, &TightBindingModel, kx, ky);
    
    char* EigenstateOutputFile;
    char* TmpExtention = new char [512];
    int QxTwoBrillouinZones = (TotalKx[0] - kx + 2*NbrSiteX) % (2*NbrSiteX);
    int QyTwoBrillouinZones = (TotalKy[0] - ky + 2*NbrSiteY) % (2*NbrSiteY);
    sprintf (TmpExtention, "_SMA_kx_%d_ky_%d.vec", QxTwoBrillouinZones, QyTwoBrillouinZones);
    EigenstateOutputFile = ReplaceExtensionToFileName(GroundStateFiles[0], ".vec", TmpExtention);
    ComplexVector EigenstateOutput(SpaceDestination->GetHilbertSpaceDimension(), true);
    
    Projector->LowLevelAddMultiply(GroundStates[0], EigenstateOutput, 0, SpaceSource->GetHilbertSpaceDimension());
    EigenstateOutput.Normalize();
    EigenstateOutput.WriteVector(EigenstateOutputFile);
    
    
    delete[] EigenstateOutputFile;
   
   }
  }
  }
  
  else
  {
    for (int Qx0 = MinKx; Qx0 <= MaxKx; ++Qx0)
    {
      for (int Qy0 = MinKy; Qy0 <= MaxKy; ++Qy0)
      {
	if (Statistics == true)
	  SpaceDestination = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, Qx0, Qy0);
	else
	  SpaceDestination = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, Qx0, Qy0);
	for (int i = 0; i < NbrSpaces; ++i)
	{
	  cout << "N = " << NbrParticles << " Nx = " << NbrSiteX << " Ny = " << NbrSiteY << " Kx = "<< TotalKx[i] << " Ky = "<< TotalKy[i] <<  endl;
	  if (Statistics == true)
	      SpaceSource = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
	  else
	      SpaceSource = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
      
	  SpaceSource->SetTargetSpace(SpaceDestination);
	  for (int kx = 0; kx < NbrSiteX; ++kx)
	    {
	      for (int ky = 0; ky < NbrSiteY; ++ky)
	      {
		int indexDagger = TightBindingModel.GetLinearizedMomentumIndexSafe(kx + Qx0 - TotalKx[i], ky + Qy0 - TotalKy[i]);
		int index = TightBindingModel.GetLinearizedMomentumIndexSafe(kx, ky);
		Projector = new ParticleOnSphereDensityOperator(SpaceSource, indexDagger, index);
// 		cout << kx << " " << ky << " " << Qx0 << " " << Qy0 << endl;
		char* EigenstateOutputFile;
		char* TmpExtention = new char [512];
		sprintf (TmpExtention, "_bilinear_kx_%d_ky_%d_qx0_%d_qy0_%d.vec", kx, ky, Qx0, Qy0);
		EigenstateOutputFile = ReplaceExtensionToFileName(GroundStateFiles[i], ".vec", TmpExtention);
		ComplexVector EigenstateOutput(SpaceDestination->GetHilbertSpaceDimension(), true);
    
		Projector->LowLevelAddMultiply(GroundStates[i], EigenstateOutput, 0, SpaceSource->GetHilbertSpaceDimension());
// 		EigenstateOutput.Normalize();
		EigenstateOutput.WriteVector(EigenstateOutputFile);
    
    
		delete[] EigenstateOutputFile;
	      }
	    }
	  }   
      }
    }
  }
  delete SpaceDestination;
  delete SpaceSource;
  delete Projector;
  
  return 0;
}
