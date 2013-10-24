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
#include "Tools/FTITightBinding/TightBindingModelAlternativeKagomeLattice.h"

#include "Operator/AbstractOperator.h"
#include "Operator/ParticleOnLatticeProjectedDensityOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "GeneralTools/FilenameTools.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);  
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
  
  if (Manager.GetString("eigenstate-file") == 0) 
    {
      cout << "error, an eigenstate state file should be provided. See man page for option syntax or type FCIGenerateSMA -h" << endl;
      return -1;
    }
    
 
    
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  int TotalKx = 0;
  int TotalKy = 0;
  bool Statistics = true;
  ComplexVector GroundState;
   
  char* GroundStateFile;
  GroundStateFile = new char [strlen(Manager.GetString("eigenstate-file")) + 1];
  strcpy (GroundStateFile, Manager.GetString("eigenstate-file"));
  
   if (GroundState.ReadVector (GroundStateFile) == false)
    {
      cout << "can't open vector file " << GroundStateFile << endl;
      return -1;      
    }
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFile, NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << GroundStateFile << endl;
      return -1;
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
  if (Statistics == true)
    SpaceSource = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
  else
    SpaceSource = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
  
  
  cout << NbrParticles << " " << NbrSiteX << " " << NbrSiteY << " "<< TotalKx << " "<< TotalKy <<  endl;
  TightBindingModelAlternativeKagomeLattice TightBindingModel(NbrSiteX, NbrSiteY, t1, t2, l1, l2, mu, gammaX, gammaY, Architecture.GetArchitecture(), ExportOneBody); 
//   TightBindingModelKagomeLattice TightBindingModel(NbrSiteX, NbrSiteY, t1, t2, l1, l2, mu, gammaX, gammaY, Architecture.GetArchitecture(), ExportOneBody); 
    
//   int nx1 = 3;
//   int ny1 = 0;
//   int nx2 = -3;
//   int ny2 = 4;
//   int offset = 3;
  
  int nx1 = 2;
  int ny1 = 0;
  int nx2 = -2;
  int ny2 = 3;
  int offset = 2;
//   TightBindingModelAlternativeKagomeLattice TightBindingModel(NbrSiteX, NbrSiteY, nx1, ny1, nx2, ny2, offset, t1, t2, l1, l2, mu, gammaX, gammaY, Architecture.GetArchitecture(), ExportOneBody); 
//   cout << MinKx << " " << MinKy << " " << MaxKx << " " << MaxKy << endl;
  AbstractOperator* Projector;
  ParticleOnSphere* SpaceDestination = 0;
  for (int kx = MinKx; kx <= MaxKx; ++kx)
  {
   for (int ky = MinKy; ky <= MaxKy; ++ky)
   {
    
    int TotalQx = (TotalKx - kx + NbrSiteX) % NbrSiteX;
    int TotalQy = (TotalKy - ky + NbrSiteY) % NbrSiteY;
//     cout << kx << " " << ky << endl;
//     cout << TotalKx << " " << TotalKy << " " << TotalQx << " " << TotalQy << endl;
    if (Statistics == true)
      SpaceDestination = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalQx, TotalQy);
    else
      SpaceDestination = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalQx, TotalQy);
    
    Projector = new ParticleOnLatticeProjectedDensityOperator(SpaceSource, SpaceDestination, &TightBindingModel, kx, ky);
    
    char* EigenstateOutputFile;
    char* TmpExtention = new char [512];
    sprintf (TmpExtention, "_SMA_kx_%d_ky_%d.vec", kx, ky);
    EigenstateOutputFile = ReplaceExtensionToFileName(GroundStateFile, ".vec", TmpExtention);
    ComplexVector EigenstateOutput(SpaceDestination->GetHilbertSpaceDimension(), true);
    
    Projector->LowLevelAddMultiply(GroundState, EigenstateOutput, 0, SpaceSource->GetHilbertSpaceDimension());
    EigenstateOutput.WriteVector(EigenstateOutputFile);
    
    
    delete[] EigenstateOutputFile;
   
   }
  }

  delete SpaceDestination;
  delete SpaceSource;
  delete Projector;
  
  return 0;
}
