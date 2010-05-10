#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphereThreeLandauLevels.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"


#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Operator/ParticleOnSphereDensityDensityOperator.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereMonomialsTimesSlaterProjectionOperation.h"

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
  OptionManager Manager ("FQHESphereLLLProjectedLLLTimesTwoLL" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  //OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  //OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  Manager += SystemGroup;
  //Manager += PrecalculationGroup;
  Manager += OutputGroup;
  //Manager += ToolsGroup;
  Manager += MiscGroup;
  ArchitectureManager Architecture;
  Architecture.AddOptionGroup(&Manager);
  
  (*SystemGroup) += new SingleStringOption  ('\n', "fermion", "name of the file corresponding to the fermion state in the Two Lowest Landau Level");
  (*SystemGroup) += new SingleStringOption  ('\n', "boson", "name of the file corresponding to the boson state in the LLL");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state (should be the one of the bosonic state)");
  (*SystemGroup) += new BooleanOption ('\n',"projection","the state will be projected into the LLL");
  (*SystemGroup) += new BooleanOption ('\n', "resume", "the last calcul will be resumed from its last save step");
  (*SystemGroup) += new BooleanOption  ('\n', "3-ll", "consider particles within three Landau levels");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetric", "Lz->-Lz symmetry used");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-component", "the min component", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-components", "the number of component computed", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n',"step", "number of time the Bosonic will be divided in",1);
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "express the projected state in the normalized basis");
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the Jack polynomial decomposition into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the Jack polynomial decomposition into a text file");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereLLLProjectedLLLTimesManyLL -h" << endl;
      return -1;
    }	
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrBoson=0;
  int LzMaxBoson=0;
  int TotalLzBoson=0;
  int NbrFermion=0;
  int LzMaxFermion=0;
  int TotalLzFermion=0;
  bool FermionFlag = false;
  int Step = Manager.GetInteger("step");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  bool Projection = Manager.GetBoolean("projection");
  char* OutputFileName = Manager.GetString("bin-output");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  bool Resume = Manager.GetBoolean("resume");
  int MinComponent = Manager.GetInteger("min-component");
  int NbrComponents = Manager.GetInteger("nbr-components");
  char* BosonFileName = Manager.GetString("boson");
  char* FermionFileName = Manager.GetString("fermion");
  bool Flag3LL =  Manager.GetBoolean("3-ll");
  bool Symmetric = Manager.GetBoolean("symmetric");
  
  
  if (BosonFileName == 0)
    {
      cout << "error, a bosonic state file should be provided. See man page for option syntax or type FQHESphereLLLProjectedLLLTimesManyLL -h" << endl;
      return -1;
    }
  
  if (FermionFileName == 0)
    {
      cout << "error, a fermionic state file should be provided. See man page for option syntax or type FQHESphereLLLProjectedLLLTimesManyLL -h" << endl;
      return -1;
    }
  
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(BosonFileName, NbrBoson,LzMaxBoson,TotalLzBoson, FermionFlag) == false)
    {		
      return -1;
    }
  
  
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(FermionFileName, NbrFermion,LzMaxFermion,TotalLzFermion, FermionFlag) == false)
    {
      return -1;
    }
  
  if(NbrBoson!=NbrFermion)
    {
      cout << "The number of bosons and of fermions must be the same" <<endl;
      return -1;
    }
  
  int Parity = TotalLzBoson & 1;
  if (Parity != ((NbrBoson * LzMaxBoson) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;
    }
  
  
  if (IsFile(BosonFileName) == false)
    {
      cout << "state " << BosonFileName << " does not exist or can't be opened" << endl;
      return -1;
    }
  if (IsFile(FermionFileName) == false)
    {
      cout << "state " << FermionFileName << " does not exist or can't be opened" << endl;
      return -1;
    }
  
  RealVector BosonState;
  if (BosonState.ReadVector (BosonFileName) == false)
    {
      cout << "can't open vector file " << BosonFileName << endl;
      return -1;      
    }
  
  BosonOnSphereShort* BosonSpace = 0;
  
  if(HaldaneBasisFlag==false)
    {
#ifdef  __64_BITS__
      if((LzMaxBoson + NbrBoson - 1) < 63)
#else
	if ((LzMaxBoson + NbrBoson - 1) < 31)	
#endif
	  {
	    BosonSpace = new BosonOnSphereShort (NbrBoson, TotalLzBoson, LzMaxBoson);
	  }
    }
  else
    {
      int* ReferenceState = 0;
      if (Manager.GetString("reference-file") == 0)
	{
	  cout << "error, a reference file is needed for bosons in Haldane basis" << endl;
	  return -1;
	}
      ConfigurationParser ReferenceStateDefinition;
      if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
	{
	  ReferenceStateDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrBoson) == false) || (NbrBoson <= 0))
	{
	  cout << "NbrParticles is not defined or as a wrong value" << endl;
	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMaxBoson) == false) || (LzMaxBoson < 0))
	{
	  cout << "LzMax is not defined or as a wrong value" << endl;
	  return 0;
	}
      int MaxNbrLz;
      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	{
	  cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
	  return -1;     
	}
      if (MaxNbrLz != (LzMaxBoson + 1))
	{
	  cout << "wrong LzMax value in ReferenceState" << endl;
	  return -1;
	}
#ifdef  __64_BITS__
      if (LzMaxBoson  < 63)
#else
	if (LzMaxBoson  < 31)	
#endif
	  BosonSpace = new BosonOnSphereHaldaneBasisShort(NbrBoson,TotalLzBoson, LzMaxBoson, ReferenceState);
    }
  
  if (BosonSpace->GetHilbertSpaceDimension() != BosonState.GetVectorDimension())
    {
      cout << "Number of rows of the bosonic vector is not equal to the Hilbert space dimension!"<<endl;;
      return -1;
    }
  if (Step>BosonSpace->GetHilbertSpaceDimension())
    {
      cout << " The bosonic space cannot be divided in less than 1 dimension space"<<endl;
      return -1;
    }
  ParticleOnSphere * SpaceLL=0;
  int LzMaxUp = LzMaxFermion + 2;
  int LzMaxDown = LzMaxFermion;
  if(Flag3LL==true)
    {
      SpaceLL = new FermionOnSphereThreeLandauLevels (NbrBoson, TotalLzFermion,LzMaxFermion);
    }
  else
    {
      SpaceLL = new FermionOnSphereTwoLandauLevels (NbrBoson, TotalLzFermion, LzMaxUp, LzMaxDown);
    }
  RealVector FermionState;
  
  if (FermionState.ReadVector (FermionFileName) == false)
    {
      cout << "can't open vector file " << FermionFileName << endl;
      return -1;
    }
  
  if(SpaceLL->GetHilbertSpaceDimension()!=FermionState.GetVectorDimension())
    {
      cout <<"Number of rows of the fermionic vector is not equal to the Hilbert space dimension!" <<endl;
      return -1;
    }
  
  ParticleOnSphere * FinalSpace;
  if(Projection)
    {
      FinalSpace = new FermionOnSphere (NbrBoson, TotalLzFermion, LzMaxDown+LzMaxBoson);
    }
  else
    FinalSpace = new FermionOnSphereTwoLandauLevels (NbrBoson, TotalLzFermion, LzMaxUp+LzMaxBoson, LzMaxDown+LzMaxBoson);
  
  RealVector * OutputVector = 0;
  if(Resume)
    {
      char * ResumeVectorName = "temporary_projection_vector.vec";
      char * LogFile = "projection.dat";
      if (IsFile(LogFile) == false)
	{
	  cout << "The calcul cannot be resume as the LogFile doesn't exist." << endl;
	  return -1;
	}
      ifstream File;
      File.open(LogFile , ios::in);
      
      if (!File.is_open())
	{
	  cout << "Cannot open the file: " << LogFile<< endl;
	  return -1;
	}
      File.seekg (0, ios::beg);
      File >> MinComponent;
      File.close();
      if (IsFile(ResumeVectorName) == false)
	{
	  cout << "The calcul cannot be resume as there is no vector saved." << endl;
	  return -1;
	}
      OutputVector = new RealVector;
      if (OutputVector->ReadVector (ResumeVectorName) == false)
	{
	  cout << "can't open vector file " << ResumeVectorName << endl;
	  return -1;
	}
      if(FinalSpace->GetHilbertSpaceDimension()!=OutputVector->GetVectorDimension())
	{
	  cout <<"Number of rows of the resume vector is not equal to the Hilbert space dimension!" <<endl;
	  return -1;
	}
    }
  else
    {
      OutputVector = new RealVector(FinalSpace->GetHilbertSpaceDimension(),true);
    }
  FQHESphereMonomialsTimesSlaterProjectionOperation Operation(SpaceLL, BosonSpace, FinalSpace, &FermionState, &BosonState,OutputVector, MinComponent, NbrComponents, Projection, 
							      Step, Flag3LL, Symmetric);																
  Operation.ApplyOperation(Architecture.GetArchitecture());
  
  if(NbrComponents+MinComponent!=FinalSpace->GetHilbertSpaceDimension())
    {
      char * OutputFileName = "temporary_projection_vector.vec";
      char * LogFile = "projection.dat";
      (*OutputVector).WriteVector(OutputFileName);
      ofstream File;
      File.open(LogFile, ios::binary | ios::out);
      File.precision(14);
      File << NbrComponents+MinComponent;
      File.close();
    }
  
  if(Projection)
    {
      if(Manager.GetBoolean("normalize"))
	{
	  FinalSpace->ConvertFromUnnormalizedMonomial((*OutputVector),0,true);
	}
    }
  
  (*OutputVector).WriteVector(OutputFileName);
  ofstream File;
  if(OutputTxtFileName!=0)
    {
      File.open(OutputTxtFileName, ios::binary | ios::out);
      File.precision(14);
      for (long i = 0; i < FinalSpace->GetLargeHilbertSpaceDimension(); ++i)
	{
	  File << OutputVector[i] << " ";
	  FinalSpace->PrintStateMonomial(File, i) << endl;
	}
    }
  File.close();
  return 0;
}
