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
  (*SystemGroup) += new SingleStringOption  ('\n', "lll-state", "name of the file corresponding to the state in the LLL");
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
  
  int LLLNbrParticles = 0;
  int LLLLzMax = 0;
  int LLLTotalLz = 0;
  int NbrFermion = 0;
  int LzMaxFermion = 0;
  int TotalLzFermion = 0;
  bool FermionFlag = false;
  int Step = Manager.GetInteger("step");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  bool Projection = Manager.GetBoolean("projection");
  char* OutputFileName = Manager.GetString("bin-output");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  bool Resume = Manager.GetBoolean("resume");
  int MinComponent = Manager.GetInteger("min-component");
  int NbrComponents = Manager.GetInteger("nbr-components");
  char* LLLFileName = Manager.GetString("lll-state");
  char* FermionFileName = Manager.GetString("fermion");
  bool Flag3LL =  Manager.GetBoolean("3-ll");
  bool Symmetric = Manager.GetBoolean("symmetric");
  bool LLLFermionFlag = true;
  
  
  if (LLLFileName == 0)
    {
      cout << "error, a lowest Landau level state file should be provided. See man page for option syntax or type FQHESphereLLLProjectedLLLTimesManyLL -h" << endl;
      return -1;
    }

  if (FermionFileName == 0)
    {
      cout << "error, a fermionic state file should be provided. See man page for option syntax or type FQHESphereLLLProjectedLLLTimesManyLL -h" << endl;
      return -1;
    }
  
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(LLLFileName, LLLNbrParticles, LLLLzMax, LLLTotalLz, LLLFermionFlag) == false)
    {
      return -1;
    }
  
  
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(FermionFileName, NbrFermion,LzMaxFermion,TotalLzFermion, FermionFlag) == false)
    {
      return -1;
    }
  
 if (LLLNbrParticles != NbrFermion)
    {
      cout << "The number of particles in the two states must be the same" <<endl;
      return -1;
    }

  int Parity = LLLTotalLz & 1;
  if (Parity != ((LLLNbrParticles * LLLLzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;
    }


  if (IsFile(LLLFileName) == false)
    {
      cout << "state " << LLLFileName << " does not exist or can't be opened" << endl;
      return -1;
    }
  if (IsFile(FermionFileName) == false)
    {
      cout << "state " << FermionFileName << " does not exist or can't be opened" << endl;
      return -1;
    }
  
  RealVector LLLState;
  if (LLLState.ReadVector (LLLFileName) == false)
    {
      cout << "can't open vector file " << LLLFileName << endl;
      return -1;
    }
 ParticleOnSphere* LLLSpace = 0;

  if (LLLFermionFlag == false)
    {
      if(HaldaneBasisFlag == false)
        {
#ifdef  __64_BITS__
          if((LLLLzMax + LLLNbrParticles - 1) < 63)
#else
            if ((LLLLzMax + LLLNbrParticles - 1) < 31)
#endif
              {
                LLLSpace = new BosonOnSphereShort (LLLNbrParticles, LLLTotalLz, LLLLzMax);
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
          if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", LLLNbrParticles) == false) || (LLLNbrParticles <= 0))
            {
              cout << "NbrParticles is not defined or as a wrong value" << endl;
              return -1;
            }
      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LLLLzMax) == false) || (LLLLzMax < 0))
        {
          cout << "LzMax is not defined or as a wrong value" << endl;
          return -1;
        }
      int MaxNbrLz;
      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
        {
          cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
          return -1;
        }
      if (MaxNbrLz != (LLLLzMax + 1))
        {
          cout << "wrong LzMax value in ReferenceState" << endl;
          return -1;
        }
#ifdef  __64_BITS__
      if (LLLLzMax  < 63)
#else
        if (LLLLzMax  < 31)
#endif
          LLLSpace = new BosonOnSphereHaldaneBasisShort (LLLNbrParticles,LLLTotalLz, LLLLzMax, ReferenceState);
        }
    }
  else
    {
      if (HaldaneBasisFlag == false)
        {
#ifdef __64_BITS__
          if (LLLLzMax <= 63)
#else
            if (LLLLzMax <= 31)
#endif
              {
                LLLSpace = new FermionOnSphere (LLLNbrParticles,LLLTotalLz,LLLLzMax);
              }
            else
              {
                cout << "This lowest Landau level space requires FermionOnSphereLong class for which this kind of calculation is not available"<<endl;
                return -1;
              }
        }
      else
        {
          int * ReferenceState = 0;
          ConfigurationParser ReferenceStateDefinition;
          if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
            {
              ReferenceStateDefinition.DumpErrors(cout) << endl;
              return -1;
            }
          if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", LLLNbrParticles) == \
               false) || (LLLNbrParticles <= 0))
            {
              cout << "NbrParticles is not defined or as a wrong value" << endl;
              return -1;
            }
          if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LLLLzMax) == false) || (LLLLzMax < 0))
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
          if (MaxNbrLz != (LLLLzMax + 1))
            {
              cout << "wrong LzMax value in ReferenceState" << endl;
              return -1;
            }
#ifdef  __64_BITS__
          if (LLLLzMax  < 63)
#else
            if (LLLLzMax  < 31)
#endif
              LLLSpace = new FermionOnSphereHaldaneBasis (LLLNbrParticles,LLLTotalLz, LLLLzMax, ReferenceState);
        }
    }

  if (LLLSpace->GetHilbertSpaceDimension() != LLLState.GetVectorDimension())
    {
      cout << "Number of rows of the LLL vector is not equal to the Hilbert space dimension!"<<endl;;
      return -1;
    }
  if (Step > LLLSpace->GetHilbertSpaceDimension())
    {
      cout << " The LLL space cannot be divided in less than 1 dimension space"<<endl;
      return -1;
    }

  ParticleOnSphere * SpaceLL=0;
  int LzMaxUp = LzMaxFermion + 2;
  int LzMaxDown = LzMaxFermion;
  if(Flag3LL==true)
    {
      SpaceLL = new FermionOnSphereThreeLandauLevels (LLLNbrParticles, TotalLzFermion, LzMaxFermion);
    }
  else
    {
      SpaceLL = new FermionOnSphereTwoLandauLevels (LLLNbrParticles, TotalLzFermion, LzMaxUp, LzMaxDown);
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
  if (LLLFermionFlag == false)
    {
      if(Projection)
        {
          FinalSpace = new FermionOnSphere (LLLNbrParticles, TotalLzFermion, LzMaxDown + LLLLzMax);
        }
      else
        FinalSpace = new FermionOnSphereTwoLandauLevels (LLLNbrParticles, TotalLzFermion, LzMaxUp + LLLLzMax, LzMaxDown + LLLLzMax);
    }
  else
    {
      FinalSpace = new BosonOnSphereShort (LLLNbrParticles, TotalLzFermion, LzMaxDown + LLLLzMax);
    }

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
  FQHESphereMonomialsTimesSlaterProjectionOperation Operation(SpaceLL, LLLSpace, FinalSpace, &FermionState, &LLLState, OutputVector, MinComponent, NbrComponents, Projection, 
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
