#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphere.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Vector/RealVector.h"

#include "Options/Options.h"
#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <cstdlib>
#include <climits>
#include <cmath>
#include <cstring>
#include <sys/time.h>
#include <stdio.h>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;



int main(int argc, char** argv)
{
  cout.precision(14);
  
  // some running options and help
  OptionManager Manager ("FQHESphereTwoLandauLevelProjection" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* SystemGroup = new OptionGroup("system options");
  
  Manager += MiscGroup;
  Manager +=SystemGroup;
  Manager +=OutputGroup;
  
  (*SystemGroup) += new SingleStringOption  ('\n', "state", "name of the file corresponding to the state to be projected");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particule (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the inital momentum projection for the system (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "max lz value a fermions in the LLL can have (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new BooleanOption  ('n', "normalize", "normalize the vector in the end");
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the result of the projection into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the result of the projection into a text file");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereTwoLandauLevelProjection -h" << endl;
      return -1;
    }

  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");	
  int LandauLevelIndexDifference = 1;
  
  int TotalLz=Manager.GetInteger("total-lz");
  bool FermionFlag = false;
  char * StateName= Manager.GetString("state");
  char* OutputFileName = ((SingleStringOption*) Manager["bin-output"])->GetString();
  char* OutputTxtFileName = ((SingleStringOption*) Manager["txt-output"])->GetString();
  bool Normalize = false; 
  if (Manager.GetBoolean("normalize") == true) 
    Normalize = true;
  
  if (StateName == 0)
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereTwoLandauLevelProjection -h" << endl;
      return -1;
    }
  
  if (IsFile(StateName) == false)
    {
      cout << "state " << StateName << " does not exist or can't be opened" << endl;
      return -1;
    }
  
  RealVector GroundState;
  if (GroundState.ReadVector (StateName) == false)
    {
      cout << "can't open vector file " << StateName << endl;
      return -1;
    }

  
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(StateName,NbrParticles, LzMax, TotalLz, FermionFlag) == false)
    {
      cout << "error while retrieving system parameters from file name " <<StateName  << endl;
      return -1;
    }  
  
  
  if ( FermionFlag ) 
    {
      int LzMaxUp = LzMax + (2 * LandauLevelIndexDifference);
      
      int LzMaxDown = LzMax;
      
      ParticleOnSphereWithSpin* Space = new FermionOnSphereTwoLandauLevels (NbrParticles, TotalLz, LzMaxUp, LzMaxDown);
      
      if (Space->GetHilbertSpaceDimension() != GroundState.GetVectorDimension())
	{
	  cout <<Space->GetHilbertSpaceDimension()<<" "<<GroundState.GetVectorDimension()<<endl;
	  cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
	  return -1;
	}
      
      FermionOnSphere * FinalSpace=new FermionOnSphere(NbrParticles,TotalLz,LzMaxDown);
      
      RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
      
      ((FermionOnSphereTwoLandauLevels*) Space)->ProjectionInTheLowestLevel(GroundState,OutputVector,FinalSpace);
      
      if ( Normalize ) 
	OutputVector.Normalize();
      
      OutputVector.WriteVector(OutputFileName);	
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
    }
  else 
    {
      int LzMaxUp = LzMax + (2 * LandauLevelIndexDifference);
      
      int LzMaxDown = LzMax;
      
      
      ParticleOnSphereWithSpin* Space = new BosonOnSphereTwoLandauLevels(NbrParticles, TotalLz, LzMaxUp, LzMaxDown);
      
      if (Space->GetHilbertSpaceDimension() != GroundState.GetVectorDimension())
	{
	  cout <<Space->GetHilbertSpaceDimension()<<" "<<GroundState.GetVectorDimension()<<endl;
	  cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
	  return -1;
	}
	
      /*for (long i = 0; i < Space->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      cout << GroundState[i] << " ";
	      ((BosonOnSphereTwoLandauLevels*)Space)->PrintStateMonomial(cout, i) << endl;
	    }*/
      
      BosonOnSphereShort * FinalSpace=new BosonOnSphereShort(NbrParticles,TotalLz,LzMaxDown);
      
      RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
      
      ((BosonOnSphereTwoLandauLevels*) Space)->ProjectionInTheLowestLevel(GroundState,OutputVector,FinalSpace);
      
      if ( Normalize ) 
	OutputVector.Normalize();
      
      OutputVector.WriteVector(OutputFileName);	
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
     
    }
    
  return 0;
}
