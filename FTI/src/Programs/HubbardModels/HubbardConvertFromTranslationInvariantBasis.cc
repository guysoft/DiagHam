#include "Vector/ComplexVector.h"

#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Options/Options.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("HubbardConvertFromTranslationInvariantBasis" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption('i', "input-state", "name of the file containing the state whose Kx momentum has to be computed");
  (*SystemGroup) += new SingleStringOption('\n', "degenerate-states", "name of the file containing a list of states (override input-state)");
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardConvertFromTranslationInvariantBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int NbrSites = 0;
  int XMomentum = 0;
  int YMomentum = 0;
  int XPeriodicity = 0;
  int YPeriodicity = 0;
  bool Statistics = true;
  bool GutzwillerFlag = false;
  int NbrInputStates = 0;

  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0))
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type HubbardSuperconductorOrderParameter -h" << endl;
      return -1;
    }

  if (Manager.GetString("input-state") != 0)
    {
      NbrInputStates = 1;
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}
      if ((FTIHubbardModelWith1DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, XMomentum, XPeriodicity, Statistics, GutzwillerFlag) == false) && (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, XMomentum, YMomentum, XPeriodicity, YPeriodicity, Statistics, GutzwillerFlag) == false))
	{
	  cout << "error while retrieving system parameters from file name " <<Manager.GetString("input-state")  << endl;
	  return -1;
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrInputStates = DegenerateFile.GetNbrLines();
      if ((FTIHubbardModelWith1DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites, XMomentum, XPeriodicity, Statistics, GutzwillerFlag) == false) && (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites, XMomentum, YMomentum,  XPeriodicity, YPeriodicity, Statistics, GutzwillerFlag) == false))
	{
	  cout << "error while retrieving system parameters from file name " << DegenerateFile(0, 0) << endl;
	  return -1;
	}
    }
  if (YPeriodicity == 0)
    {
      cout << "Convert from kx basis to real space basis" << endl;
      cout << "Nbr particles=" << NbrParticles << " Nbr sites=" << NbrSites << " Kx=" << XMomentum << " Tx=" << XPeriodicity << endl;
    }
  else
    {
      cout << "Convert from kx,ky basis to real space basis" << endl;
      cout << "Nbr particles=" << NbrParticles << " Nbr sites=" << NbrSites << " Kx=" << XMomentum << " Tx=" << XPeriodicity << " Ky=" << YMomentum << " Ty=" << YPeriodicity << endl;
    }
  ComplexVector* InputStates = new ComplexVector[NbrInputStates];
  char** InputStateNames = new char*[NbrInputStates];
  if (Manager.GetString("input-state") != 0)
    {
      InputStateNames[0] = new char [strlen(Manager.GetString("input-state")) + 1];
      strcpy (InputStateNames[0], Manager.GetString("input-state"));
      if (InputStates[0].ReadVector (Manager.GetString("input-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	  return -1;      
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      if (InputStates[0].ReadVector (DegenerateFile(0, 0)) == false)
	{
	  cout << "can't open vector file " << DegenerateFile(0, 0) << endl;
	  return -1;      
	}	  
      InputStateNames[0] = new char [strlen(DegenerateFile(0, 0)) + 1];
      strcpy (InputStateNames[0], DegenerateFile(0, 0));
      for (int i = 1; i < NbrInputStates; ++i)
	{
	  if (InputStates[i].ReadVector (DegenerateFile(0, i)) == false)
	    {
	      cout << "can't open vector file " << DegenerateFile(0, i) << endl;
	      return -1;      
	    }	  
	  if (InputStates[0].GetVectorDimension() != InputStates[i].GetVectorDimension())
	    {
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << "don't have the same  dimension (" << InputStates[0].GetVectorDimension() << " and " << InputStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	  InputStateNames[i] = new char [strlen(DegenerateFile(0, i)) + 1];
	  strcpy (InputStateNames[i], DegenerateFile(0, i));
	}
    }


  ParticleOnTorusWithSpinAndMagneticTranslations* InputSpace = 0;
  if (Statistics == true)
    {
      if (GutzwillerFlag == false)
	if (YPeriodicity == 0)
	  InputSpace = new FermionOnLatticeWithSpinRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, XPeriodicity);
	else
	  InputSpace = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
      else
	if (YPeriodicity == 0)
	  InputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, XPeriodicity);
	else
	  InputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
    }
  if (InputSpace->GetHilbertSpaceDimension() != InputStates[0].GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("left-state")  << " has a wrong dimension (" << InputStates[0].GetVectorDimension() << ", should be " << InputSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  

  ParticleOnSphereWithSpin* OutputSpace = 0;
  if (Statistics == true)
    {
      if (GutzwillerFlag == false)
	OutputSpace = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
      else
	OutputSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
    }

  char* XPeriodicityString = new char[256];
  char* XMomentumString = new char[256];
  if (YPeriodicity == 0)
  {
    sprintf (XPeriodicityString, "_xmomentum_%d", XPeriodicity);
    sprintf (XMomentumString, "_kx_%d", XMomentum);
  }
  else
  {
    sprintf (XPeriodicityString, "_xymomentum_%d_%d", XPeriodicity, YPeriodicity);
    sprintf (XMomentumString, "_kx_%d_ky_%d", XMomentum, YMomentum);
  }
  
  for (int i = 0; i < NbrInputStates; ++i)
    {
      char* TmpVectorOutputName = ReplaceString(InputStateNames[i], XPeriodicityString, "");
      char* VectorOutputName = ReplaceString(TmpVectorOutputName, XMomentumString, "");
      ComplexVector TmpVector = InputSpace->ConvertFromKxKyBasis(InputStates[i], OutputSpace);
      if (TmpVector.WriteVector(VectorOutputName) == false)
	{
	  cout << "error, can't write vector " << VectorOutputName << endl;
	}
      delete[] TmpVectorOutputName;
      delete[] VectorOutputName;
    }
}
