#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"

#include "Operator/ParticleOnSphereWithSpinSuperconductorOrderParameterOperator.h"


#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("HubbardSuperconductorOrderParameter" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "left-state", "name of the file corresponding to the state |Psi_L> (the order parameter being <Psi_L|c^+c^+|Psi_R>");
  (*SystemGroup) += new SingleStringOption  ('\n', "right-state", "name of the file corresponding to the state |Psi_R> (the order parameter being <Psi_L|c^+c^+|Psi_R>");   (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardSuperconductorOrderParameter -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int LeftNbrParticles = 0;
  int LeftNbrSites = 0;
  bool LeftStatistics = true;
  bool LeftGutzwillerFlag = false;
  if (Manager.GetString("left-state") == 0)
    {
      cout << "error, a left state file should be provided. See man page for option syntax or type HubbardSuperconductorOrderParameter -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("left-state")) == false)
    {
      cout << "can't open file " << Manager.GetString("left-state") << endl;
      return -1;
    }
  if (FTIHubbardModelFindSystemInfoFromVectorFileName(Manager.GetString("left-state"), LeftNbrParticles, LeftNbrSites, LeftStatistics, LeftGutzwillerFlag) == false)
    {
      cout << "error while retrieving system parameters from file name " <<Manager.GetString("left-state")  << endl;
      return -1;
    }

  int RightNbrParticles = 0;
  int RightNbrSites = 0;
  bool RightStatistics = true;
  bool RightGutzwillerFlag = false;
  if (Manager.GetString("right-state") == 0)
    {
      cout << "error, a right state file should be provided. See man page for option syntax or type HubbardSuperconductorOrderParameter -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("right-state")) == false)
    {
      cout << "can't open file " << Manager.GetString("right-state") << endl;
      return -1;
    }
  if (FTIHubbardModelFindSystemInfoFromVectorFileName(Manager.GetString("right-state"), RightNbrParticles, RightNbrSites, RightStatistics, RightGutzwillerFlag) == false)
    {
      cout << "error while retrieving system parameters from file name " <<Manager.GetString("right-state")  << endl;
      return -1;
    }

  if (RightNbrSites != LeftNbrSites)
    {
      cout << "error, " << Manager.GetString("left-state")  << " and " << Manager.GetString("right-state") << "don't have the same number of sites" << endl;
    }

  if (LeftNbrParticles != (RightNbrParticles + 2))
    {
      cout << "error, " << Manager.GetString("left-state")  << " and " << Manager.GetString("right-state") << "don't have the proper number of particles" << endl;
    }

  ComplexVector LeftState;
  if (LeftState.ReadVector (Manager.GetString("left-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("left-state") << endl;
      return -1;      
    }
  ParticleOnSphereWithSpin* LeftSpace = 0;
  if (LeftStatistics == true)
    {
      if (LeftGutzwillerFlag == false)
	LeftSpace = new FermionOnLatticeWithSpinRealSpace (LeftNbrParticles, LeftNbrSites);
      else
	LeftSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (LeftNbrParticles, LeftNbrSites);
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
    }
  if (LeftSpace->GetHilbertSpaceDimension() != LeftState.GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("left-state")  << " has a wrong dimension (" <<LeftState.GetVectorDimension() << ", should be " << LeftSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  

  ComplexVector RightState;
  if (RightState.ReadVector (Manager.GetString("right-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("right-state") << endl;
      return -1;      
    }
  ParticleOnSphereWithSpin* RightSpace = 0;
  if (RightStatistics == true)
    {
      if (RightGutzwillerFlag == false)
	RightSpace = new FermionOnLatticeWithSpinRealSpace (RightNbrParticles, RightNbrSites);
      else
	RightSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (RightNbrParticles, RightNbrSites);
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
    }
  if (RightSpace->GetHilbertSpaceDimension() != RightState.GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("right-state")  << " has a wrong dimension (" <<RightState.GetVectorDimension() << ", should be " << RightSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }

  RightSpace->SetTargetSpace(LeftSpace);
    
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(Manager.GetString("left-state"), "vec", "orderparam.dat");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("left-state") << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }

  File.precision(14);
  cout.precision(14);
  File << "# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R> with sigma,sigma' = 0 (down) or 1 (up)" << endl
       << "# i j sigma sigma' <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R> Norm(<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>) Arg(<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>)" << endl;
  for (int i = 0; i < RightNbrSites; ++i)
    {
      int Index = i;
      if (RightStatistics == true)
	{
	  if ((RightGutzwillerFlag == false) && (LeftGutzwillerFlag == false))
	    {
	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorDownUpDiag (RightSpace, i, 0, i, 1);
	      OperatorMatrixElementOperation OperationDownUp(&OperatorDownUpDiag, LeftState, RightState, RightState.GetVectorDimension());
	      OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
	      Complex Tmp = OperationDownUp.GetScalar();
	      File << i << " " << i << " 0 1 " << Tmp << " " << Norm(Tmp) << " " << Arg(Tmp) << endl;
	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorUpDownDiag (RightSpace, i, 1, i, 0);
	      OperatorMatrixElementOperation OperationUpDown(&OperatorUpDownDiag, LeftState, RightState, RightState.GetVectorDimension());
	      OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
	      Tmp = OperationUpDown.GetScalar();
	      File << i << " " << i << " 1 0 " << Tmp << " " << Norm(Tmp) << " " << Arg(Tmp) << endl;
	    }
	  ++Index;
	}
      for (int j = Index; j < RightNbrSites; ++j)
	{
	  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorDownDown (RightSpace, i, 0, j, 0);
	  OperatorMatrixElementOperation OperationDownDown(&OperatorDownDown, LeftState, RightState, RightState.GetVectorDimension());
	  OperationDownDown.ApplyOperation(Architecture.GetArchitecture());
	  Complex Tmp = OperationDownDown.GetScalar();
	  File << i << " " << j << " 0 0 " << Tmp << " " << Norm(Tmp) << " " << Arg(Tmp) << endl;
	  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorDownUp (RightSpace, i, 0, j, 1);
	  OperatorMatrixElementOperation OperationDownUp(&OperatorDownUp, LeftState, RightState, RightState.GetVectorDimension());
	  OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
	  Tmp = OperationDownUp.GetScalar();
	  File << i << " " << j << " 0 1 " << Tmp << " " << Norm(Tmp) << " " << Arg(Tmp) << endl;
	  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorUpDown (RightSpace, i, 1, j, 0);
	  OperatorMatrixElementOperation OperationUpDown(&OperatorUpDown, LeftState, RightState, RightState.GetVectorDimension());
	  OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
	  Tmp = OperationUpDown.GetScalar();
	  File << i << " " << j << " 1 0 " << Tmp << " " << Norm(Tmp) << " " << Arg(Tmp) << endl;
	  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorUpUp (RightSpace, i, 1, j, 1);
	  OperatorMatrixElementOperation OperationUpUp(&OperatorUpUp, LeftState, RightState, RightState.GetVectorDimension());
	  OperationUpUp.ApplyOperation(Architecture.GetArchitecture());
	  Tmp = OperationUpUp.GetScalar();
	  File << i << " " << j << " 1 1 " << Tmp << " " << Norm(Tmp) << " " << Arg(Tmp) << endl;
	}
    }
  File.close();

  return 0;
}
