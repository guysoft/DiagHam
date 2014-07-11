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


// Diagionalize the superconductor order parameter matrix
//
// orderParameter = reference to the superconductor order parameter matrix
// output = reference on the output stream
void HubbardSuperconductorOrderParameterMatrixDiagonalize(ComplexMatrix& orderParameter, ostream& output);


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
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerate-leftstates", "single column file describing a set of degenerate left states");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerate-rightstates", "single column file describing a set of degenerate right states");
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
  int NbrLeftStates = 0;
  if ((Manager.GetString("left-state") == 0) && (Manager.GetString("degenerate-leftstates") == 0))
    {
      cout << "error, a left state file should be provided. See man page for option syntax or type HubbardSuperconductorOrderParameter -h" << endl;
      return -1;
    }
  if (Manager.GetString("left-state") != 0)
    {
      NbrLeftStates = 1;
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
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-leftstates")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrLeftStates = DegenerateFile.GetNbrLines();
      if (FTIHubbardModelFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), LeftNbrParticles, LeftNbrSites, LeftStatistics, LeftGutzwillerFlag) == false)
	{
	  cout << "error while retrieving system parameters from file name " <<Manager.GetString("left-state")  << endl;
	  return -1;
	}
    }

  int RightNbrParticles = 0;
  int RightNbrSites = 0;
  bool RightStatistics = true;
  bool RightGutzwillerFlag = false;
  int NbrRightStates = 0;
  if ((Manager.GetString("right-state") == 0) && (Manager.GetString("degenerate-rightstates") == 0))
    {
      cout << "error, a right state file should be provided. See man page for option syntax or type HubbardSuperconductorOrderParameter -h" << endl;
      return -1;
    }
  if (Manager.GetString("right-state") != 0)
    {
      NbrRightStates = 1;
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
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-rightstates")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrRightStates = DegenerateFile.GetNbrLines();
      if (FTIHubbardModelFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), RightNbrParticles, RightNbrSites, RightStatistics, RightGutzwillerFlag) == false)
	{
	  cout << "error while retrieving system parameters from file name " <<  DegenerateFile(0, 0) << endl;
	  return -1;
	}
    }

  if (RightNbrSites != LeftNbrSites)
    {
      cout << "error, left and right states don't have the same number of sites" << endl;
    }

  if (LeftNbrParticles != (RightNbrParticles + 2))
    {
      cout << "error, left and right states don't have the proper number of particles" << endl;
    }

  ComplexVector* LeftStates = new ComplexVector[NbrLeftStates];
  if (Manager.GetString("left-state") != 0)
    {
      if (LeftStates[0].ReadVector (Manager.GetString("left-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("left-state") << endl;
	  return -1;      
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-leftstates")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      if (LeftStates[0].ReadVector (DegenerateFile(0, 0)) == false)
	{
	  cout << "can't open vector file " << DegenerateFile(0, 0) << endl;
	  return -1;      
	}	  
      for (int i = 1; i < NbrLeftStates; ++i)
	{
	  if (LeftStates[i].ReadVector (DegenerateFile(0, i)) == false)
	    {
	      cout << "can't open vector file " << DegenerateFile(0, i) << endl;
	      return -1;      
	    }	  
	  if (LeftStates[0].GetVectorDimension() != LeftStates[i].GetVectorDimension())
	    {
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << "don't have the same  dimension (" << LeftStates[0].GetVectorDimension() << " and " << LeftStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	}
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
  if (LeftSpace->GetHilbertSpaceDimension() != LeftStates[0].GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("left-state")  << " has a wrong dimension (" <<LeftStates[0].GetVectorDimension() << ", should be " << LeftSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  

  ComplexVector* RightStates = new ComplexVector[NbrRightStates];
  if (Manager.GetString("right-state") != 0)
    {
      if (RightStates[0].ReadVector (Manager.GetString("right-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("right-state") << endl;
	  return -1;      
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-rightstates")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      if (RightStates[0].ReadVector (DegenerateFile(0, 0)) == false)
	{
	  cout << "can't open vector file " << DegenerateFile(0, 0) << endl;
	  return -1;      
	}	  
      for (int i = 1; i < NbrRightStates; ++i)
	{
	  if (RightStates[i].ReadVector (DegenerateFile(0, i)) == false)
	    {
	      cout << "can't open vector file " << DegenerateFile(0, i) << endl;
	      return -1;      
	    }	  
	  if (RightStates[0].GetVectorDimension() != RightStates[i].GetVectorDimension())
	    {
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << "don't have the same  dimension (" << RightStates[0].GetVectorDimension() << " and " << RightStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	}
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
  if (RightSpace->GetHilbertSpaceDimension() != RightStates[0].GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("right-state")  << " has a wrong dimension (" << RightStates[0].GetVectorDimension() << ", should be " << RightSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }

  RightSpace->SetTargetSpace(LeftSpace);
    
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      if (Manager.GetString("left-state") != 0)
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
      else
	{
	  MultiColumnASCIIFile DegenerateFile;
	  if (DegenerateFile.Parse(Manager.GetString("degenerate-leftstates")) == false)
	    {
	      DegenerateFile.DumpErrors(cout);
	      return -1;
	    }
	  char* TmpFileName = ReplaceExtensionToFileName(DegenerateFile(0, 0), "vec", "orderparam.dat");
	  if (TmpFileName == 0)
	    {
	      cout << "no vec extension was find in " << DegenerateFile(0, 0) << " file name" << endl;
	      return 0;
	    }
	  File.open(TmpFileName, ios::binary | ios::out);
	  delete[] TmpFileName;
	}
    }

  File.precision(14);
  cout.precision(14);
  if ((NbrLeftStates == 1) && (NbrRightStates == 1))
    {
      File << "# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R> with sigma,sigma' = 0 (down) or 1 (up)" << endl
	   << "# i j sigma sigma' |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>|^2 |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>| Arg(<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>)" << endl;
    }
  else
    {
      File << "# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R> with sigma,sigma' = 0 (down) or 1 (up)" << endl
	   << "# i j sigma sigma' |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>|^2" << endl;
    }
  double NormalizationFactor = 1.0 / sqrt((double) (NbrLeftStates * NbrRightStates));
  for (int i = 0; i < RightNbrSites; ++i)
    {
      int Index = i;
      if (RightStatistics == true)
	{
	  if ((RightGutzwillerFlag == false) && (LeftGutzwillerFlag == false))
	    {
	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorDownUpDiag (RightSpace, i, 0, i, 1);
	      ComplexMatrix TmpMatrix (NbrLeftStates, NbrRightStates);
	      for (int k = 0; k < NbrLeftStates; ++k)
		for (int l = 0; l < NbrRightStates; ++l)
		  {
		    OperatorMatrixElementOperation OperationDownUp(&OperatorDownUpDiag, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		    OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
		    TmpMatrix[l][k] = OperationDownUp.GetScalar() * NormalizationFactor;
		  }
	      File << i << " " << i << " 0 1";
	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	      File << endl;
	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorUpDownDiag (RightSpace, i, 1, i, 0);
	      for (int k = 0; k < NbrLeftStates; ++k)
		for (int l = 0; l < NbrRightStates; ++l)
		  {
		    OperatorMatrixElementOperation OperationUpDown(&OperatorUpDownDiag, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		    OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
		    TmpMatrix[l][k] = OperationUpDown.GetScalar() * NormalizationFactor;
		  }
	      File << i << " " << i << " 1 0";
	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	      File << endl;
	    }
	  ++Index;
	}
      for (int j = Index; j < RightNbrSites; ++j)
	{
	  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorDownDown (RightSpace, i, 0, j, 0);
	  ComplexMatrix TmpMatrix (NbrLeftStates, NbrRightStates);
	  for (int k = 0; k < NbrLeftStates; ++k)
	    for (int l = 0; l < NbrRightStates; ++l)
	      {
		OperatorMatrixElementOperation OperationDownDown(&OperatorDownDown, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		OperationDownDown.ApplyOperation(Architecture.GetArchitecture());
		TmpMatrix[l][k] = OperationDownDown.GetScalar() * NormalizationFactor;
	      }
	  File << i << " " << j << " 0 0 ";
	  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	  File << endl;
	  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorDownUp (RightSpace, i, 0, j, 1);
	  for (int k = 0; k < NbrLeftStates; ++k)
	    for (int l = 0; l < NbrRightStates; ++l)
	      {
		OperatorMatrixElementOperation OperationDownUp(&OperatorDownUp, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
		TmpMatrix[l][k] = OperationDownUp.GetScalar() * NormalizationFactor;
	      }
	  File << i << " " << j << " 0 1";
	  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	  File << endl;
	  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorUpDown (RightSpace, i, 1, j, 0);
	  for (int k = 0; k < NbrLeftStates; ++k)
	    for (int l = 0; l < NbrRightStates; ++l)
	      {
		OperatorMatrixElementOperation OperationUpDown(&OperatorUpDown, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
		TmpMatrix[l][k] = OperationUpDown.GetScalar() * NormalizationFactor;
	      }
	  File << i << " " << j << " 1 0";
	  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	  File << endl;
	  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator OperatorUpUp (RightSpace, i, 1, j, 1);
	  for (int k = 0; k < NbrLeftStates; ++k)
	    for (int l = 0; l < NbrRightStates; ++l)
	      {
		OperatorMatrixElementOperation OperationUpUp(&OperatorUpUp, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
		OperationUpUp.ApplyOperation(Architecture.GetArchitecture());
		TmpMatrix[l][k] = OperationUpUp.GetScalar() * NormalizationFactor;
	      }
	  File << i << " " << j << " 1 1";
	  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
	  File << endl;
	}
    }
  File.close();

  return 0;
}

// Diagionalize the superconductor order parameter matrix
//
// orderParameter = reference to the superconductor order parameter matrix
// output = reference on the output stream

void HubbardSuperconductorOrderParameterMatrixDiagonalize(ComplexMatrix& orderParameter, ostream& output)
{
  if ((orderParameter.GetNbrRow() > 1) || (orderParameter.GetNbrColumn() > 1))
    {
      ComplexMatrix TmpConjugate;
      TmpConjugate.Copy(orderParameter);
      TmpConjugate.HermitianTranspose();
      if (orderParameter.GetNbrRow() >= orderParameter.GetNbrColumn())
	{      
	  HermitianMatrix TmpMatrix (TmpConjugate * orderParameter);
	  RealDiagonalMatrix TmpDiag (TmpMatrix.GetNbrRow());
#ifdef __LAPACK__
	  TmpMatrix.LapackDiagonalize(TmpDiag);
#else
	  TmpMatrix.Diagonalize(TmpDiag);
#endif		  
	  TmpDiag.SortMatrixDownOrder();
	  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
	    output << " " << TmpDiag[i];
	}
      else
	{
	  HermitianMatrix TmpMatrix2 (orderParameter * TmpConjugate);
	  RealDiagonalMatrix TmpDiag2 (TmpMatrix2.GetNbrRow());
#ifdef __LAPACK__
	  TmpMatrix2.LapackDiagonalize(TmpDiag2);
#else
	  TmpMatrix2.Diagonalize(TmpDiag2);
#endif		  
	  TmpDiag2.SortMatrixDownOrder();
	  for (int i = 0; i < TmpDiag2.GetNbrRow(); ++i)
	    output << " " << TmpDiag2[i];
	}
    }
  else
    {
      output << " " << SqrNorm(orderParameter[0][0]) << " " << Norm(orderParameter[0][0]) << Arg(orderParameter[0][0]);
    }
}
