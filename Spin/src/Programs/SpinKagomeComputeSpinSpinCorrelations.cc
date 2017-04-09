#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "HilbertSpace/Spin1_2ChainNewAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainNewSzSymmetryAnd2DTranslation.h"

#include "Operator/AbstractOperator.h"
#include "Operator/SpinWith2DTranslationSpinSpinCorrelationOperator.h"


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


// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

int GetLinearizedIndex(int xPosition, int yPosition, int atomicIndex, int nbrSpinX, int nbrSpinY);

int main(int argc, char** argv)
{
  OptionManager Manager ("SpinKagomeComputeSpinSpinCorrelations" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  Manager += PrecalculationGroup;

  
  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "name of the file corresponding to the state |Psi_R> (the order parameter being <Psi_L|c^+c^+|Psi_R>");   (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
//   (*SystemGroup) += new BooleanOption ('\n', "only-cc", "compute only the parameters c^+_sigma c^+_sigma' instead of their linear combinations");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Sz symmetry)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Sz symmetry)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinKagomeComputeSpinSpinCorrelations -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSites = 0;
  int XMomentum = 0;
  int YMomentum = 0;
  int InversionSector = 0;
  int SzValue = 0;
  int SzParitySector = 0;
  bool TotalSpinConservedFlag = true;
  bool InversionFlag = false;
  bool SzSymmetryFlag = false;
  int XPeriodicity = 0;
  int YPeriodicity = 0;
  bool Statistics = true;
  int NbrInputStates = 0;
  int SpinValue = 1;
  int Offset = 0;
  
  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type SpinSystemConvertFromTranslationInvariantBasis -h" << endl;
      return -1;
    }

  else
    {
      InversionSector = 0;
      SzParitySector = 0;
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}   
	
      SpinFindSystemInfoFromVectorFileName (Manager.GetString("input-state"), NbrSites, SzValue, SpinValue, XMomentum, InversionSector, SzParitySector, Offset);
      
//       (char* filename, int& nbrSpins, int& sz, int& spin, int& momentum, int& inversion, int& szSymmetry, int& offset)
      cout << Offset << endl;
      if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue, SpinValue, XMomentum, XPeriodicity,
								YMomentum, YPeriodicity) == false)
	{
	  TotalSpinConservedFlag = false;
	  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum, XPeriodicity, 
								    YMomentum, YPeriodicity) == false)
	    {
	      cout << "error while retrieving system parameters from file name " <<Manager.GetString("input-state")  << endl;
	      return -1;
	    }
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum, XPeriodicity, 
											 YMomentum, YPeriodicity, InversionSector);
	  if (InversionFlag == false)
	    cout << "error" << endl;
	  else
	    cout << "OK" << " " << InversionSector << endl;
	}
      else
	{
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue, SpinValue, XMomentum, XPeriodicity,
											 YMomentum, YPeriodicity, InversionSector);
	  if (SzValue == 0)
	  {
	    SpinFindSystemInfoFromVectorFileName((Manager.GetString("input-state")), NbrSites, SzValue, SpinValue, InversionSector, SzParitySector);
	    cout << NbrSites << " " << SzValue << " " << SpinValue << " " << XMomentum << " " << XPeriodicity << " " << YMomentum << " " <<  YPeriodicity << " " << SzParitySector << endl;
	    if (SzParitySector != 0)
	      SzSymmetryFlag = true;
	  }
	}
    }
   
      
  ComplexVector State;
  if (Manager.GetString("input-state") != 0)
    {
      if (State.ReadVector (Manager.GetString("input-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	  return -1;      
	}
    }

  Spin1_2ChainNewAnd2DTranslation* Space = 0;
      
   if (SzSymmetryFlag == false)
    Space = new Spin1_2ChainNewAnd2DTranslation (NbrSites, SzValue, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
   else
   {
      cout << "Create HilbertSpace" << endl;
      if (Manager.GetString("load-hilbert") != 0)
	{
	  Space = new Spin1_2ChainNewSzSymmetryAnd2DTranslation(Manager.GetString("load-hilbert"));
	}
      else
	Space = new Spin1_2ChainNewSzSymmetryAnd2DTranslation (NbrSites, SzValue, (1 - SzParitySector)/2, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
   }

  

  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("input-state")  << " has a wrong dimension (" << State.GetVectorDimension() << ", should be " << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  
  if (Manager.GetString("save-hilbert") != 0)
  {
    Space->WriteHilbertSpace(Manager.GetString("save-hilbert"));
    return 0;
  }
 
  ofstream File;
  char* OutputFileName;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      File.open(OutputFileName, ios::binary | ios::out);
    }
  else
    {
      if (Manager.GetString("input-state") != 0)
	{
	  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("input-state"), "vec", "neighboring_spin_correlation.dat");
	  if (OutputFileName == 0)
	    {
	      cout << "no vec extension was find in " << Manager.GetString("input-state") << " file name" << endl;
	      return 0;
	    }
	  File.open(OutputFileName, ios::binary | ios::out);
	}
    }
  File.precision(14);
  cout.precision(14);
  
  cout << "NbrSites = " << NbrSites << " " << " XPeriodicity = " << XPeriodicity << " YPeriodicity  = " << YPeriodicity << " Sz = " << SzValue << " SzParitySector = " << SzParitySector << " XMomentum = " << XMomentum << " YMomentum = " << YMomentum << endl;
  
  AbstractOperator* Operator;
  Complex* NeighborSpinSpinCorrelation = new Complex [6];
  int** TmpIndices = new int* [6];
  for (int i = 0; i < 6; ++i)
  {
    TmpIndices[i] = new int[2];
    NeighborSpinSpinCorrelation[i] = 0.0;
  }
  
  File << "# E_{AB} E_{BA} E_{BC} E_{CB} E_{AC} E_{CA}" << endl;
  cout << "# E_{AB} E_{BA} E_{BC} E_{CB} E_{AC} E_{CA}" << endl;
  for (int i = 0; i < XPeriodicity; ++i)
  {
    for (int j = 0; j < YPeriodicity; ++j)
    {
      TmpIndices[0][0] = GetLinearizedIndex(i, j, 0, XPeriodicity, YPeriodicity);
      TmpIndices[0][1] = GetLinearizedIndex(i, j, 1, XPeriodicity, YPeriodicity);
      TmpIndices[1][0] = GetLinearizedIndex(i, j, 1, XPeriodicity, YPeriodicity);
      TmpIndices[1][1] = GetLinearizedIndex(i + 1, j, 0, XPeriodicity, YPeriodicity);
      TmpIndices[2][0] = GetLinearizedIndex(i, j, 1, XPeriodicity, YPeriodicity);
      TmpIndices[2][1] = GetLinearizedIndex(i, j, 2, XPeriodicity, YPeriodicity);
      TmpIndices[3][0] = GetLinearizedIndex(i, j, 2, XPeriodicity, YPeriodicity);
      TmpIndices[3][1] = GetLinearizedIndex(i - 1 + Offset, j + 1, 1, XPeriodicity, YPeriodicity);
      TmpIndices[4][0] = GetLinearizedIndex(i, j, 0, XPeriodicity, YPeriodicity);
      TmpIndices[4][1] = GetLinearizedIndex(i, j, 2, XPeriodicity, YPeriodicity);
      TmpIndices[5][0] = GetLinearizedIndex(i, j, 2, XPeriodicity, YPeriodicity);
      TmpIndices[5][1] = GetLinearizedIndex(i + Offset, j + 1, 0, XPeriodicity, YPeriodicity);
  
//       cout << i << " " << j << " " ;
//       File << i << " " << j << " " ;
      for (int l = 0; l < 6; ++l)
      {
	Operator = new SpinWith2DTranslationSpinSpinCorrelationOperator(Space, XMomentum, XPeriodicity, YMomentum, YPeriodicity, TmpIndices[l][0], TmpIndices[l][1]);
	OperatorMatrixElementOperation Operation(Operator, State, State, State.GetVectorDimension());
	Operation.ApplyOperation(Architecture.GetArchitecture());
	NeighborSpinSpinCorrelation[l] += Operation.GetScalar();
// 	cout << NeighborSpinSpinCorrelation[l] << " " ;
	delete Operator;
      }
//       cout << endl;
      
    }
  }
  
  Complex NematicOP = (-NeighborSpinSpinCorrelation[0] + NeighborSpinSpinCorrelation[1] + NeighborSpinSpinCorrelation[2] * Phase( - M_PI / 3.0 ) + NeighborSpinSpinCorrelation[3] *  Phase( 2.0 * M_PI / 3.0 ) + NeighborSpinSpinCorrelation[4]  * Phase( M_PI / 3.0 ) + NeighborSpinSpinCorrelation[5]  * Phase( -2.0 * M_PI / 3.0 )) / (XPeriodicity * YPeriodicity);
  for (int l = 0; l < 6; ++l)
  {
    File << NeighborSpinSpinCorrelation[l] / (XPeriodicity * YPeriodicity) << " " ;
    cout << NeighborSpinSpinCorrelation[l] / (XPeriodicity * YPeriodicity) << " " ;
  }
  File << endl;
  cout << endl;
  File << "C3 OP = " << sqrt(NematicOP.Re*NematicOP.Re + NematicOP.Im*NematicOP.Im) << endl;
  File.close();
  cout << "C3 OP = " << sqrt(NematicOP.Re*NematicOP.Re + NematicOP.Im*NematicOP.Im) << endl;

//   ofstream FileFourierTransform;
//   if ((RightMomentumFlag == true) && (LeftMomentumFlag == true))
//     {
//       char* TmpFileName = ReplaceExtensionToFileName(OutputFileName, "dat", "fourier.dat");
//       if (TmpFileName == 0)
// 	{
// 	  cout << "no dat extension was find in " << OutputFileName << " file name" << endl;
// 	  return 0;
// 	}
//       FileFourierTransform.open(TmpFileName, ios::binary | ios::out);
//       FileFourierTransform.precision(14);
//     }

//   if (Manager.GetBoolean("only-cc") == true)
//     {
//       if ((NbrLeftStates == 1) && (NbrRightStates == 1))
// 	{
// 	  File << "# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R> with sigma,sigma' = 0 (down) or 1 (up)" << endl
// 	       << "# i j sigma sigma' |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>|^2 |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>| Arg(<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>)" << endl;
// 	}
//       else
// 	{
// 	  File << "# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R> with sigma,sigma' = 0 (down) or 1 (up)" << endl
// 	       << "# i j sigma sigma' |<Psi_L| c^+_{i,sigma} c^+_{j,sigma'} |Psi_R>|^2" << endl;
// 	}
//       double NormalizationFactor = 1.0 / sqrt((double) (NbrLeftStates * NbrRightStates));
//       for (int i = 0; i < RightNbrSites; ++i)
// 	{
// 	  int Index = i;
// 	  if (RightStatistics == true)
// 	    {
// 	      if ((RightGutzwillerFlag == false) && (LeftGutzwillerFlag == false))
// 		{
// 		  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownUpDiag = 0;
// 		  if (RightMomentumFlag == false)
// 		    OperatorDownUpDiag = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator (RightSpace, i, 0, i, 1);
// 		  else
// 		    OperatorDownUpDiag = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 0, i, 1);
// 		  ComplexMatrix TmpMatrix (NbrLeftStates, NbrRightStates);
// 		  for (int k = 0; k < NbrLeftStates; ++k)
// 		    for (int l = 0; l < NbrRightStates; ++l)
// 		      {
// 			OperatorMatrixElementOperation OperationDownUp(OperatorDownUpDiag, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 			OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
// 			TmpMatrix[l][k] = OperationDownUp.GetScalar() * NormalizationFactor;
// 		      }
// 		  File << i << " " << i << " 0 1";
// 		  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 		  File << endl;
// 		  delete OperatorDownUpDiag;
// 		  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpDownDiag = 0;
// 		  if (RightMomentumFlag == false)
// 		    OperatorUpDownDiag = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator (RightSpace, i, 1, i, 0);
// 		  else
// 		    OperatorUpDownDiag = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, i, 0);
// 		  for (int k = 0; k < NbrLeftStates; ++k)
// 		    for (int l = 0; l < NbrRightStates; ++l)
// 		      {
// 			OperatorMatrixElementOperation OperationUpDown(OperatorUpDownDiag, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 			OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
// 			TmpMatrix[l][k] = OperationUpDown.GetScalar() * NormalizationFactor;
// 		      }
// 		  File << i << " " << i << " 1 0";
// 		  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 		  File << endl;
// 		  delete OperatorUpDownDiag;
// 		}
// 	      ++Index;
// 	    }
// 	  for (int j = Index; j < RightNbrSites; ++j)
// 	    {
// 	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownDown = 0;
// 	      if (RightMomentumFlag == false)
// 		OperatorDownDown = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 0, j, 0);
// 	      else
// 		OperatorDownDown = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 0, j, 0);
// 	      ComplexMatrix TmpMatrix (NbrLeftStates, NbrRightStates);
// 	      for (int k = 0; k < NbrLeftStates; ++k)
// 		for (int l = 0; l < NbrRightStates; ++l)
// 		  {
// 		    OperatorMatrixElementOperation OperationDownDown(OperatorDownDown, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 		    OperationDownDown.ApplyOperation(Architecture.GetArchitecture());
// 		    TmpMatrix[l][k] = OperationDownDown.GetScalar() * NormalizationFactor;
// 		  }
// 	      File << i << " " << j << " 0 0";
// 	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 	      File << endl;
// 	      delete OperatorDownDown;
// 	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownUp = 0;
// 	      if (RightMomentumFlag == false)
// 		OperatorDownUp = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 0, j, 1);
// 	      else
// 		OperatorDownUp = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 0, j, 1);
// 	      for (int k = 0; k < NbrLeftStates; ++k)
// 		for (int l = 0; l < NbrRightStates; ++l)
// 		  {
// 		    OperatorMatrixElementOperation OperationDownUp(OperatorDownUp, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 		    OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
// 		    TmpMatrix[l][k] = OperationDownUp.GetScalar() * NormalizationFactor;
// 		  }
// 	      File << i << " " << j << " 0 1";
// 	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 	      File << endl;
// 	      delete OperatorDownUp;
// 	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpDown = 0;
// 	      if (RightMomentumFlag == false)
// 		OperatorUpDown = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 0);
// 	      else
// 		OperatorUpDown = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 0);
// 	      for (int k = 0; k < NbrLeftStates; ++k)
// 		for (int l = 0; l < NbrRightStates; ++l)
// 		  {
// 		    OperatorMatrixElementOperation OperationUpDown(OperatorUpDown, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 		    OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
// 		    TmpMatrix[l][k] = OperationUpDown.GetScalar() * NormalizationFactor;
// 		  }
// 	      File << i << " " << j << " 1 0";
// 	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 	      File << endl;
// 	      delete OperatorUpDown;
// 
// 	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpUp = 0;
// 	      if (RightMomentumFlag == false)
// 		OperatorUpUp = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 1);
// 	      else
// 		OperatorUpUp = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 1);
// 	      for (int k = 0; k < NbrLeftStates; ++k)
// 		for (int l = 0; l < NbrRightStates; ++l)
// 		  {
// 		    OperatorMatrixElementOperation OperationUpUp(OperatorUpUp, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 		    OperationUpUp.ApplyOperation(Architecture.GetArchitecture());
// 		    TmpMatrix[l][k] = OperationUpUp.GetScalar() * NormalizationFactor;
// 		  }
// 	      File << i << " " << j << " 1 1";
// 	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 	      File << endl;
// 	      delete OperatorUpUp;
// 	    }
// 	}
//     }
//   else
//     {
//       File << "# <Psi_L| c^+_{i,sigma} c^+_{j,sigma'} +/- c^+_{i,sigma} c^+_{j,sigma'}|Psi_R> with sigma,sigma' = 0 (down) or 1 (up)" << endl
// 	   << "# for each case is given the (norm)^2, the norm and the argument" << endl;
//       File << "# i j ";
//       if ((NbrLeftStates == 1) && (NbrRightStates == 1))
// 	{
// 	  File << "<Psi_L| (c^+_{i,up} c^+_{j,up} + c^+_{i,down} c^+_{j,down} |Psi_R> <Psi_L| (c^+_{i,up} c^+_{j,up} - c^+_{i,down} c^+_{j,down} |Psi_R> <Psi_L| (c^+_{i,up} c^+_{j,down} + c^+_{i,up} c^+_{j,down} |Psi_R> <Psi_L| (c^+_{i,down} c^+_{j,up} - c^+_{i,down} c^+_{j,up} |Psi_R>" << endl;
// 	}
//       else
// 	{
// 	  File << "|<Psi_L| (c^+_{i,up} c^+_{j,up} + c^+_{i,down} c^+_{j,down} |Psi_R> |^2 |<Psi_L| (c^+_{i,up} c^+_{j,up} - c^+_{i,down} c^+_{j,down} |Psi_R> |^2 |<Psi_L| (c^+_{i,up} c^+_{j,down} + c^+_{i,down} c^+_{j,up} |Psi_R> |^2 |<Psi_L| (c^+_{i,up} c^+_{j,down} - c^+_{i,down} c^+_{j,up} |Psi_R> |^2" << endl;
// 	}
//       double NormalizationFactor = 1.0 / sqrt((double) (NbrLeftStates * NbrRightStates));
//       for (int i = 0; i < RightNbrSites; ++i)
// 	{
// 	  int Index = i;
// 	  if (RightStatistics == true)
// 	    {
// 	      if ((RightGutzwillerFlag == false) && (LeftGutzwillerFlag == false))
// 		{
// 		  ComplexMatrix TmpMatrix (NbrLeftStates, NbrRightStates);
// 		  File << i << " " << i << " (0,0) 0 0 (0,0) 0 0";
// 
// 		  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownUpDiag = 0;
// 		  if (RightMomentumFlag == false)
// 		    OperatorDownUpDiag = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 0, i, 1, 1, 0, 1.0);
// 		  else
// 		    OperatorDownUpDiag = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 0, i, 1, 1, 0, 1.0);
// 		  for (int k = 0; k < NbrLeftStates; ++k)
// 		    for (int l = 0; l < NbrRightStates; ++l)
// 		      {
// 			OperatorMatrixElementOperation OperationDownUp(OperatorDownUpDiag, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 			OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
// 			TmpMatrix[l][k] = OperationDownUp.GetScalar() * NormalizationFactor;
// 		      }
// 		  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 		  delete OperatorDownUpDiag;
// 
// 		  ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpDownDiag = 0;
// 		  if (RightMomentumFlag == false)
// 		    OperatorUpDownDiag = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, i, 0, 1, 0, 1.0);
// 		  else
// 		    OperatorUpDownDiag = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, i, 0, 1, 0, 1.0);
// 		  for (int k = 0; k < NbrLeftStates; ++k)
// 		    for (int l = 0; l < NbrRightStates; ++l)
// 		      {
// 			OperatorMatrixElementOperation OperationUpDown (OperatorUpDownDiag, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 			OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
// 			TmpMatrix[l][k] = OperationUpDown.GetScalar() * NormalizationFactor;
// 		      }
// 		  HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 		  delete OperatorUpDownDiag;
// 		  File << endl;
// 		}
// 	      ++Index;
// 	    }
// 	  for (int j = Index; j < RightNbrSites; ++j)
// 	    {
// 	      ComplexMatrix TmpMatrix (NbrLeftStates, NbrRightStates);
// 	      File << i << " " << j << " ";
// 
// 	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpUp = 0;
// 	      if (RightMomentumFlag == false)
// 		OperatorUpUp = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 1, 0, 0, 1.0);
// 	      else
// 		OperatorUpUp = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 1, 0, 0, 1.0);
// 	      for (int k = 0; k < NbrLeftStates; ++k)
// 		for (int l = 0; l < NbrRightStates; ++l)
// 		  {
// 		    OperatorMatrixElementOperation OperationUpUp(OperatorUpUp, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 		    OperationUpUp.ApplyOperation(Architecture.GetArchitecture());
// 		    TmpMatrix[l][k] = OperationUpUp.GetScalar() * NormalizationFactor;
// 		  }
// 	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 	      delete OperatorUpUp;
// 
// 	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownDown = 0;
// 	      if (RightMomentumFlag == false)
// 		OperatorDownDown = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 1, 0, 0, -1.0);
// 	      else
// 		OperatorDownDown = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 1, 0, 0, -1.0);
// 	      for (int k = 0; k < NbrLeftStates; ++k)
// 		for (int l = 0; l < NbrRightStates; ++l)
// 		  {
// 		    OperatorMatrixElementOperation OperationDownDown(OperatorDownDown, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 		    OperationDownDown.ApplyOperation(Architecture.GetArchitecture());
// 		    TmpMatrix[l][k] = OperationDownDown.GetScalar() * NormalizationFactor;
// 		  }
// 	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 	      delete OperatorDownDown;
// 
// 	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorDownUp = 0;
// 	      if (RightMomentumFlag == false)
// 		OperatorDownUp = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 0, 0, 1, 1.0);
// 	      else
// 		OperatorDownUp = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 0, 0, 1, 1.0);
// 	      for (int k = 0; k < NbrLeftStates; ++k)
// 		for (int l = 0; l < NbrRightStates; ++l)
// 		  {
// 		    OperatorMatrixElementOperation OperationDownUp(OperatorDownUp, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 		    OperationDownUp.ApplyOperation(Architecture.GetArchitecture());
// 		    TmpMatrix[l][k] = OperationDownUp.GetScalar() * NormalizationFactor;
// 		  }
// 	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 	      delete OperatorDownUp;
// 
// 	      ParticleOnSphereWithSpinSuperconductorOrderParameterOperator* OperatorUpDown = 0;
// 	      if (RightMomentumFlag == false)
// 		OperatorUpDown = new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(RightSpace, i, 1, j, 0, 0, 1, -1.0);
// 	      else
// 		OperatorUpDown = new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) RightSpace, i, 1, j, 0, 0, 1, -1.0);
// 	      for (int k = 0; k < NbrLeftStates; ++k)
// 		for (int l = 0; l < NbrRightStates; ++l)
// 		  {
// 		    OperatorMatrixElementOperation OperationUpDown(OperatorUpDown, LeftStates[k], RightStates[l], RightStates[0].GetVectorDimension());
// 		    OperationUpDown.ApplyOperation(Architecture.GetArchitecture());
// 		    TmpMatrix[l][k] = OperationUpDown.GetScalar() * NormalizationFactor;
// 		  }
// 	      HubbardSuperconductorOrderParameterMatrixDiagonalize(TmpMatrix, File);
// 	      delete OperatorUpDown;
// 	      File << endl;
// 	    }
// 	}
//     }
// 
//   File.close();
//   if ((RightMomentumFlag == true) && (LeftMomentumFlag == true))
//     {
//       FileFourierTransform.close();
//     }
  return 0;
}


// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

int GetLinearizedIndex(int xPosition, int yPosition, int atomicIndex , int nbrSpinX, int nbrSpinY)
{
  if (xPosition < 0)
    xPosition += nbrSpinX;
  if (xPosition >= nbrSpinX)
    xPosition -= nbrSpinX;
  if (yPosition < 0)
    yPosition += nbrSpinY;
  if (yPosition >= nbrSpinY)
    yPosition -= nbrSpinY;
  return (3 * ((xPosition * nbrSpinY) + yPosition) + atomicIndex);
}