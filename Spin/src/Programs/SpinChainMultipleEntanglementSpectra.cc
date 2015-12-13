#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Hamiltonian/SpinChainHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "Options/Options.h"

#include "GeneralTools/Endian.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(4); 

  // some running options and help
  OptionManager Manager ("SpinChainMultipleEntanglementSpectra" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 1);
  (*SystemGroup) += new SingleStringOption  ('\n', "multiple-states", "provide as a matrix a series of states whose entanglement spectrum/entropy have to be computed");  
  (*SystemGroup) += new SingleIntegerOption ('\n', "min-multiplestates", "index of the first state to consider in the matrix provided by the --multiple-groundstate option", 0);  
  (*SystemGroup) += new SingleIntegerOption ('\n', "max-multiplestates", "index of the last state to consider in the matrix provided by the --multiple-groundstate option (negative if this is the last available state)", -1);  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "consider complex wave function");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "la", "subsystem size (negative if half of the system has to be considered)", -1);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainMultipleEntanglementSpectra -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("multiple-states") == 0)
    {
      cout << "error, an eigenstate file should be provided. See man page for option syntax or type SpinChainMultipleEntanglementSpectra -h" << endl;
      return -1;
    }

  int SpinValue = 0;
  int NbrSpins = 0;
  int TotalSz = 0;
  int SubsystemSize = Manager.GetInteger("la");
  bool SzFlag = true;
  if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, TotalSz, SpinValue) == false)
    {
      SzFlag = false;
      if (SpinFindSystemInfoFromFileName(Manager.GetString("multiple-states"), NbrSpins, SpinValue) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("multiple-states") << endl;
	  return -1;
	}
    }
  if (SzFlag == true)
    cout << "N=" << NbrSpins << " Sz=" <<  TotalSz << " 2s=" << SpinValue << endl;
  else
    cout << "N=" << NbrSpins << " 2s=" << SpinValue << endl;

  if (SubsystemSize < 0)
    {
      SubsystemSize = NbrSpins / 2;
    }
  char* TmpExtenstion = new char[32];
  sprintf(TmpExtenstion, "_la_%d.full.ent", SubsystemSize);
  char* DensityMatrixFileName = ReplaceExtensionToFileName(Manager.GetString("multiple-states"), ".mat", TmpExtenstion);

  int MinIndex = Manager.GetInteger("min-multiplestates");
  int MaxIndex = Manager.GetInteger("max-multiplestates");
  RealMatrix RealEigenstates;
  ComplexMatrix ComplexEigenstates;

  if (Manager.GetBoolean("complex") == false)
    {
      if (RealEigenstates.ReadMatrix(Manager.GetString("multiple-states")) == false)
	{
	  cout << "can't read " << Manager.GetString("multiple-states") << endl;
	}
      if (MaxIndex < 0)
	MaxIndex = RealEigenstates.GetNbrColumn() - 1;
    }
  else
    {
      if (ComplexEigenstates.ReadMatrix(Manager.GetString("multiple-states")) == false)
	{
	  cout << "can't read " << Manager.GetString("multiple-states") << endl;
	}
      if (MaxIndex < 0)
	MaxIndex = ComplexEigenstates.GetNbrColumn() - 1;
    }
  int NbrStates = (MaxIndex - MinIndex + 1);

  AbstractSpinChain* Space;

  if (SzFlag == true)
    {
      switch (SpinValue)
	{
	case 1 :
	  Space = new Spin1_2Chain (NbrSpins, TotalSz, 1000000);
	  break;
	case 2 :
	  Space = new Spin1Chain (NbrSpins, TotalSz, 1000000);
	  break;
	default :
	  {
	    if ((SpinValue & 1) == 0)
	      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
	    else 
	  cout << "spin " << SpinValue << "/2 are not available" << endl;
	    return -1;
	  }
	}
    }
  else
    {
      switch (SpinValue)
	{
	case 1 :
	  Space = new Spin1_2ChainFull (NbrSpins);
	  break;
	default :
	  {
	    if ((SpinValue & 1) == 0)
	      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
	    else 
	  cout << "spin " << SpinValue << "/2 are not available" << endl;
	    return -1;
	  }
	}
    }
 
  int MaxSzA = (SubsystemSize * SpinValue);
  int MinSzA = -MaxSzA;
  int MaxSzB = ((NbrSpins - SubsystemSize) * SpinValue);
  int MinSzB = -MaxSzB;
  if (SzFlag == false)
    {
      MaxSzA = 0;
      MinSzA = 0;
      MaxSzB = 0;
      MinSzB =0;
      TotalSz = 0;
    }
  double*** EntanglementSpectra = new double**[((MaxSzA - MinSzA) / 2) + 1];
  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
    {
      EntanglementSpectra[(TmpSzA - MinSzA) / 2] = new double*[NbrStates];
    }
  int* EntanglementSpectrumDimension = new int[((MaxSzA - MinSzA) / 2) + 1];
  
  if (Manager.GetBoolean("complex") == false)
    {
      for (int TmpSzA = MinSzA; TmpSzA <= MaxSzA; TmpSzA += 2)
	{
	  int SzB = TotalSz - TmpSzA;
	  if ((SzB <= MaxSzB) && (SzB >= MinSzB))
	    {
	      if (SzFlag == true)
		{
		  cout << "processing subsytem size " << SubsystemSize << " SzA=" << TmpSzA << endl;
		}
	      else
		{
		  cout << "processing subsytem size " << SubsystemSize << endl;
		}
	      for (int Index = MinIndex; Index <= MaxIndex; ++Index)
		{
		  RealMatrix PartialEntanglementMatrix = Space->EvaluatePartialEntanglementMatrix(SubsystemSize, TmpSzA, RealEigenstates[Index]);
		  if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
		    {
		      double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
		      int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
		      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
			{
			  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
			}
		      for (int i = 0; i < TmpDimension; ++i)
			TmpValues[i] *= TmpValues[i];   
		      SortArrayDownOrdering(TmpValues, TmpDimension);
		      EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] = TmpDimension;
		      EntanglementSpectra[(TmpSzA - MinSzA) / 2][Index - MinIndex] = TmpValues;
		    }
		  else
		    {
		      double* TmpValues = new double[1];
		      TmpValues[0] = 0.0;
		      if (PartialEntanglementMatrix.GetNbrRow() == 1)
			{
			  for (int i = 0; i < PartialEntanglementMatrix.GetNbrColumn(); ++i)
			    TmpValues[0] += PartialEntanglementMatrix[i][0] * PartialEntanglementMatrix[i][0];
			}
		      else
			{
			  for (int i = 0; i < PartialEntanglementMatrix.GetNbrRow(); ++i)
			    TmpValues[0] += PartialEntanglementMatrix[0][i] * PartialEntanglementMatrix[0][i];				  
			}
		      EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] = 1;
		      EntanglementSpectra[(TmpSzA - MinSzA) / 2][Index - MinIndex] = TmpValues;
		    }
		}
	    }
	  else
	    {
	      EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] = 0;
	    }  
	}
    }
  else
    {
      for (int TmpSzA = MinSzA; TmpSzA <= MaxSzA; TmpSzA += 2)
	{
	  int SzB = TotalSz - TmpSzA;
	  if ((SzB <= MaxSzB) && (SzB >= MinSzB))
	    {
	      if (SzFlag == true)
		{
		  cout << "processing subsytem size " << SubsystemSize << " SzA=" << TmpSzA << endl;
		}
	      else
		{
		  cout << "processing subsytem size " << SubsystemSize << endl;
		}
	      for (int Index = MinIndex; Index <= MaxIndex; ++Index)
		{
		  ComplexMatrix PartialEntanglementMatrix = Space->EvaluatePartialEntanglementMatrix(SubsystemSize, TmpSzA, ComplexEigenstates[Index]);
		  if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
		    {
		      double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
		      int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
		      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
			{
			  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
			}
		      for (int i = 0; i < TmpDimension; ++i)
			TmpValues[i] *= TmpValues[i];   
		      SortArrayDownOrdering(TmpValues, TmpDimension);
		      EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] = TmpDimension;
		      EntanglementSpectra[(TmpSzA - MinSzA) / 2][Index - MinIndex] = TmpValues;
		    }
		  else
		    {
		      double* TmpValues = new double[1];
		      TmpValues[0] = 0.0;
		      if (PartialEntanglementMatrix.GetNbrRow() == 1)
			{
			  for (int i = 0; i < PartialEntanglementMatrix.GetNbrColumn(); ++i)
			    TmpValues[0] += SqrNorm(PartialEntanglementMatrix[i][0]);
			}
		      else
			{
			  for (int i = 0; i < PartialEntanglementMatrix.GetNbrRow(); ++i)
			    TmpValues[0] += SqrNorm(PartialEntanglementMatrix[0][i]);				  
			}
		      EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] = 1;
		      EntanglementSpectra[(TmpSzA - MinSzA) / 2][Index - MinIndex] = TmpValues;
		    }
		}
	    }
	  else
	    {
	      EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] = 0;
	    }  
	}
    }
  
  int TotalDimensionPerState = 0;
  int NbrSzASectors = 0;
  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
    {
      TotalDimensionPerState += EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2];
      if (EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] > 0)
	{
	  NbrSzASectors++;
	}
    }
  ofstream File;
  File.open(DensityMatrixFileName, ios::binary | ios::out);
  WriteLittleEndian(File, NbrStates);  
  WriteLittleEndian(File, NbrSzASectors);  
  for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
    {
      if (EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] > 0)
	{
	  WriteLittleEndian(File, TmpSzA);  
	}
    }
  
  for (int Index = MinIndex; Index <= MaxIndex; ++Index)
    {
      for (int TmpSzA = MinSzA;  TmpSzA <= MaxSzA; TmpSzA += 2)
	{
	  if (EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2] > 0)
	    {
	      WriteLittleEndian(File, EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2]);
	      WriteBlockLittleEndian(File, EntanglementSpectra[(TmpSzA - MinSzA) / 2][Index - MinIndex], EntanglementSpectrumDimension[(TmpSzA - MinSzA) / 2]);
	    }
	}
    }
  File.close();
  return 0;
}
