#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Hamiltonian/SpinChainHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("SpinChainEntanglementEntropy" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-la", "minimum size of the subsystem whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-la", "maximum size of the subsystem whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the reduced density matrices in the a given file");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainEntanglementEntropy -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type SpinChainEntanglementEntropy -h" << endl;
      return -1;
    }

  int SpinValue = 0;
  int NbrSpins = 0;
  int SzValue = 0;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  bool SVDFlag = false;//Manager.GetBoolean("use-svd");
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  
  RealVector* GroundStates = 0;
  double* Weights =0;
  bool WeightFlag = false;
  int NbrSpaces = 1;
  char** GroundStateFiles = 0;
  AbstractSpinChain** Spaces = 0;
  int* TotalSz = 0;
  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalSz = new int[1];
      Weights = new double[1];
      Weights[0] = 1.0;
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));      
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       NbrSpaces = DegeneratedFile.GetNbrLines();
       GroundStateFiles = new char* [NbrSpaces];
       TotalSz = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
       if (DegeneratedFile.GetNbrColumns() > 1)
	 {
	   Weights = DegeneratedFile.GetAsDoubleArray(1);
	   WeightFlag = true;
	 }
       else
	 {
	   Weights = new double[NbrSpaces];
	   for (int i = 0; i < NbrSpaces; ++i)
	     Weights[i] = 1.0;
	 }
    }
 
  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalSz[i] = 0;
      if (SpinFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrSpins, TotalSz[i], SpinValue) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
    }

  GroundStates = new RealVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i] << endl;
	  return -1;      
	}
    }

  Spaces = new AbstractSpinChain* [NbrSpaces];

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "# l_a    Sz    lambda" << endl;
      DensityMatrixFile.close();
    }

  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName;
      TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "ent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);

  for (int i = 0; i < NbrSpaces; ++i)
    {
      switch (SpinValue)
	{
	case 1 :
	  Spaces[i] = new Spin1_2Chain (NbrSpins, TotalSz[i], 1000000);
	  break;
	case 2 :
	  Spaces[i] = new Spin1Chain (NbrSpins, TotalSz[i], 1000000);
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
  
  int SubsystemSize = Manager.GetInteger("min-la");
  if (SubsystemSize < 1)
    SubsystemSize = 1;
  int MeanSubsystemSize = NbrSpins >> 1;
  if (Manager.GetInteger("max-la") > 0)
    {
      MeanSubsystemSize = Manager.GetInteger("max-la");
      if (MeanSubsystemSize > NbrSpins)
	MeanSubsystemSize = NbrSpins;
    }
  for (; SubsystemSize <= MeanSubsystemSize; ++SubsystemSize)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      int MaxSzA = (SubsystemSize * SpinValue);
      int MinSzA = -MaxSzA;
      int MaxSzB = ((NbrSpins - SubsystemSize) * SpinValue);
      int MinSzB = -MaxSzB;
      for (; MinSzA <= MaxSzA; MinSzA += 2)
	{
	  RealSymmetricMatrix PartialDensityMatrix;
	  RealMatrix PartialEntanglementMatrix;
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      int SzB = SzValue - MinSzA;
	      if ((SzB <= MaxSzB) && (SzB >= MinSzB))
		{
		  cout << "processing subsytem size " << SubsystemSize << " SzA=" << MinSzA << endl;
		  RealSymmetricMatrix TmpPartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrix(SubsystemSize, MinSzA, GroundStates[0]);;
		  if (WeightFlag == true)
		    TmpPartialDensityMatrix *= Weights[i];
		  if (PartialDensityMatrix.GetNbrRow() == 0)
		    PartialDensityMatrix = TmpPartialDensityMatrix;
		  else
		    PartialDensityMatrix += TmpPartialDensityMatrix;
		}
	    }
	  if(SVDFlag == false)
	    {
	      if ((NbrSpaces > 1) && (WeightFlag == false))
		PartialDensityMatrix /= ((double) NbrSpaces);
	    }
	  else
	    {
	      if ((NbrSpaces > 1) && (WeightFlag == false))
		PartialEntanglementMatrix /= sqrt(((double) NbrSpaces));
	    }
	  if ((PartialDensityMatrix.GetNbrRow() > 1) || (PartialEntanglementMatrix.GetNbrRow() >= 1))
	    {
	      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
	      if (SVDFlag == false)
		{
#ifdef __LAPACK__
		  if (LapackFlag == true)
		    {
		      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		      TmpDiag.SortMatrixDownOrder();
		    }
		  else
		    {
		      PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
		      TmpDiag.SortMatrixDownOrder();
		    }
#else
		  PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
		  TmpDiag.SortMatrixDownOrder();
#endif		  				      
		  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
		    {
		      if (TmpDiag[i] > 1e-14)
			{
			  EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
			  DensitySum += TmpDiag[i];
			}
		    }
		  if (DensityMatrixFileName != 0)
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      for (int i = 0; i <TmpDiag.GetNbrRow(); ++i)
			DensityMatrixFile << SubsystemSize << " " << MinSzA << " " << TmpDiag[i] << endl;
		      DensityMatrixFile.close();
		    }
		}
	      else
		{
		  cout << "error : SVD not implemented for spin chains" << endl; 
		}
	    }
	  else
	    {
	      if (PartialDensityMatrix.GetNbrRow() == 1)
		{
		  double TmpValue = PartialDensityMatrix(0,0);
		  if (TmpValue > 1e-14)
		    {
		      EntanglementEntropy += TmpValue * log(TmpValue);
		      DensitySum += TmpValue;
		    }
		  if (DensityMatrixFileName != 0)
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      DensityMatrixFile << SubsystemSize << " " << MinSzA << " " << TmpValue << endl;
		      DensityMatrixFile.close();
		    }		  
		}
	    }
	}
      File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << endl;
    }
  File.close();
  return 0;
}
