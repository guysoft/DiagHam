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

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"


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
  OptionManager Manager ("HubbardEntanglementEntropyRealSpacePartition" , "0.01");
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

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");  
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardEntanglementEntropyRealSpacePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSpaces = 1;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  ComplexVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  int NbrParticles = 0;
  int NbrSites = 0;
  bool Statistics = true;
  bool GutzwillerFlag = false;
  double* Coefficients = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type HubbardEntanglementEntropyRealSpacePartition -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) && 
      (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  if ((Manager.GetString("degenerated-groundstate") != 0) && 
      (IsFile(Manager.GetString("degenerated-groundstate")) == false))
    {
      cout << "can't open file " << Manager.GetString("degenerated-groundstate") << endl;
      return -1;
    }

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      Coefficients = new double[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));
      Coefficients[0] = 1.0;
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
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));		   
	 }
       if (DegeneratedFile.GetNbrColumns() == 1)
	 {
	   Coefficients = new double[NbrSpaces];
	   for (int i = 0; i < NbrSpaces; ++i)
	     Coefficients[i] = 1.0 / ((double) NbrSpaces);
	 }
       else
	 {
	   double TmpSum = 0.0;
	   Coefficients = DegeneratedFile.GetAsDoubleArray(1);
	   for (int i = 0; i < NbrSpaces; ++i)
	     TmpSum += Coefficients[i];
	   TmpSum = 1.0 / TmpSum;
	   for (int i = 0; i < NbrSpaces; ++i)
	     Coefficients[i] *= TmpSum;
	 }
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (FTIHubbardModelFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrParticles, NbrSites, Statistics, GutzwillerFlag) == false)
 	{
 	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
 	  return -1;
 	}
    }

  GroundStates = new ComplexVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i] << endl;
	  return -1;      
	}
    }

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  N    lambda";
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  ParticleOnSphere* Space = 0;
  if (Statistics == true)
    {
      if (GutzwillerFlag == false)
	Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
      else
	Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
    }
  else
    {
      cout << "not available for bosons" << endl;
      return -1;
    }
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "ent");
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
  
  int NbrKeptOrtbitals = 2;
  int* KeptOrbitals = new int [NbrKeptOrtbitals];
  KeptOrbitals[0] = 0;
  KeptOrbitals[2] = 0;
  int MaxSubsystemNbrParticles = 2 * NbrKeptOrtbitals;
  if (GutzwillerFlag == true)
    {
      MaxSubsystemNbrParticles = NbrKeptOrtbitals;
    }
  int SubsystemNbrParticles = 0;
  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << endl;      
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      if (ShowTimeFlag == true)
	{
	  gettimeofday (&(TotalStartingTime), 0);
	}
      ComplexMatrix PartialEntanglementMatrix;
      PartialEntanglementMatrix = ((FermionOnLatticeWithSpinRealSpace*) Space)->EvaluatePartialEntanglementMatrix(SubsystemNbrParticles, NbrKeptOrtbitals, KeptOrbitals, GroundStates[0], Architecture.GetArchitecture());
      PartialEntanglementMatrix *= Coefficients[0];
      for (int i = 1; i < NbrSpaces; ++i)
	{
	  ComplexMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Space)->EvaluatePartialEntanglementMatrix(SubsystemNbrParticles, NbrKeptOrtbitals, KeptOrbitals, GroundStates[i], Architecture.GetArchitecture());
	  TmpMatrix *= Coefficients[i];
	  PartialEntanglementMatrix += TmpMatrix;
	}
      if (ShowTimeFlag == true)
	{
	  gettimeofday (&(TotalEndingTime), 0);
	  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
	}
      if (PartialEntanglementMatrix.GetNbrRow() > 1)
	{
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalStartingTime), 0);
	    }
	  double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
	  int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
	  if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
	    {
	      TmpDimension = PartialEntanglementMatrix.GetNbrRow();
	    }
	  for (int i = 0; i < TmpDimension; ++i)
	    TmpValues[i] *= TmpValues[i];
	  RealDiagonalMatrix TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
	  TmpDiag.SortMatrixDownOrder();
	  if (DensityMatrixFileName != 0)
	    {
	      ofstream DensityMatrixFile;
	      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
	      DensityMatrixFile.precision(14);
	      for (int i = 0; i < TmpDimension; ++i)
		DensityMatrixFile << SubsystemNbrParticles << " " << TmpDiag[i] << endl;
	      DensityMatrixFile.close();
	    }
	  for (int i = 0; i < TmpDimension; ++i)
	    {
	      if (TmpDiag[i] > 1e-14)
		{
		  EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
		  DensitySum +=TmpDiag[i];
		}
	    }
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "diagonalization done in " << Dt << "s" << endl;
	    }
	}
          File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
    }
  File.close();

  return 0;
}
