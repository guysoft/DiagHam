#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"

#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslations.h"

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

#include "MathTools/BinomialCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

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
  OptionManager Manager ("PEPSEntanglementSpectrum" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-sites", "number of sites", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "spin", "spin quantum number", 0);
  (*SystemGroup) += new SingleIntegerOption  ('k', "momentum", "momentum of the state", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new BooleanOption  ('c', "complex", "states of the density matrix are complex");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PEPSEntanglementSpectrum -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type PEPSEntanglementSpectrum -h" << endl;
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


  int NbrSites = Manager.GetInteger("nbr-sites"); 
  int K = Manager.GetInteger("momentum"); 
  int Sz = Manager.GetInteger("spin"); 
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");

  int NbrSpaces = 1;
  RealVector* GroundStates = 0;
  ComplexVector * ComplexGroundStates = 0;
  char** GroundStateFiles = 0;
  double EigenvalueError = 1e-10;
  DoubledSpin0_1_2_ChainWithTranslations Space (NbrSites,Sz,10000,10000);
  bool ComplexFlag = Manager.GetBoolean("complex");

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
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
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
    }
  
  
  if (ComplexFlag == false)
    {
      GroundStates = new RealVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }
	  
	  if (Space.GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and ground state" << endl;
	      return 0;
	    }
	}
    }
  else
    {
      ComplexGroundStates = new ComplexVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	{	
	  if (ComplexGroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }
	  
	  if (Space.GetLargeHilbertSpaceDimension() != ComplexGroundStates[i].GetLargeVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and ground state" << endl;
	      return 0;
	    }
	}
    }
  
  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  N    Sz    lambda";
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
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

  int MaxSubsystemSz = Sz;
  
  
  for (int SubSystemSz; SubSystemSz <=  MaxSubsystemSz; ++ SubSystemSz)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      
      cout << "processing subsystem Sz =" << SubSystemSz << endl;
      if(ComplexFlag == false)
	{
	  RealSymmetricMatrix PartialDensityMatrix = Space.EvaluatePartialDensityMatrix(SubSystemSz,GroundStates[0]);
	  for (int i = 1; i < NbrSpaces; ++i)
	    {
	      RealSymmetricMatrix TmpMatrix = Space.EvaluatePartialDensityMatrix(SubSystemSz,GroundStates[i]);
	      PartialDensityMatrix += TmpMatrix;
	    }
	  if (NbrSpaces > 1)
	    PartialDensityMatrix /= ((double) NbrSpaces);

	  if (PartialDensityMatrix.GetNbrRow() > 1)
		{
		  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
		  if (LapackFlag == true)
		    {
		      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		    }
		  else
		    {
		      PartialDensityMatrix.Diagonalize(TmpDiag);
		    }
#else
		  PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
		  TmpDiag.SortMatrixDownOrder();
		  if (DensityMatrixFileName != 0)
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			DensityMatrixFile << SubSystemSz << " " << TmpDiag[i] << endl;
		      DensityMatrixFile.close();
		    }
		  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		    {
		      if (TmpDiag[i] > 1e-14)
			{
			  EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
			  DensitySum +=TmpDiag[i];
			}
		    }
		}
	      else
		{
		  if (PartialDensityMatrix.GetNbrRow() == 1)
		    {
		      double TmpValue = PartialDensityMatrix(0,0);
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  DensityMatrixFile << SubSystemSz << " " << TmpValue << endl;
			  DensityMatrixFile.close();
			}		  
		      if (TmpValue > 1e-14)
			{
			  EntanglementEntropy += TmpValue * log(TmpValue);
			  DensitySum += TmpValue;
			}
		    }
		}
	}
      else
	{
	  HermitianMatrix PartialDensityMatrix = Space.EvaluatePartialDensityMatrix(SubSystemSz, ComplexGroundStates[0]);
	      for (int i = 1; i < NbrSpaces; ++i)
		{
		  HermitianMatrix TmpMatrix =  Space.EvaluatePartialDensityMatrix(SubSystemSz, ComplexGroundStates[i]);
		  PartialDensityMatrix += TmpMatrix;
		}
	      if (NbrSpaces > 1)
		PartialDensityMatrix /= ((double) NbrSpaces);
	      if (PartialDensityMatrix.GetNbrRow() > 1)
		{
		  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
		  if (LapackFlag == true)
		    {
		      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		    }
		  else
		    {
		      PartialDensityMatrix.Diagonalize(TmpDiag);
		    }
#else
		  PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
		  TmpDiag.SortMatrixDownOrder();
		  if (DensityMatrixFileName != 0)
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			DensityMatrixFile <<  SubSystemSz << " " << TmpDiag[i] << endl;
		      DensityMatrixFile.close();
		    }
		  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		    {
		      if (TmpDiag[i] > 1e-14)
			{
			  EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
			  DensitySum +=TmpDiag[i];
			}
		    }
		}
	      else
		{
		  if (PartialDensityMatrix.GetNbrRow() == 1)
		    {
		      double TmpValue;
		      PartialDensityMatrix.GetMatrixElement(0,0,TmpValue);
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  DensityMatrixFile << SubSystemSz << " " << TmpValue << endl;
			  DensityMatrixFile.close();
			}		  
		      if (TmpValue > 1e-14)
			{
			  EntanglementEntropy += TmpValue * log(TmpValue);
			  DensitySum += TmpValue;
			}
		    }
		}
	}
      File << SubSystemSz << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
    }
  File.close();
}
