#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
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
  OptionManager Manager ("ExtractLinearlyIndependentVectors" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  Manager += SystemGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new  SingleStringOption ('b', "basis", "name of the file that contains the vector files used to describe the basis");
  (*SystemGroup) += new BooleanOption ('\n', "check-only", "check how many vectors are linearly independent without extracting the basis");
  (*SystemGroup) += new SingleDoubleOption ('\n', "error", "bound above which vectors are consider as linearly independent", 1e-10);
  (*SystemGroup) += new  SingleStringOption ('\n', "vector-prefix", "prefix to use for each vector of the basis", "vector");
  (*SystemGroup) += new  SingleStringOption ('\n', "export-overlap", "optional file name to export the overlap matrix (in text mode)");
  (*SystemGroup) += new  SingleStringOption ('\n', "export-binoverlap", "optional file name to export the overlap matrix (in binary mode)");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type ExtractLinearlyIndependentVectors -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  char* BasisDescription = Manager.GetString("basis");  
  char* VectorPrefix = Manager.GetString("vector-prefix");
  double Error = ((SingleDoubleOption*) Manager["error"])->GetDouble();

  ConfigurationParser ReducedBasis;
  if (ReducedBasis.Parse(BasisDescription) == false)
    {
      ReducedBasis.DumpErrors(cout) << endl;
      return -1;
    }
  int NbrVectors;
  char** VectorFileNames;
  if (ReducedBasis.GetAsStringArray("Basis", ' ', VectorFileNames, NbrVectors) == false)
    {
      cout << "Vectors are not defined or have a wrong value in " << BasisDescription << endl;
      return -1;
    }

  if(Manager.GetBoolean("complex") == true)
    {
      ComplexVector * Basis = new ComplexVector[NbrVectors];
      HermitianMatrix HRep(NbrVectors);
      char* DirectoryName = ReducedBasis["Directory"];
      char* TmpName;
      for (int i = 0; i < NbrVectors; ++i)
	{
	  TmpName = VectorFileNames[i];
	  if (DirectoryName != 0)
	    {
	      TmpName = ConcatenatePathAndFileName(DirectoryName, TmpName);
	    }
	  cout << "reading vector " << TmpName << endl;
	  if (Basis[i].ReadVector(TmpName) == false)
	    {
	      cout << "error while reading " << TmpName << endl;
	      if (DirectoryName != 0)
		delete[] DirectoryName;
	      for (int j = 0; j < NbrVectors; ++j)
		delete[] VectorFileNames[j];
	      delete[] VectorFileNames;
	      return -1;
	    }
	  if (DirectoryName != 0)
	    delete[] TmpName;
	}   
      for (int i = 0; i < NbrVectors; ++i)
	for (int j = 0; j < NbrVectors; ++j)
	  HRep.SetMatrixElement(i ,j, Basis[j] * Basis[i]);
      
      if (Manager.GetString("export-overlap") != 0)
	{
	  HRep.WriteAsciiMatrix(Manager.GetString("export-overlap"));
	}
      if (Manager.GetString("export-binoverlap") != 0)
	{
	  HRep.WriteMatrix(Manager.GetString("export-binoverlap"));
	}

      RealDiagonalMatrix TmpDiag (NbrVectors);
      
      
      if (Manager.GetBoolean("check-only") == true)	
	{
#ifdef __LAPACK__
	  if (Manager.GetBoolean("use-lapack") == true)
	    HRep.LapackDiagonalize(TmpDiag);
	  else
	    HRep.Diagonalize(TmpDiag);
#else
	  HRep.Diagonalize(TmpDiag);
#endif		  
	  int Count = 0;
	  for (int i = 0; i < NbrVectors; ++i)
	    {
	      cout << TmpDiag[i] << " ";
	      if (fabs(TmpDiag[i]) > Error)
		Count++;
	    }
	  cout << endl;
	  cout << Count << " linearly independent vectors" << endl;
	}
      else
	{
	  ComplexMatrix TmpEigenvector (NbrVectors, NbrVectors, true);	      
	  for (int l = 0; l < NbrVectors; ++l)
	    TmpEigenvector(l, l) = 1.0;
#ifdef __LAPACK__
	  if (Manager.GetBoolean("use-lapack") == true)
	    HRep.LapackDiagonalize(TmpDiag, TmpEigenvector);
	  else
	    HRep.Diagonalize(TmpDiag, TmpEigenvector);
#else
	  HRep.Diagonalize(TmpDiag, TmpEigenvector);
#endif
	  char* OutputVectorFileName = new char [strlen(VectorPrefix) + 32];
	  int TmpDimension = Basis[0].GetVectorDimension();
	  int Count = 0;
	  for (int i = 0; i < NbrVectors; ++i)
	    {
	      cout << TmpDiag[i] << " ";
	      if (fabs(TmpDiag[i]) > Error)
		{
		  ComplexVector TmpVector (TmpDimension, true);
		  for (int j = 0; j < NbrVectors; ++j)
		    for (int k = 0; k < TmpDimension; ++k)
		      {
			TmpVector[k] += TmpEigenvector[i][j] * Basis[j][k];
		      }
		  TmpVector /= TmpVector.Norm();
		  sprintf (OutputVectorFileName, "%s%d.vec", VectorPrefix, Count);
		  TmpVector.WriteVector(OutputVectorFileName);
		  Count++;
		}
	    }
	  cout << endl;
	  cout << Count << " linearly independent vectors" << endl;
	  delete[] OutputVectorFileName;
	}
      
      if (DirectoryName != 0)
	delete[] DirectoryName;
      for (int j = 0; j < NbrVectors; ++j)
	delete[] VectorFileNames[j];
      delete[] VectorFileNames;
      return 0;
      
    }
  else
    {
      RealVector * Basis = new RealVector[NbrVectors];
      
      char* DirectoryName = ReducedBasis["Directory"];
      char* TmpName;
      for (int i = 0; i < NbrVectors; ++i)
	{
	  TmpName = VectorFileNames[i];
	  if (DirectoryName != 0)
	    {
	      TmpName = ConcatenatePathAndFileName(DirectoryName, TmpName);
	    }
	  cout << "reading vector " << TmpName << endl;
	  if (Basis[i].ReadVector(TmpName) == false)
	    {
	      cout << "error while reading " << TmpName << endl;
	      if (DirectoryName != 0)
		delete[] DirectoryName;
	      for (int j = 0; j < NbrVectors; ++j)
		delete[] VectorFileNames[j];
	      delete[] VectorFileNames;
	      return -1;
	    }
	  if (DirectoryName != 0)
	    delete[] TmpName;
	}
      
      RealSymmetricMatrix HRep(NbrVectors);
      for (int i = 0; i < NbrVectors; ++i)
	for (int j = 0; j < NbrVectors; ++j)
	  HRep(i ,j) = Basis[j] * Basis[i];
      if (Manager.GetString("export-overlap") != 0)
	{
	  HRep.WriteAsciiMatrix(Manager.GetString("export-overlap"));
	}
      if (Manager.GetString("export-binoverlap") != 0)
	{
	  HRep.WriteMatrix(Manager.GetString("export-binoverlap"));
	}
      RealDiagonalMatrix TmpDiag (NbrVectors);
      
      
      if (Manager.GetBoolean("check-only") == true)	
	{
#ifdef __LAPACK__
	  if (Manager.GetBoolean("use-lapack") == true)
	    HRep.LapackDiagonalize(TmpDiag);
	  else
	    HRep.Diagonalize(TmpDiag);
#else
	  HRep.Diagonalize(TmpDiag);
#endif		  
	  int Count = 0;
	  for (int i = 0; i < NbrVectors; ++i)
	    {
	      cout << TmpDiag[i] << " ";
	      if (fabs(TmpDiag[i]) > Error)
		Count++;
	    }
	  cout << endl;
	  cout << Count << " linearly independent vectors" << endl;
	}
      else
	{
	  RealMatrix TmpEigenvector (NbrVectors, NbrVectors, true);	      
	  for (int l = 0; l < NbrVectors; ++l)
	    TmpEigenvector(l, l) = 1.0;
#ifdef __LAPACK__
	  if (Manager.GetBoolean("use-lapack") == true)
	    HRep.LapackDiagonalize(TmpDiag, TmpEigenvector);
	  else
	    HRep.Diagonalize(TmpDiag, TmpEigenvector);
#else
	  HRep.Diagonalize(TmpDiag, TmpEigenvector);
#endif
	  char* OutputVectorFileName = new char [strlen(VectorPrefix) + 32];
	  int TmpDimension = Basis[0].GetVectorDimension();
	  int Count = 0;
	  for (int i = 0; i < NbrVectors; ++i)
	    {
	      cout << TmpDiag[i] << " ";
	      if (fabs(TmpDiag[i]) > Error)
		{
		  RealVector TmpVector (TmpDimension, true);
		  for (int j = 0; j < NbrVectors; ++j)
		    for (int k = 0; k < TmpDimension; ++k)
		      {
			TmpVector[k] += TmpEigenvector[i][j] * Basis[j][k];
		      }
		  TmpVector /= TmpVector.Norm();
		  sprintf (OutputVectorFileName, "%s%d.vec", VectorPrefix, Count);
		  TmpVector.WriteVector(OutputVectorFileName);
		  Count++;
		}
	    }
	  cout << endl;
	  cout << Count << " linearly independent vectors" << endl;
	  delete[] OutputVectorFileName;
	}
      
      if (DirectoryName != 0)
	delete[] DirectoryName;
      for (int j = 0; j < NbrVectors; ++j)
	delete[] VectorFileNames[j];
      delete[] VectorFileNames;
      return 0;
    }
  return 0;
}
