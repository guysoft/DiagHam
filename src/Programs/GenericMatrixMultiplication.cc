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
#include "GeneralTools/MultiColumnASCIIFile.h"

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


// load a matrix from a given line of the matrix product description file
//
// productFile = reference on the matrix product description file
// index = index of the line to parse in productFile
// return value = corresponding complex matrix
RealMatrix GetRealMatrixFromProductFile(MultiColumnASCIIFile& productFile, int index);

// load a matrix from a given line of the matrix product description file
//
// productFile = reference on the matrix product description file
// index = index of the line to parse in productFile
// return value = corresponding complex matrix
ComplexMatrix GetComplexMatrixFromProductFile(MultiColumnASCIIFile& productFile, int index);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("GenericMatrixMultiplication" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup  = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleStringOption ('\0', "product", "name of the file that the matrix product description");
  (*SystemGroup) += new  SingleStringOption ('\n', "export-transformation", "optional file name to export the transformation matrix that convert th original basis into the new one");
  (*SystemGroup) += new  SingleStringOption ('\n', "export-bintransformation", "optional file name to export the transformation matrix that convert th original basis into the new one");
#ifdef __LAPACK__
  (*OutputGroup) += new  SingleStringOption ('\n', "export-matrix", "export the result as a binary matrix");
  (*OutputGroup) += new  SingleStringOption ('\n', "export-asciimatrix", "export the result as an ascii matrix");
  (*OutputGroup) += new  SingleStringOption ('\n', "export-vectors", "export the result as a series of column vectors, with a given file prefix");
  (*OutputGroup) += new BooleanOption  ('s', "std-output", "use standard output instead of an output file");
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericMatrixMultiplication -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("product") == 0)
    {
      cout << "a matrix product file description has to be provided" << endl;
      return -1;
    }
  
  MultiColumnASCIIFile MatrixProductFile;
  if (MatrixProductFile.Parse(Manager.GetString("product")) == false)
    {
      MatrixProductFile.DumpErrors(cout);
      return -1;
    } 
  int NbrMatrices = MatrixProductFile.GetNbrLines();
  bool ComplexFlag = false;
  for (int i = 0; i < NbrMatrices; ++i)
    {
      if ((strcmp(MatrixProductFile(1, i), "ComplexMatrix") == 0) || (strcmp(MatrixProductFile(1, i), "ComplexVector") == 0))
	{
	  ComplexFlag = true;
	}
    }
  
  if (ComplexFlag == true)
    {
      ComplexMatrix TmpMatrix1 = GetComplexMatrixFromProductFile(MatrixProductFile, 0);
      if ((TmpMatrix1.GetNbrRow() == 0) && (TmpMatrix1.GetNbrColumn() == 0))
	{
	  return -1;
	}
      for (int i = 1; i < NbrMatrices; ++i)
	{
	  ComplexMatrix TmpMatrix2 = GetComplexMatrixFromProductFile(MatrixProductFile, i);
	  if ((TmpMatrix2.GetNbrRow() == 0) && (TmpMatrix2.GetNbrColumn() == 0))
	    {
	      return -1;
	    }
	  if (TmpMatrix2.GetNbrRow() != TmpMatrix1.GetNbrColumn())
	    {
	      cout << "error, the number of column in " << MatrixProductFile(0, i - 1) << " does not match the number of row in " 
		   << MatrixProductFile(0, i) << " ( " << TmpMatrix1.GetNbrColumn()<< " vs " << TmpMatrix2.GetNbrRow()<< " )" << endl;
	      return -1;
	    }
	  TmpMatrix1.Multiply(TmpMatrix2);
	}
      if (Manager.GetBoolean("std-output") == true)
	{
	  cout << TmpMatrix1 << endl;
	}
      else
	{
	  if (Manager.GetString("export-vectors") != 0)
	    {
	      for (int i = 0; i < TmpMatrix1.GetNbrColumn(); ++i)
		{
		  char* OutputFile = new char [strlen(Manager.GetString("export-vectors")) + 256];
		  sprintf (OutputFile, "%s.%d.vec", Manager.GetString("export-vectors"), i);
		  if (TmpMatrix1[i].WriteVector(OutputFile) == false)
		    {
		      cout << "error, can't write " << OutputFile << endl;
		      return -1;
		    }
		  delete[] OutputFile;
		}
	    }
	  if (Manager.GetString("export-matrix") != 0)
	    {
	      TmpMatrix1.WriteMatrix(Manager.GetString("export-matrix"));
	    }
	  if (Manager.GetString("export-asciimatrix") != 0)
	    {
	      TmpMatrix1.WriteAsciiMatrix(Manager.GetString("export-asciimatrix"));
	    }
	}
    }
  else
    {
      RealMatrix TmpMatrix1 = GetRealMatrixFromProductFile(MatrixProductFile, 0);
      if ((TmpMatrix1.GetNbrRow() == 0) && (TmpMatrix1.GetNbrColumn() == 0))
	{
	  return -1;
	}
      for (int i = 1; i < NbrMatrices; ++i)
	{
	  RealMatrix TmpMatrix2 = GetRealMatrixFromProductFile(MatrixProductFile, i);
	  if ((TmpMatrix2.GetNbrRow() == 0) && (TmpMatrix2.GetNbrColumn() == 0))
	    {
	      return -1;
	    }
	  if (TmpMatrix2.GetNbrRow() != TmpMatrix1.GetNbrColumn())
	    {
	      cout << "error, the number of column in " << MatrixProductFile(0, i - 1) << " does not match the number of row in " 
		   << MatrixProductFile(0, i) << " ( " << TmpMatrix1.GetNbrColumn()<< " vs " << TmpMatrix2.GetNbrRow()<< " )" << endl;
	      return -1;
	    }
	  TmpMatrix1.Multiply(TmpMatrix2);
	}
      if (Manager.GetBoolean("std-output") == true)
	{
	  cout << TmpMatrix1 << endl;
	}
      else
	{
	  if (Manager.GetString("export-vectors") != 0)
	    {
	      for (int i = 0; i < TmpMatrix1.GetNbrColumn(); ++i)
		{
		  char* OutputFile = new char [strlen(Manager.GetString("export-vectors")) + 256];
		  sprintf (OutputFile, "%s.%d.vec", Manager.GetString("export-vectors"), i);
		  if (TmpMatrix1[i].WriteVector(OutputFile) == false)
		{
		  cout << "error, can't write " << OutputFile << endl;
		  return -1;
		}
		  delete[] OutputFile;
		}
	    }
	  if (Manager.GetString("export-matrix") != 0)
	    {
	      TmpMatrix1.WriteMatrix(Manager.GetString("export-matrix"));
	    }
	  if (Manager.GetString("export-asciimatrix") != 0)
	    {
	      TmpMatrix1.WriteAsciiMatrix(Manager.GetString("export-asciimatrix"));
	    }
	}
    }

  return 0;
}

// load a matrix from a given line of the matrix product description file
//
// productFile = reference on the matrix product description file
// index = index of the line to parse in productFile
// return value = corresponding complex matrix

RealMatrix GetRealMatrixFromProductFile(MultiColumnASCIIFile& productFile, int index)
{
  if (strcmp(productFile(1, index), "RealMatrix") == 0)
    {
      RealMatrix TmpMatrix;
      TmpMatrix.ReadMatrix(productFile(0, index));
      return TmpMatrix;
    }
  if (strcmp(productFile(1, index), "RealVector") == 0)
    {
      MultiColumnASCIIFile VectorFile;
      if (VectorFile.Parse(productFile(0, index)) == false)
	{
	  VectorFile.DumpErrors(cout);
	  return RealMatrix();
	} 
      int NbrVectors = VectorFile.GetNbrLines();
      RealVector* TmpVectors = new RealVector[NbrVectors];
      for (int i = 0; i < NbrVectors; ++i)
	{
	  if (TmpVectors[i].ReadVector(VectorFile(0, i)) == false)
	    {
	      cout << "can't read " << VectorFile(0, i) << endl;
	      return RealMatrix();
	    }
	}
      return RealMatrix(TmpVectors, NbrVectors);
    }
  return RealMatrix();
}

// load a matrix from a given line of the matrix product description file
//
// productFile = reference on the matrix product description file
// index = index of the line to parse in productFile
// return value = corresponding complex matrix

ComplexMatrix GetComplexMatrixFromProductFile(MultiColumnASCIIFile& productFile, int index)
{
  if (strcmp(productFile(1, index), "ComplexMatrix") == 0)
    {
      ComplexMatrix TmpMatrix;
      TmpMatrix.ReadMatrix(productFile(0, index));
      return TmpMatrix;
    }
  if (strcmp(productFile(1, index), "RealMatrix") == 0)
    {
      RealMatrix TmpMatrix;
      TmpMatrix.ReadMatrix(productFile(0, index));
      return ComplexMatrix(TmpMatrix);
    }
  if (strcmp(productFile(1, index), "ComplexVector") == 0)
    {
      MultiColumnASCIIFile VectorFile;
      if (VectorFile.Parse(productFile(0, index)) == false)
	{
	  VectorFile.DumpErrors(cout);
	  return ComplexMatrix();
	} 
      int NbrVectors = VectorFile.GetNbrLines();
      ComplexVector* TmpVectors = new ComplexVector[NbrVectors];
      for (int i = 0; i < NbrVectors; ++i)
	{
	  if (TmpVectors[i].ReadVector(VectorFile(0, i)) == false)
	    {
	      cout << "can't read " << VectorFile(0, i) << endl;
	      return ComplexMatrix();
	    }
	}
      return ComplexMatrix(TmpVectors, NbrVectors);
    }
  if (strcmp(productFile(1, index), "RealVector") == 0)
    {
      MultiColumnASCIIFile VectorFile;
      if (VectorFile.Parse(productFile(0, index)) == false)
	{
	  VectorFile.DumpErrors(cout);
	  return ComplexMatrix();
	} 
      int NbrVectors = VectorFile.GetNbrLines();
      ComplexVector* TmpVectors = new ComplexVector[NbrVectors];
      for (int i = 0; i < NbrVectors; ++i)
	{
	  RealVector TmpVector;
	  if (TmpVector.ReadVector(VectorFile(0, i)) == false)
	    {
	      cout << "can't read " << VectorFile(0, i) << endl;
	      return ComplexMatrix();
	    }
	  TmpVectors[i] = TmpVector;
	}
      return ComplexMatrix(TmpVectors, NbrVectors);
    }
  return ComplexMatrix();
}
