#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/ArrayTools.h"

#include "Operator/QHEOperator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "HilbertSpace/QHEHilbertSpace/BosonOnSphere.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


void SortOverlaps (double* bestOverlaps, double** overlaps, int nbrOverlaps, int nbrTestStates);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHENBodyQuasiHoleOverlap" , "0.01");
  OptionGroup* MainGroup = new OptionGroup ("main options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += MiscGroup;
  Manager += MainGroup;

  (*MainGroup) += new SingleStringOption  ('\n', "input-file", "name of the file which contains definition of overlaps to evaluate");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHENBodyQuasiHoleOverlap -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  ConfigurationParser OverlapDefinition;
  if (OverlapDefinition.Parse(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      OverlapDefinition.DumpErrors(cout) << endl;
      return -1;
    }

  if (OverlapDefinition["Degeneracy"] == 0)
    {
      cout << "no Degeneracy defined in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  char* TmpString = OverlapDefinition["Degeneracy"];
  int MaxNbrLz = (strlen(TmpString) >> 1) + 1;
  if (MaxNbrLz == 0)
     {
      cout << "error while parsing Degeneracy in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }

  int* Degeneracy = new int [MaxNbrLz];
  MaxNbrLz = 0;
  char* End = TmpString;
  bool ErrorFlag = false;
  while ((*TmpString) != '\0')
    {
      while (((*End) != '\0') && ((*End) != ' ') && ((*End) != '\t'))
	{
	  if ((((*End) < '0') || ((*End) > '9')))
	    ErrorFlag = true;
	  ++End;
	}
      if (ErrorFlag == true)
	{
	  cout << "error while parsing Degeneracy in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
	  return -1;     
	}
      Degeneracy[MaxNbrLz] = atoi (TmpString);
      ++MaxNbrLz;
      while ((((*End) == ' ') || ((*End) == '\t')) && ((*End) != '\0'))
	{
	  ++End;
	}
      TmpString = End;
    }

  if (OverlapDefinition["InputVectors"] == 0)
    {
      cout << "no InputVectors defined in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  if (OverlapDefinition["OutputVectors"] == 0)
    {
      cout << "no OutputVectors defined in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  
  char* InputVectors = new char [strlen(OverlapDefinition["InputVectors"]) + 24];
  char* OutputVectors = new char [strlen(OverlapDefinition["OutputVectors"]) + 24];
  strcpy (InputVectors, OverlapDefinition["InputVectors"]);
  strcpy (OutputVectors, OverlapDefinition["OutputVectors"]);
  char* InputVectors2 = InputVectors + strlen(InputVectors);
  char* OutputVectors2 = OutputVectors + strlen(OutputVectors);

  for (int i = 6; i < MaxNbrLz; ++i)
    {
      int Lz = i << 1;
      cout << "Lz = " << i << endl;
      double* BestOverlaps = new double[Degeneracy[i]];
      double** Overlaps = new double*[Degeneracy[i]];      
      RealVector* TestVectors = new RealVector[Degeneracy[i]];
 
      for (int j = 0; j < Degeneracy[i]; ++j)
	{
	  sprintf (OutputVectors2, "%d.%d.vec", Lz, j);	      
	  if (TestVectors[j].ReadVector(OutputVectors) == false)
	    {
	      cout << "error while reading " << OutputVectors << endl;
	      return -1;
	    }
     
	}
      for (int j = 0; j < Degeneracy[i]; ++j)
	{
	  Overlaps[j] = new double [Degeneracy[i]];
	  sprintf (InputVectors2, "%d.%d.vec", Lz, j);
	  cout << InputVectors << ":" << endl;
	  RealVector ReferenceVector;	  
	  if (ReferenceVector.ReadVector(InputVectors) == false)
	    {
	      cout << "error while reading " << InputVectors << endl;
	      return -1;
	    }
	  BestOverlaps[j] = 0.0;
	  
	  for (int k = 0; k < Degeneracy[i]; ++k)
	    {
	      if (ReferenceVector.GetVectorDimension() != TestVectors[k].GetVectorDimension())
		{
		  sprintf (OutputVectors2, "%d.%d.vec", Lz, k);
		  cout << "dimension mismatch between " << InputVectors2 << " and " << OutputVectors2 << endl;
		  return -1;
		}
	      
	      double Scalar = ReferenceVector * TestVectors[k];
	      Overlaps[j][k] = Scalar;
	      cout << "    " << OutputVectors << "  = " << Scalar << endl;
	      BestOverlaps[j] += Scalar * Scalar;
	    }
	  BestOverlaps[j] = sqrt(BestOverlaps[j]);
	  cout << endl << "best overlap = " << BestOverlaps[j] << endl;

	  BosonOnSphere Space (8, i, 8);
	  RealSymmetricMatrix HRep (Degeneracy[i], Degeneracy[i]);
	  ParticleOnSphereSquareTotalMomentumOperator oper(&Space, i, 8);
	  for (int k = 0; k < Degeneracy[i]; ++k)
	    {
	      for (int l = k; l < Degeneracy[i]; ++l)
		{	
		  HRep(l, k) = oper.MatrixElement(TestVectors[k], TestVectors[l]).Re;
		}
	    }
	  cout << endl << endl;
	  RealTriDiagonalSymmetricMatrix TmpTriDiag (Degeneracy[i]);
	  cout << HRep << endl;
	  HRep.Householder(TmpTriDiag, 1e-7);
	  TmpTriDiag.Diagonalize();
	  TmpTriDiag.SortMatrixUpOrder();
	  cout << "eigenvalues = " << endl;
	  for (int k = 0; k < Degeneracy[i]; ++k)	    
	    {
	      cout << TmpTriDiag.DiagonalElement(k) << endl;
	    }
	  return 0;

	}
      SortOverlaps(BestOverlaps, Overlaps, Degeneracy[i], Degeneracy[i]);
      cout << endl << "overlaps = " << endl;
      for (int j = 0; j < Degeneracy[i]; ++j)
	{
	  cout << BestOverlaps[j] << endl;
	}
      cout << "----------------------------------------------------" << endl;
      delete[] BestOverlaps;
      for (int j = 0; j < Degeneracy[i]; ++j)
	{
	  delete[] Overlaps[j];
	}
      delete[] Overlaps;
      delete[] TestVectors;
    }

  delete[] Degeneracy;
  delete[] InputVectors;
  delete[] OutputVectors;

  return 0;
}


void SortOverlaps (double* bestOverlaps, double** overlaps, int nbrOverlaps, int nbrTestStates)
{
  if (nbrOverlaps <= 1)
    return;
  SortArrayUpOrdering(bestOverlaps, overlaps, nbrOverlaps);
  for (int i = 0; i < (nbrOverlaps - 1); ++i)
    {
      for (int j = 0; j < nbrTestStates; ++j)
	{	  
	  double Tmp = 0.0;
	  for (int k = 0; k < nbrTestStates; ++k)
	    Tmp += overlaps[nbrOverlaps - 1][k] * overlaps[i][k];
	  overlaps[i][j] -= Tmp * overlaps[i][j];	  
	}
      bestOverlaps[i] = 0.0;
      for (int j = 0; j < nbrTestStates; ++j)
	{	  
	  bestOverlaps[i] += overlaps[i][j] * overlaps[i][j];
	}
    }
  SortOverlaps (bestOverlaps, overlaps, nbrOverlaps - 1, nbrTestStates);
}
