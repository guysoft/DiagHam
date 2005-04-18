#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

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
  (*MainGroup) += new SingleDoubleOption  ('\n', "ortho-error", "scalar product value below which two states are considered as orthogonal", 1e-12);
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

  double OrthogonalityError = ((SingleDoubleOption*) Manager["ortho-error"])->GetDouble();

  ConfigurationParser OverlapDefinition;
  if (OverlapDefinition.Parse(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      OverlapDefinition.DumpErrors(cout) << endl;
      return -1;
    }

  int LzMax = 0;
  int NbrParticles = 0;
  if ((OverlapDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax < 0))
    {
      cout << "LzMax is not defined or as a wrong value" << endl;
      return -1;
    }
  if ((OverlapDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
    {
      cout << "NbrParticles is not defined or as a wrong value" << endl;
      return -1;
    }

  if (OverlapDefinition["Degeneracy"] == 0)
    {
      cout << "no Degeneracy defined in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }

  int MaxNbrLz;
  int* Degeneracy;
  if (OverlapDefinition.GetAsIntegerArray("Degeneracy", ' ', Degeneracy, MaxNbrLz) == false)
    {
      cout << "error while parsing Degeneracy in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  
  
  char* InputVectors = new char [strlen(OverlapDefinition["InputVectors"]) + 24];
  char* OutputVectors = new char [strlen(OverlapDefinition["OutputVectors"]) + 24];
  strcpy (InputVectors, OverlapDefinition["InputVectors"]);
  strcpy (OutputVectors, OverlapDefinition["OutputVectors"]);
  char* InputVectors2 = InputVectors + strlen(InputVectors);
  char* OutputVectors2 = OutputVectors + strlen(OutputVectors);

  int TotalMaxLz = LzMax * NbrParticles;
  RealVector** SortedTestVectors = new RealVector*[MaxNbrLz];
  int* NbrSortedTestVectors = new int [TotalMaxLz];
  bool Parity = true;
  if ((TotalMaxLz & 1) != 0)
    Parity = false;
  double** BestOverlaps = new double* [MaxNbrLz];
  int* NbrBestOverlaps = new int [MaxNbrLz];
  double** BestOverlapsPerL = new double* [MaxNbrLz];
  int** BestOverlapsPerLPosition = new int* [MaxNbrLz];
  int* NbrBestOverlapsPerL = new int [MaxNbrLz];
  int* TestVectorPosition = new int [TotalMaxLz];
  for (int i = 0; i < MaxNbrLz; ++i)
    {
      BestOverlapsPerL[i] = new double [TotalMaxLz];
      BestOverlapsPerLPosition[i] = new int [TotalMaxLz];
      NbrBestOverlapsPerL[i] = 0;
    }

  for (int i = 0; i < MaxNbrLz; ++i)
    {
      int Lz = i << 1;
      if (Parity == false)
	{
	  ++Lz;
	  cout << "Lz = " << Lz << "/2" << endl;
	}
      else
	{
	  cout << "Lz = " << i << endl;
	}
      BestOverlaps[i] = new double[Degeneracy[i]];
      double*** Overlaps = new double**[TotalMaxLz]; 
      int* NbrVectorWithMomentum = new int [TotalMaxLz]; 
      RealVector* TestVectors = new RealVector[Degeneracy[i]]; 
      for (int j = 0; j < TotalMaxLz; ++j)    
	{ 
	  TestVectorPosition[j] = TotalMaxLz;
	  NbrSortedTestVectors[j] = 0;
	  Overlaps[j] = new double*[Degeneracy[i]];
	  for (int l = 0; l < Degeneracy[i]; ++l)
	    Overlaps[j][l] = new double [Degeneracy[i]];
	  NbrVectorWithMomentum[j] = 0;
	}

      for (int j = 0; j < Degeneracy[i]; ++j)
	{
	  sprintf (OutputVectors2, "%d.%d.vec", Lz, j);	      
	  if (TestVectors[j].ReadVector(OutputVectors) == false)
	    {
	      cout << "error while reading " << OutputVectors << endl;
	      return -1;
	    }
	}

      RealMatrix DiagonalBasis (TestVectors, Degeneracy[i]);
      BosonOnSphere Space (NbrParticles, Lz, LzMax);
      if (Degeneracy[i] > 1)
	{
	  RealSymmetricMatrix HRep (Degeneracy[i], Degeneracy[i]);
	  ParticleOnSphereSquareTotalMomentumOperator oper(&Space, Lz, LzMax);
	  for (int k = 0; k < Degeneracy[i]; ++k)
	    {
	      for (int l = k; l < Degeneracy[i]; ++l)
		{	
		  HRep.SetMatrixElement(l, k,  oper.MatrixElement(TestVectors[k], TestVectors[l]).Re);
		}
	    }
	  cout << endl << endl;
	  
	  RealMatrix TmpEigenvector (Degeneracy[i], Degeneracy[i], true);
	  for (int l = 0; l < Degeneracy[i]; ++l)
	    TmpEigenvector(l, l) = 1.0;
	  RealTriDiagonalSymmetricMatrix TmpTriDiag (Degeneracy[i]);
	  HRep.Householder(TmpTriDiag, 1e-7, TmpEigenvector);
	  TmpTriDiag.Diagonalize(TmpEigenvector);
	  TmpTriDiag.SortMatrixUpOrder(TmpEigenvector);
	  cout << "angular momentum of test eigenvectors = ";
	  int TmpAngularMomentum;
	  for (int k = 0; k < Degeneracy[i]; ++k)	    
	    {
	      TmpAngularMomentum = ((int) round((sqrt ((4.0 * TmpTriDiag.DiagonalElement(k)) + 1.0) - 1.0)));
	      if ((TmpAngularMomentum & 1) == 0)
		{
		  TmpAngularMomentum >>= 1;
		  cout << TmpAngularMomentum << "(" << TmpTriDiag.DiagonalElement(k) << ")" << " ";	      
		}
	      else
		{		  
		  cout << TmpAngularMomentum << "/2 ";	  
		  TmpAngularMomentum -= 1;
		  TmpAngularMomentum >>= 1;		  
		}
	      NbrSortedTestVectors[TmpAngularMomentum]++;
	      if (TestVectorPosition[TmpAngularMomentum] > k)
		TestVectorPosition[TmpAngularMomentum] = k;
	    }
	  
	  cout << endl;
	  DiagonalBasis.Multiply(TmpEigenvector);
	}
      else
	if (Degeneracy[i] == 1)
	  {
	    ParticleOnSphereSquareTotalMomentumOperator oper(&Space, Lz, LzMax);
	    int TmpAngularMomentum = ((int) round((sqrt ((4.0 * oper.MatrixElement(TestVectors[0], TestVectors[0]).Re) + 1.0) - 1.0)));
	    cout << "angular momentum of test eigenvectors = ";
	    if ((TmpAngularMomentum & 1) == 0)
	      {
		TmpAngularMomentum >>= 1;
		cout << TmpAngularMomentum << " " << endl;	      
	      }
	    else
	      {		  
		cout << TmpAngularMomentum << "/2" << endl;	  
		TmpAngularMomentum -= 1;
		TmpAngularMomentum >>= 1;		  
	      }
	    NbrSortedTestVectors[TmpAngularMomentum]++;	    
	    TestVectorPosition[TmpAngularMomentum] = 0;
	  }

      int NbrNonZeroOverlap = 0;
      for (int j = 0; j < Degeneracy[i]; ++j)
	{
	  sprintf (InputVectors2, "%d.%d.vec", Lz, j);
	  cout << InputVectors << ":" << endl;
	  RealVector ReferenceVector;	  
	  if (ReferenceVector.ReadVector(InputVectors) == false)
	    {
	      cout << "error while reading " << InputVectors << endl;
	      return -1;
	    }
	  
	  ParticleOnSphereSquareTotalMomentumOperator oper(&Space, Lz, LzMax);
	  int AngularMomentum = ((int) round((sqrt ((4.0 * oper.MatrixElement(ReferenceVector, ReferenceVector).Re) + 1.0) - 1.0)));
	  cout << "angular momentum = ";
	    if ((AngularMomentum & 1) == 0)
	      {
		AngularMomentum >>= 1;
		cout << AngularMomentum << " " << endl;	      
	      }
	    else
	      {		  
		cout << AngularMomentum << "/2" << endl;	  
		AngularMomentum -= 1;
		AngularMomentum >>= 1;		  
	      }
	 
	  if ((AngularMomentum < MaxNbrLz) && (NbrSortedTestVectors[AngularMomentum] > 0))
	    {
	      BestOverlaps[i][NbrNonZeroOverlap] = 0.0;
	      int Shift = TestVectorPosition[AngularMomentum];
	      for (int k = 0; k < NbrSortedTestVectors[AngularMomentum]; ++k)
		{
		  if (ReferenceVector.GetVectorDimension() != DiagonalBasis[k + Shift].GetVectorDimension())
		    {
		      sprintf (OutputVectors2, "%d.%d.vec", Lz, k + Shift);
		      cout << "dimension mismatch between " << InputVectors << " and " << OutputVectors << endl;
		      return -1;
		    }
		  
		  double Scalar = ReferenceVector * DiagonalBasis[k + Shift];
		  Overlaps[AngularMomentum][NbrVectorWithMomentum[AngularMomentum]][k] = Scalar;
		  if (fabs(Scalar) > OrthogonalityError)
		    cout << "    " << OutputVectors << "  = " << Scalar << endl;
		  BestOverlaps[i][NbrNonZeroOverlap] += Scalar * Scalar;
		}
	      NbrVectorWithMomentum[AngularMomentum]++;
	      BestOverlaps[i][NbrNonZeroOverlap] = sqrt(BestOverlaps[i][j]);
	      BestOverlapsPerL[AngularMomentum][NbrBestOverlapsPerL[AngularMomentum]] = BestOverlaps[i][NbrNonZeroOverlap];
	      if ((2 * AngularMomentum) == Lz)
		{
		  BestOverlapsPerLPosition[AngularMomentum][NbrBestOverlapsPerL[AngularMomentum]] = j;
		}
	      else
		{
                  BestOverlapsPerLPosition[AngularMomentum][NbrBestOverlapsPerL[AngularMomentum]] = -1;
		}
	      NbrBestOverlapsPerL[AngularMomentum]++;
	      cout << endl << "best overlap = " << BestOverlaps[i][NbrNonZeroOverlap] << endl;
	      ++NbrNonZeroOverlap;
	    }
	  else
	    {
	      cout << "no possible overlap " << endl;
	    }
	  SortArrayDownOrdering(BestOverlaps[i], NbrNonZeroOverlap);
	  NbrBestOverlaps[i] = NbrNonZeroOverlap;
	}

/*      SortOverlaps(BestOverlaps, Overlaps, Degeneracy[i], Degeneracy[i]);
      cout << endl << "overlaps = " << endl;
      for (int j = 0; j < Degeneracy[i]; ++j)
	{
	  cout << BestOverlaps[j] << endl;
	}*/


      cout << "----------------------------------------------------" << endl;

       for (int j = 0; j < TotalMaxLz; ++j)    
	 if (NbrVectorWithMomentum[j] > 0)
	   {
	     cout << j << ": ";
	     for (int k = 0; k < NbrVectorWithMomentum[j]; ++k)
	       for (int l = k + 1; l < NbrVectorWithMomentum[j]; ++l)
		 {
		   double Scalar = 0.0;
		   for (int m = 0; m < NbrSortedTestVectors[j]; ++m)
		     {
		       Scalar += Overlaps[j][k][m] * Overlaps[j][l][m];
		     }
		   cout << Scalar << " ";
		 }
	     cout << endl;
	   }
      cout << "----------------------------------------------------" << endl;

      for (int j = 0; j < Degeneracy[i]; ++j)
	{
	  for (int l = 0; l < Degeneracy[i]; ++l)
	    delete[] Overlaps[j][l];
	  delete[] Overlaps[j];
	}
      delete[] Overlaps;
//      delete[] TestVectors;
    }
  
  
  cout << "----------------------------------------------------" << endl;
  cout << "final results" << endl;
  cout << "----------------------------------------------------" << endl;
  cout << "Lz sorted results " << endl << endl;
  for (int i = 0; i < MaxNbrLz; ++i)
    {
      int Lz = i << 1;
      if (Parity == false)
	{
	  ++Lz;
	  cout << "Lz = " << Lz << "/2 : ";
	}
      else
	{
	  cout << "Lz = " << i <<" : ";
	}
      for (int j = 0; j < NbrBestOverlaps[i]; ++j)
	cout << BestOverlaps[i][j] << " ";
      cout << endl;
   }
  cout << "----------------------------------------------------" << endl;
  cout << "L sorted results " << endl << endl;
  for (int i = 0; i < MaxNbrLz; ++i)
    {
      if (NbrBestOverlapsPerL[i] > 0)
	{
	  SortArrayDownOrdering(BestOverlapsPerL[i], BestOverlapsPerLPosition[i], NbrBestOverlapsPerL[i]);
	  int Lz = i << 1;
	  if (Parity == false)
	    {
	      ++Lz;
	      cout << "L = " << Lz << "/2 : ";
	    }
	  else
	    {
	      cout << "L = " << i <<" : ";
	    }
	  double CurrentOverlap = BestOverlapsPerL[i][0];
	  cout << CurrentOverlap;
	  int NbrOverlaps = 1;
	  int Pos = BestOverlapsPerLPosition[i][0];
	  for (int j = 1; j <  NbrBestOverlapsPerL[i]; ++j)
	    {
	      if (fabs (BestOverlapsPerL[i][j] - CurrentOverlap) > OrthogonalityError)
		{
		  cout << "(" << ((2 * (NbrOverlaps - 1)) + 1 ) << ", " << Pos <<  ") ";
		  CurrentOverlap = BestOverlapsPerL[i][j];
		  Pos = BestOverlapsPerLPosition[i][j];
		  cout << CurrentOverlap;
		  NbrOverlaps = 1;
		}
	      else
		{
		  ++NbrOverlaps;
		  if (Pos < BestOverlapsPerLPosition[i][j])
		    Pos = BestOverlapsPerLPosition[i][j];
		}
	    }
	  cout << "(" << ((2 * (NbrOverlaps - 1)) + 1 ) << ", " << Pos << ") " << endl;
	}
   }


  delete[] NbrBestOverlaps;
  delete[] NbrBestOverlapsPerL;
  for (int i = 0; i < MaxNbrLz; ++i)
    {
      delete[] BestOverlapsPerLPosition[i];
      delete[] BestOverlapsPerL[i];
    }
  delete[] BestOverlapsPerLPosition;
  delete[] BestOverlapsPerL;
  delete[] Degeneracy;
  delete[] InputVectors;
  delete[] OutputVectors;
  delete[] SortedTestVectors;
  delete[] NbrSortedTestVectors;
  delete[] TestVectorPosition;
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
