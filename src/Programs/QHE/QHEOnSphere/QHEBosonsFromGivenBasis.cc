#include "config.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

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

#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereDeltaHamiltonian.h"

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

  ArchitectureManager Architecture;

  Manager += MainGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*MainGroup) += new SingleStringOption  ('\n', "input-file", "name of the file which contains definition of overlaps to evaluate");
  (*MainGroup) += new BooleanOption  ('\n', "lsort", "");
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
  
  
  char* VectorBasisName = new char [strlen(OverlapDefinition["Vectors"]) + 24];
  strcpy (VectorBasisName, OverlapDefinition["Vectors"]);
  char* VectorBasisName2 = VectorBasisName + strlen(VectorBasisName);

  int TotalMaxLz = LzMax * NbrParticles;
  RealVector** SortedTestVectors = new RealVector*[MaxNbrLz];
  bool Parity = true;
  if ((TotalMaxLz & 1) != 0)
    Parity = false;
  int* NbrSortedTestVectors = new int [TotalMaxLz];
  int* TestVectorPosition = new int [TotalMaxLz];

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
      RealVector* TestVectors = new RealVector[Degeneracy[i]]; 
      for (int j = 0; j < TotalMaxLz; ++j)  
	{  
	  TestVectorPosition[j] = TotalMaxLz;
	  NbrSortedTestVectors[j] = 0;
	}
      
      for (int j = 0; j < Degeneracy[i]; ++j)
       {
	 sprintf (VectorBasisName2, "%d.%d.vec", Lz, j);	      
	 if (TestVectors[j].ReadVector(VectorBasisName) == false)
	   {
	     cout << "error while reading " << VectorBasisName << endl;
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
	      TmpAngularMomentum = ((int) round(0.5 * (sqrt ((4.0 * TmpTriDiag.DiagonalElement(k)) + 1.0) - 1.0)));
	      cout << TmpAngularMomentum << " ";	      
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
	    int TmpAngularMomentum = ((int) round(0.5 * (sqrt ((4.0 * oper.MatrixElement(TestVectors[0], TestVectors[0]).Re) + 1.0) - 1.0)));
	    cout << "angular momentum of test eigenvectors = " << TmpAngularMomentum << endl;
	    NbrSortedTestVectors[TmpAngularMomentum]++;	    
	    TestVectorPosition[TmpAngularMomentum] = 0;
	  }

     Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());
     AbstractQHEOnSphereHamiltonian* Hamiltonian = new ParticleOnSphereDeltaHamiltonian(&Space, NbrParticles, LzMax, Architecture.GetArchitecture(), 1000000000l);
     for (int j = 0; j < MaxNbrLz; ++j)
       {
	 if (NbrSortedTestVectors[j] > 1)
	   {
	     int Shift = TestVectorPosition[j];
	     RealSymmetricMatrix HRep (NbrSortedTestVectors[j], NbrSortedTestVectors[j]);
	     for (int k = 0; k < NbrSortedTestVectors[j]; ++k)
	       {
		 for (int l = k; l < NbrSortedTestVectors[j]; ++l)
		   {	
		     HRep.SetMatrixElement(l, k, Hamiltonian->MatrixElement(DiagonalBasis[k + Shift], DiagonalBasis[l + Shift]).Re);
		   }
	       }
	     RealTriDiagonalSymmetricMatrix TmpTriDiag (NbrSortedTestVectors[j]);
	     HRep.Householder(TmpTriDiag, 1e-7);
	     TmpTriDiag.Diagonalize();
	     TmpTriDiag.SortMatrixUpOrder();
	     for (int k = 0; k < NbrSortedTestVectors[j]; ++k)
	       cout << j << " " << TmpTriDiag.DiagonalElement(k) << endl;
	   }
	 else
	   if (NbrSortedTestVectors[j] == 1)
	     {
	       int Shift = TestVectorPosition[j];
	       cout << j << " " << Hamiltonian->MatrixElement(DiagonalBasis[Shift], DiagonalBasis[Shift]).Re << endl;
	     }
       }

    }


  delete[] Degeneracy;
  delete[] VectorBasisName;
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
