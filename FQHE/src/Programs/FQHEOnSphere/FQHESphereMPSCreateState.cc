#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnDiskShort.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"
#include "MathTools/FactorialCoefficient.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHEMPSCreateStateOperation.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Matrix/SparseComplexMatrix.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


void CreateLaughlinBMatrices (int laughlinIndex, ComplexMatrix* bMatrices, BosonOnDiskShort** u1BosonBasis, int pLevel);
Complex CreateLaughlinAMatrixElement (int laughlinIndex, unsigned long* partition1, unsigned long* partition2, int p1Level, int p2Level, int nValue, FactorialCoefficient& coef);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereMPSCreateState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

//   (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
//   (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 20);
//   (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "p-truncation", "truncation level", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-index", "index of the Laughlin state to generate", 3);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the MPS state into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the MPS state into a text file");
  (*OutputGroup) += new BooleanOption ('n', "normalize-sphere", "express the MPS in the normalized sphere basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMPSCreateState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  int TotalLz = 0;
  int* ReferenceState = 0;
  char* OutputFileName = Manager.GetString("bin-output");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  if ((OutputTxtFileName == 0) && (OutputFileName == 0))
    {
      cout << "error, an output file (binary or text) has to be provided" << endl;
      return 0;
    }


  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
    return -1;

//  ParticleOnSphere* Space = 0;
  FermionOnSpherePTruncated* Space = 0;
  if (Manager.GetBoolean("boson") == true)
    {
      cout << "bosons are not yet implemented" << endl;
      return 0;
    }
  else
    {
#ifdef __64_BITS__
	  if (NbrFluxQuanta <= 62)
#else
	  if (NbrFluxQuanta <= 30)
#endif
	    {
	      Space = new FermionOnSpherePTruncated(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
	    }
	  else
	    {
#ifdef __128_BIT_LONGLONG__
	      if (NbrFluxQuanta <= 126)
#else
		if (NbrFluxQuanta <= 62)
#endif
		  {
		    Space = 0;//new FermionOnSpherePTruncatedLong(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
		  }
		else
		  {
#ifdef __128_BIT_LONGLONG__
		    cout << "cannot generate an Hilbert space when nbr-flux > 126" << endl;
#else
		    cout << "cannot generate an Hilbert space when nbr-flux > 62" << endl;
#endif
		    return 0;
		  }
	    }

   }

  cout << "Hilbert space dimension : " << Space->GetLargeHilbertSpaceDimension() << endl;
  
  int LaughlinIndex = Manager.GetInteger("laughlin-index");
  int NbrBMatrices = 2;
  ComplexMatrix* BMatrices = new ComplexMatrix[NbrBMatrices];
  SparseComplexMatrix* SparseBMatrices = new SparseComplexMatrix[NbrBMatrices];
  SparseComplexMatrix* SparseConjugateBMatrices = new SparseComplexMatrix[NbrBMatrices];
  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [Manager.GetInteger("p-truncation") + 1];
  for (int i = 0; i <= Manager.GetInteger("p-truncation"); ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort(i, i, Manager.GetInteger("p-truncation") + 1);
    }
  CreateLaughlinBMatrices (LaughlinIndex, BMatrices, U1BosonBasis, Manager.GetInteger("p-truncation"));

  for (int i = 0; i < NbrBMatrices; ++i)
    {
      SparseBMatrices[i] = BMatrices[i];
      SparseConjugateBMatrices[i] = SparseBMatrices[i].HermitianTranspose();
    }
  cout << "B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;

  SparseComplexMatrix** SparseTensorProductBMatrices = new SparseComplexMatrix*[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      SparseTensorProductBMatrices[i] = new SparseComplexMatrix[NbrBMatrices];
      for (int j = 0; j < NbrBMatrices; ++j)
	SparseTensorProductBMatrices[i][j] = TensorProduct(SparseBMatrices[i], SparseConjugateBMatrices[j]);
    }
  
  RealVector State (Space->GetHilbertSpaceDimension(), true);


  FQHEMPSCreateStateOperation Operation(Space, SparseBMatrices, &State, Manager.GetInteger("p-truncation") + ((LaughlinIndex - 1) / 2));
  Operation.ApplyOperation(Architecture.GetArchitecture());
  //    Space->CreateStateFromMPSDescription(SparseBMatrices, State, Manager.GetInteger("p-truncation") + ((LaughlinIndex - 1) / 2));

  if (Manager.GetBoolean("normalize-sphere"))
    Space->ConvertFromUnnormalizedMonomial(State);
  
  if (OutputTxtFileName != 0)
    {
      ofstream File;
      File.open(OutputTxtFileName, ios::binary | ios::out);
      File.precision(14);	
      for (long i = 0; i < Space->GetLargeHilbertSpaceDimension(); ++i)
	{
	  State.PrintComponent(File, i) << " ";
	  Space->PrintStateMonomial(File, i) << endl;
	}
      File.close();
    }
  if (OutputFileName != 0)
    {
      State.WriteVector(OutputFileName);
    }
  
  return 0;
}

void CreateLaughlinBMatrices (int laughlinIndex, ComplexMatrix* bMatrices, BosonOnDiskShort** u1BosonBasis, int pLevel)
{
  int* StartingIndexPerPLevel = new int [pLevel + 1];
  int* NbrIndicesPerPLevel = new int [pLevel + 1];
  StartingIndexPerPLevel[0] = 0;
  int NbrNValue = ((2 * pLevel) + laughlinIndex);
  int NValueShift = NbrNValue - 1;
  NbrIndicesPerPLevel[0] = u1BosonBasis[0]->GetHilbertSpaceDimension() * NbrNValue;
  for (int i = 1; i <= pLevel; ++i)
    {
      StartingIndexPerPLevel[i] = StartingIndexPerPLevel[i - 1] + NbrIndicesPerPLevel[i - 1];
      NbrIndicesPerPLevel[i] = u1BosonBasis[i]->GetHilbertSpaceDimension()  * NbrNValue;
    }
  int MatrixSize = NbrIndicesPerPLevel[pLevel] + StartingIndexPerPLevel[pLevel];

  bMatrices[0] = ComplexMatrix(MatrixSize, MatrixSize, true);
  for (int i = 0; i <= pLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = u1BosonBasis[i];
      for (int j = 1; j < NbrNValue; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    {
	      bMatrices[0].SetMatrixElement(StartingIndexPerPLevel[i] + ((k * NbrNValue) + j - 1), StartingIndexPerPLevel[i] + ((k * NbrNValue) + j), 1.0);
	    }
	}
    }

  bMatrices[1] = ComplexMatrix(MatrixSize, MatrixSize, true);
  unsigned long* Partition1 = new unsigned long [pLevel + 2];
  unsigned long* Partition2 = new unsigned long [pLevel + 2];
  FactorialCoefficient Coef;

  for (int i = 0; i <= pLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace1 = u1BosonBasis[i];
      int MaxN1 = (2 * i) + laughlinIndex;
      for (int j = 0; j <= pLevel; ++j)
	{
	  BosonOnDiskShort* TmpSpace2 = u1BosonBasis[j];
	  int MaxN2 = (2 * 2) + laughlinIndex;
	  int N1 = (2 * (j - i) + laughlinIndex - 1 + NValueShift) / 2;
	  int N2 = (2 * (j - i) - laughlinIndex + 1 + NValueShift) / 2;
	  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
	    {
	      TmpSpace1->GetOccupationNumber(k1, Partition1);
	      for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		{
		  TmpSpace2->GetOccupationNumber(k2, Partition2);
		  Complex Tmp = CreateLaughlinAMatrixElement(laughlinIndex, Partition1, Partition2, i, j, - (N1 + N2 - NValueShift) / 2, Coef);
		  bMatrices[1].SetMatrixElement(StartingIndexPerPLevel[i] + ((k1 * NbrNValue) + N1), StartingIndexPerPLevel[j] + ((k2 * NbrNValue) + N2), Tmp);
//		  cout << i << " " << j << " | " << k1 << " " << k2 << " | " << N1 << " " << N2 << " " << (-(N1 + N2 - NValueShift) / 2) << " : " << Tmp << endl;
		}
	    }
	}
    }
  delete[] Partition1;
  delete[] Partition2;
}

Complex CreateLaughlinAMatrixElement (int laughlinIndex, unsigned long* partition1, unsigned long* partition2, int p1Level, int p2Level, int nValue, FactorialCoefficient& coef)
{
  Complex Tmp = 1.0;
  if (nValue != (p1Level - p2Level))
    {
      Tmp = 0.0;
      return Tmp;
    }
  int PMax = p1Level;
  if (p2Level > p1Level)
    PMax = p2Level;
  for (int i = 1; i <= PMax; ++i)
    {
      Complex Tmp2 = 0.0;
//      cout << partition1[i] << " " << partition2[i] << " " << i << " " << PMax << endl;
      for (int j = 0; j <= partition1[i]; ++j)
	{
	  int k = partition2[i] + j - partition1[i];
	  if ((k >= 0) && (k <= partition2[i]))
	    {
	      int Sum = k + j;
	      coef.SetToOne();
	      coef.PartialFactorialMultiply(partition1[i] - j + 1, partition1[i]);
	      coef.PartialFactorialMultiply(partition2[i] - k + 1, partition2[i]);
	      coef.FactorialDivide(j);
	      coef.FactorialDivide(k);
	      coef.FactorialDivide(j);
	      coef.FactorialDivide(k);
	      coef.PowerNMultiply(laughlinIndex, Sum);
	      coef.PowerNDivide(i, Sum);
//	      cout << "Sum=" << Sum << endl;
	      switch  (Sum & 0x3)
		{
		case 0:
		  Tmp2.Re += sqrt(coef.GetNumericalValue());
		  break;
		case 1:
		  Tmp2.Im += sqrt(coef.GetNumericalValue());
		  break;
		case 2:
		  Tmp2.Re -= sqrt(coef.GetNumericalValue());
		  break;
		case 3:
		  Tmp2.Im -= sqrt(coef.GetNumericalValue());
		  break;
		}
	    }
	}
//      cout << Tmp2 << endl;
      Tmp *= Tmp2;
    }
  return Tmp;
}
