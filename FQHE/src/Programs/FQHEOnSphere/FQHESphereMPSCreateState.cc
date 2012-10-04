#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnDiskShort.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

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
Complex CreateLaughlinAMatrixElement (int laughlinIndex, unsigned long* partition1, unsigned long* partition2, int p1Level, int p2Level, int nValue);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereMPSCreateState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
//   (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
//   (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 20);
//   (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "p-truncation", "truncation level", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*OutputGroup) += new BooleanOption ('\n', "show-basis", "display the  truncated Hilbert space");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_interaction-name_n_nbrparticles_2s_nbrfluxquanta_lz_totallz.0.vec");
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

  int NbrParticles = 0;//Manager.GetInteger("nbr-particles"); 
  int NbrFluxQuanta = 0;//Manager.GetInteger("nbr-flux"); 
  int TotalLz = 0;//Manager.GetInteger("lz-value");
  int* ReferenceState = 0;
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

  if (Manager.GetBoolean("show-basis") == true)
    {
      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	{
	  cout << i << " : ";
	  Space->PrintState(cout, i);
	  cout << endl;
	}
      return 0;
    }
  
  ComplexMatrix* BMatrices = new ComplexMatrix[2];
  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [Manager.GetInteger("p-truncation") + 1];
  for (int i = 0; i <= Manager.GetInteger("p-truncation"); ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort(i, i, i + 1);
//       for (int j = 0; j < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++j)
// 	{
// 	  cout << j << " : ";
// 	  U1BosonBasis[i]->PrintState(cout, j);
// 	  cout << endl;
// 	}
    }
  CreateLaughlinBMatrices (3, BMatrices, U1BosonBasis, Manager.GetInteger("p-truncation"));
  ComplexVector State (Space->GetHilbertSpaceDimension(), true);
  Space->CreateStateFromMPSDescription(BMatrices, State, 0);

  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
    {
      cout << i << " : ";
      Space->PrintState(cout, i) << " = " << State[i] << endl;
    }

  return 0;
}

void CreateLaughlinBMatrices (int laughlinIndex, ComplexMatrix* bMatrices, BosonOnDiskShort** u1BosonBasis, int pLevel)
{
  int* StartingIndexPerPLevel = new int [pLevel + 1];
  int* NbrIndicesPerPLevel = new int [pLevel + 1];
  StartingIndexPerPLevel[0] = 0;
  int NbrNValue = ((2 * pLevel) + laughlinIndex);
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
  unsigned long* Partition1 = new unsigned long [pLevel + 1];
  unsigned long* Partition2 = new unsigned long [pLevel + 1];
  for (int i = 0; i <= pLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace1 = u1BosonBasis[i];
      int MaxN1 = (2 * i) + laughlinIndex;
      for (int j = 0; j <= pLevel; ++j)
	{
	  BosonOnDiskShort* TmpSpace2 = u1BosonBasis[j];
	  int MaxN2 = (2 * 2) + laughlinIndex;
	  int N1 = (2 * (i - j) - laughlinIndex + 1 + NbrNValue) / 2;
	  int N2 = (2 * (i - j) + laughlinIndex - 1 + NbrNValue) / 2;
	  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
	    {
	      TmpSpace1->GetMonomial(k1, Partition1);
	      for (int k2 = 0; k2 < TmpSpace1->GetHilbertSpaceDimension(); ++k2)
		{
		  TmpSpace2->GetMonomial(k2, Partition2);
		  bMatrices[1].SetMatrixElement(StartingIndexPerPLevel[i] + ((k1 * NbrNValue) + N1), StartingIndexPerPLevel[j] + ((k2 * NbrNValue) + N2), 
						CreateLaughlinAMatrixElement(laughlinIndex, Partition1, Partition2, k1, k2, - (N1 + N2 - NbrNValue) / 2));
		}
	    }
	}
    }
  delete[] Partition1;
  delete[] Partition2;
}

Complex CreateLaughlinAMatrixElement (int laughlinIndex, unsigned long* partition1, unsigned long* partition2, int p1Level, int p2Level, int nValue)
{
  Complex Tmp = 0.0;
  if (nValue != (p2Level - p1Level))
    return Tmp;
  return Tmp;
}
