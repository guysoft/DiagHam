#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Polynomial/IntegerPolynomial.h"
#include "Polynomial/QDeformedBinomialCoefficients.h"

#include "Tools/QHE/QHESpectrum/MomentumMultipletSet.h"


#include <iostream>
#include <stdlib.h>
#include <math.h>

using std::cout;
using std::endl;


// maximum k value size in bits
#define K_MAXSIZE 4
// shift and mask associated to the maximum k value
#ifdef __64_BITS__
#define K_SHIFT 4
#define K_MASK 0x15l
#else
#define K_SHIFT 3
#define K_MASK 0x7l
#endif

// evaluate Hilbert space dimension for k-type parafermion in a given number of states according to the rules  given by Gurarie
//
// nbrParafermions = number of parafermions
// kType = parafermion type
// nbrStates = number of states that can be occupied by the parafermions
long EvaluateHilbertSpaceDimension(int nbrParafermions, int kType, int nbrStates);

// evaluate Hilbert space dimension for k-type parafermion without any constraint
//
// nbrParafermions = number of parafermions
// kType = parafermion type
// nbrStates = number of states that can be occupied by the parafermions
// return value = Hilbert space dimension
long EvaluateUnconstraintHilbertSpaceDimension (int nbrParafermions, int kType, int nbrStates);

// evaluate Hilbert space dimension for k-type parafermion without any constraint
//
// stateDescription = array that contains state description
// nbrParafermions = number of parafermions
// kType = parafermion type
// position = index of the current state description
// statePosition = index of the one-particle state that has to be currently filled (number of state - 1 except during recursion)
// return value = new index of the current state description
int EvaluateUnconstraintStates (long** stateDescription, int nbrParafermions, int kType, int position, int statePosition);

// get the number of ways to put F parafermions in n box
// 
// nbrParafermions = number of parafermions
// kType = number of parafermion species
// nbrBoxes = number of boxes
// binomialCoeffients = array that contains all usefull binomial coefficients
// return value = number of way
long GetParafermionPartitionNumber (int nbrParafermions, int nbrBoxes, int kType, long** binomialCoeffients);

// get the number of ways to write sum_i=1^n i a_i = F where F and n (aka number of possibilities to have a state of total momentum F with a 
// free number of bosons and n one-body state carrying each a Lz momentum ranging from 1 to n)
// 
// momentum = total momemtum (aka F)
// nbrStates = number of states (aka n)
// return value = number of way
long GetFixedLzFreeNumberBosonPartitionNumber (int momentum, int nbrStates);

// get all ways to write sum_i=1^n i a_i = F where F and n (aka number of possibilities to have a state of total momentum F with a 
// free number of bosons and n one-body state carrying each a Lz momentum ranging from 1 to n)
// 
// momentum = total momemtum (aka F)
// nbrStates = number of states (aka n)
// permutations = reference on the array where permutations will be stored
// partitionNumber = reference on the integer where the corresponding partition will be stored
void GetFixedLzFreeNumberBosonPermutations (int momentum, int nbrStates, int**& permutations, long& partitionNumber);

// recursive function associated to GetFixedLzFreeNumberBosonPermutations
// 
// momentum = total momemtum
// nbrStates = number of states
// permutations = reference on the array where permutations will be stored
// position = current position in the array
// return value = new current position in the array
int RecursiveFixedLzFreeNumberBosonPermutations (int momentum, int nbrStates, int**& permutations, int position);

// get the number of ways to put F parafermions in n box
// 
// nbrParafermions = number of parafermions
// kType = number of parafermion species
// nbrBoxes = number of boxes
// binomialCoeffients = array that contains all usefull binomial coefficients
// return value = number of way
void GetMomentumMultipletDecomposition (int nbrParticles, int nbrQuasiholes, int kType, QDeformedBinomialCoefficients& binomialCoeffients);


int main(int argc, char** argv)
{
  OptionManager Manager ("ParafermionQuasiholeDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('k', "k-value", "number of particles", 2);
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 10);
  (*SystemGroup) += new SingleIntegerOption  ('q', "nbr-quasiholes", "number of quasiholes", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "l-dimension", "get dimensions of all subspaces with fixed total l value");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type ParafermionQuasiholeDimension -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int KValue = ((SingleIntegerOption*) Manager["k-value"])->GetInteger(); 
  int NbrParticles  = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrQuasiholes = ((SingleIntegerOption*) Manager["nbr-quasiholes"])->GetInteger(); 

  if ((NbrParticles % KValue) != 0)
    {
      cout << "the number of particles must be an multiple of k" << endl;
      return 1;
    }
  if ((NbrQuasiholes % KValue) != 0)
    {
      cout << "the number of quasiholes must be an multiple of k" << endl;
      return 1;
    }
  
  int MaxBinomial = (NbrParticles / KValue) + NbrQuasiholes;
  int Tmp = NbrParticles + KValue - 2;
  if (Tmp > MaxBinomial)
    MaxBinomial = Tmp;
  Tmp = NbrQuasiholes + 2 * NbrParticles * ((KValue - 1) * (KValue - 1));
  if (Tmp > MaxBinomial)
    MaxBinomial = Tmp;
  long** BinomialCoefficients = GetBinomialCoefficients(MaxBinomial);

//  MaxBinomial = 8;
  QDeformedBinomialCoefficients QDeformed(MaxBinomial);
  
//   for (int i = 0; i <= MaxBinomial; ++i)
//     {
//       for (int j = 0; j <= i; ++j)
// 	cout << QDeformed(i, j) << " | ";
//       cout << endl;
//     }

  IntegerPolynomial TmpPoly (3);
//  TmpPoly[4] = 1l;
  TmpPoly[3] = 1l;
  TmpPoly[2] = 2l;
  TmpPoly[1] = 2l;
  TmpPoly[0] = 1l;
  cout << TmpPoly << endl;
  MomentumMultipletSet Multiplet1(TmpPoly);
  cout << Multiplet1 << endl;
  MomentumMultipletSet Multiplet2(12);
  Multiplet2.FindMultipletsForBosons(6, 2) ;
  cout << Multiplet2 << endl;
 
//  GetMomentumMultipletDecomposition(NbrParticles, NbrQuasiholes, KValue, QDeformed);

  
  long HilbertSpaceDimension = 0;
  int NbrUnclusteredParafermions = 0;
  while (NbrUnclusteredParafermions <= NbrParticles)
    {
//       cout << GetParafermionPartitionNumber(NbrUnclusteredParafermions, NbrQuasiholes / KValue, KValue, BinomialCoefficients) << " " 
// 	   << BinomialCoefficients[NbrQuasiholes + ((NbrParticles - NbrUnclusteredParafermions) / KValue)][NbrQuasiholes]<< " " 
// 	   << (NbrQuasiholes + ((NbrParticles - NbrUnclusteredParafermions) / KValue)) << " " << NbrQuasiholes << endl;
      HilbertSpaceDimension += (GetParafermionPartitionNumber(NbrUnclusteredParafermions, NbrQuasiholes / KValue, KValue, BinomialCoefficients) * 
				BinomialCoefficients[NbrQuasiholes + ((NbrParticles - NbrUnclusteredParafermions) / KValue)][NbrQuasiholes]);
      NbrUnclusteredParafermions += KValue;
    }
  cout << HilbertSpaceDimension << endl;
}


// evaluate Hilbert space dimension for k-type parafermion in a given number of states according to the rules  given by Gurarie
//
// nbrParafermions = number of parafermions
// kType = parafermion type
// nbrStates = number of states that can be occupied by the parafermions

long EvaluateHilbertSpaceDimension(int nbrParafermions, int kType, int nbrStates)
{
  long MaximumDimension = EvaluateUnconstraintHilbertSpaceDimension(nbrParafermions, kType, nbrStates);
  long** StateDescription = new long* [MaximumDimension];
  int NbrBlock = nbrStates >> K_SHIFT;
  if ((NbrBlock & K_MASK) != 0)
    ++NbrBlock;
  for (int i = 0; i < MaximumDimension; ++i)
    {
      StateDescription[i] = new long [NbrBlock];
      for (int j = 0; j < NbrBlock; ++j)
	StateDescription[i][j] = 0l;
    }
  EvaluateUnconstraintStates(StateDescription, 0, nbrParafermions, kType, nbrStates - 1);
  long Dimension = 0;
  for (int i = 0; i < MaximumDimension; ++i)
    {
      ++Dimension;
      delete[] StateDescription[i];      
    }
  return Dimension;
}

// evaluate Hilbert space dimension for k-type parafermion without any constraint
//
// nbrParafermions = number of parafermions
// kType = parafermion type
// nbrStates = number of states that can be occupied by the parafermions
// return value = Hilbert space dimension

long EvaluateUnconstraintHilbertSpaceDimension (int nbrParafermions, int kType, int nbrStates)
{
  if (nbrParafermions == 0)
    return 0l;
  if (nbrParafermions == 1)
    return (long) nbrStates;
  if (nbrStates == 1)
    if (nbrParafermions >= kType)
      return 0l;
    else
      return 1l;

  int Max = nbrParafermions;
  if (kType < Max)
    Max = kType - 1;
  long Tmp = 0;
  for (int i = 0; i <= Max; ++i)
    Tmp += EvaluateUnconstraintHilbertSpaceDimension(nbrParafermions - i, kType, nbrStates - 1);
  return Tmp;
}

// evaluate Hilbert space dimension for k-type parafermion without any constraint
//
// stateDescription = array that contains state description
// nbrParafermions = number of parafermions
// kType = parafermion type
// position = index of the current state description
// statePosition = index of the one-particle state that has to be currently filled (number of state - 1 except during recursion)
// return value = new index of the current state description

int EvaluateUnconstraintStates (long** stateDescription, int nbrParafermions, int kType, int position, int statePosition)
{
  if (statePosition == 0)
    if (nbrParafermions >= kType)
      return position;
    else
      {
	stateDescription[position][0] |= ((long) nbrParafermions);
	return (position + 1);
      }
  if (nbrParafermions == 0)
    return position;
  if (nbrParafermions == 1)
    {
      for (; statePosition >= 0; --statePosition)
	{
	  stateDescription[position][statePosition >> K_SHIFT] |= 1l << ((statePosition & K_MASK) * K_MAXSIZE);
	  ++position;
	}
      return position;
    }

  int Max = nbrParafermions;
  if (kType < Max)
    Max = kType - 1;
  int TmpPos = 0;
  int Shift = (statePosition & K_MASK) * K_MAXSIZE;
  int Index = statePosition >> K_SHIFT;
  for (int i = 0; i <= Max; ++i)
    {
      TmpPos = EvaluateUnconstraintStates(stateDescription, nbrParafermions - i, kType, position, statePosition - 1);
      for (; position < TmpPos; ++position)
	stateDescription[position][Index] |= ((long) i) << Shift;
    }
  return position;
}


// get the number of ways to put F parafermions in n box
// 
// nbrParafermions = number of parafermions
// kType = number of parafermion species
// nbrBoxes = number of boxes
// binomialCoeffients = array that contains all usefull binomial coefficients
// return value = number of way

long GetParafermionPartitionNumber (int nbrParafermions, int nbrBoxes, int kType, long** binomialCoeffients)
{
  if ((nbrParafermions == 0) || (kType == 1))
    return 1l;
  if (kType == 2)
    {
      if (nbrParafermions <= nbrBoxes)
        return binomialCoeffients[nbrBoxes][nbrParafermions];
      else
	return 0l;
    }
  int** Permutations;
  long NbrPermutations;
  GetFixedLzFreeNumberBosonPermutations(nbrParafermions, kType - 1, Permutations, NbrPermutations);
  long NbrWays = 0l;
  long Tmp = 1l;
  int* TmpPermutation;
  int TmpCoefficient;
  for (int i = 0; i < NbrPermutations; ++i)
    {
      TmpPermutation = Permutations[i];
      Tmp = 1l;
      for (int j = 1; ((j < kType) && (Tmp != 0l)); ++j)
	{
	  int k = 1;
	  TmpCoefficient = (j * nbrBoxes * kType) + ((kType - (2 * (kType - j) * j)) * TmpPermutation[j - 1]);
	  for (; k < j; ++k)
	    TmpCoefficient -= (2 * (kType - j) * k) * TmpPermutation[k - 1];
	  ++k;
	  for (; k < kType; ++k)
	    TmpCoefficient -= (2 * (kType - k) * j) * TmpPermutation[k - 1];
	  if ((TmpCoefficient >= 0) && ((TmpCoefficient % kType) == 0) && ((TmpPermutation[j - 1] * kType) <= TmpCoefficient))
	    {
	      Tmp *= binomialCoeffients[TmpCoefficient / kType][TmpPermutation[j - 1]];
	      cout << binomialCoeffients[TmpCoefficient / kType][TmpPermutation[j - 1]] << endl;
	    }
	  else
	    Tmp = 0l;
	}
      cout << "final tmp=" << Tmp << endl;
      NbrWays += Tmp;
    }
  cout << "final NbrWays =" << Tmp << endl;
  for (int i = 0; i < NbrPermutations; ++i)
    delete[] Permutations[i];
  delete[] Permutations;
  return NbrWays;
}

// get the number of ways to write sum_i=1^n i a_i = F where F and n (aka number of possibilities to have a state of total momentum F with a 
// free number of bosons and n one-body state carrying each a Lz momentum ranging from 1 to n)
// 
// momentum = total momemtum (aka F)
// nbrStates = number of states (aka n)
// return value = number of way

long GetFixedLzFreeNumberBosonPartitionNumber (int momentum, int nbrStates)
{
  if ((nbrStates == 1) || (momentum == 0))
    return 1l;
  int Max = momentum / nbrStates;
  long Tmp = 0l;
  for (int i = 0; i <= Max; ++i)
    Tmp += GetFixedLzFreeNumberBosonPartitionNumber(momentum - (i * nbrStates), nbrStates - 1);
  return Tmp;
}

// get all ways to write sum_i=1^n i a_i = F where F and n (aka number of possibilities to have a state of total momentum F with a 
// free number of bosons and n one-body state carrying each a Lz momentum ranging from 1 to n)
// 
// momentum = total momemtum (aka F)
// nbrStates = number of states (aka n)
// permutations = reference on the array where permutations will be stored
// partitionNumber = reference on the integer where the corresponding partition will be stored

void GetFixedLzFreeNumberBosonPermutations (int momentum, int nbrStates, int**& permutations, long& partitionNumber)
{
  partitionNumber = GetFixedLzFreeNumberBosonPartitionNumber(momentum, nbrStates);
  permutations = new int* [partitionNumber];
  for (int i = 0; i < partitionNumber; ++i)
    {
      permutations[i] = new int [nbrStates];
      for (int j = 0; j < nbrStates; ++j)
	permutations[i][j] = 0;
    }
  RecursiveFixedLzFreeNumberBosonPermutations(momentum, nbrStates, permutations, 0);
}

// recursive function associated to GetFixedLzFreeNumberBosonPermutations
// 
// momentum = total momemtum
// nbrStates = number of states
// permutations = reference on the array where permutations will be stored
// position = current position in the array
// return value = new current position in the array

int RecursiveFixedLzFreeNumberBosonPermutations (int momentum, int nbrStates, int**& permutations, int position)
{
  if (nbrStates == 1)
    {
      permutations[position][nbrStates - 1] = momentum;
      return position + 1;
    }
  if (momentum == 0)
    {
      return position + 1;
    }    
  int TmpPosition;
  int TmpMomentum = momentum;
  int Tmp = 0;
  while (TmpMomentum >= 0)
    {
      TmpPosition = RecursiveFixedLzFreeNumberBosonPermutations(TmpMomentum, nbrStates - 1, permutations, position);
      for (; position < TmpPosition; ++position)
	permutations[position][nbrStates - 1] = Tmp;
      TmpMomentum -= nbrStates;
      ++Tmp;
    }
  return position;
}

// get the number of ways to put F parafermions in n box
// 
// nbrParafermions = number of parafermions
// kType = number of parafermion species
// nbrBoxes = number of boxes
// binomialCoeffients = array that contains all usefull binomial coefficients
// return value = number of way

void GetMomentumMultipletDecomposition (int nbrParticles, int nbrQuasiholes, int kType, QDeformedBinomialCoefficients& binomialCoeffients)
{
  int NbrUnclusteredParafermions = 0;
  int* TmpPermutation;
  int TmpCoefficient;
  int** Permutations;
  long NbrPermutations;
  long Constant = 1l;
   while (NbrUnclusteredParafermions <= nbrParticles)
    {
      GetFixedLzFreeNumberBosonPermutations(NbrUnclusteredParafermions, kType - 1, Permutations, NbrPermutations);
      IntegerPolynomial TotalP;
      for (int i = 0; i < NbrPermutations; ++i)
	{
	  TmpPermutation = Permutations[i];
	  IntegerPolynomial Tmp (0, &Constant, false);
	  bool Flag = true;
	  for (int j = 1; ((j < kType) && (Flag == true)); ++j)
	    {
	      int k = 1;
	      TmpCoefficient = (j * nbrQuasiholes) + ((kType - (2 * (kType - j) * j)) * TmpPermutation[j - 1]);
	      for (; k < j; ++k)
		TmpCoefficient -= (2 * (kType - j) * k) * TmpPermutation[k - 1];
	      ++k;
	      for (; k < kType; ++k)
		TmpCoefficient -= (2 * (kType - k) * j) * TmpPermutation[k - 1];
	      if ((TmpCoefficient >= 0) && ((TmpCoefficient % kType) == 0) && ((TmpPermutation[j - 1] * kType) <= TmpCoefficient))
		{
		  Tmp *= binomialCoeffients(TmpCoefficient / kType, TmpPermutation[j - 1]);
		  cout << TmpPermutation[j - 1] << " " << (TmpCoefficient / kType) << " " << binomialCoeffients(TmpCoefficient / kType, TmpPermutation[j - 1]) << endl;
		}
	      else
		Flag = false;
	    }
	  if (Flag == true)
	    TotalP += Tmp;
	}
      cout << TotalP << endl << endl << "P(1) = " << TotalP.PolynomialEvaluate(1l) << endl;
      for (int i = 0; i < NbrPermutations; ++i)
	delete[] Permutations[i];
      delete[] Permutations;
      NbrUnclusteredParafermions += kType;
    }
}
