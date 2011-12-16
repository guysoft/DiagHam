////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        		      Permutations on array tools   		      //
//                                                                            //
//                        last modification : 08/09/2011                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////



#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include <algorithm>

// evaluate all the different way to split the permutations of a given in two sub groups (not compulsorily the same number of elements in the two subgroups)
//
// nbrPermutations = number of different permutations
// nbrElement = number of element in the permutations
// nbrElementPerColor = number of element in the first group
// permutations1 = array where will be stored the permutations of the first group
// permutations2 = array where will be stored the permutations of the second group

inline void EvaluatePermutationsOfSubGroups(unsigned long nbrPermutations, int nbrElement, int nbrElementPerColor, unsigned long * permutations1, unsigned long * permutations2)
{
	unsigned long NbrElementOtherColor = nbrElement - nbrElementPerColor;
  unsigned long MinValue = (0x1ul << nbrElementPerColor) - 0x1ul;
  unsigned long MaxValue = MinValue << (NbrElementOtherColor);
  unsigned long* TmpArrayPerm = new unsigned long [nbrElement];
  nbrPermutations = 0;
  for (; MinValue <= MaxValue; ++MinValue)
    {
      int Count = 0;
      int Pos = 0;
      while ((Pos < nbrElement) && (Count <= nbrElementPerColor))
	{
	  if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
	    ++Count;
	  ++Pos;
	}
      if (Count == nbrElementPerColor)
	{
	  int Pos1 = 0;
	  int Pos2 = nbrElementPerColor;
	  for (Pos = 0; Pos < nbrElement; ++Pos)
	    {
	      if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
		{
		  TmpArrayPerm[Pos1] = (unsigned long) Pos;
		  ++Pos1;
		}
	      else
		{
		  TmpArrayPerm[Pos2] = (unsigned long) Pos;
		  ++Pos2;
		}
	    }
	  unsigned long TmpPerm2 = 0ul;
	  unsigned long TmpPerm3 = 0ul;
	  for (int i = 0; i < nbrElementPerColor; ++i)
	    {
	      TmpPerm2 |= TmpArrayPerm[i] << (i * 5);
			}
	    for (int i = 0; i < NbrElementOtherColor; ++i)
	    {
				TmpPerm3 |= TmpArrayPerm[i + nbrElementPerColor]<< (i*5);
			}
	  permutations1[nbrPermutations] = TmpPerm2;
	  permutations2[nbrPermutations] = TmpPerm3;	      
	  ++nbrPermutations;
	}
    }
  delete[] TmpArrayPerm;
  return;
}

// evaluate all the different way to split the permutations of a given in two sub groups but exclude ones which are symmetric with higher index (not compulsorily the same number of elements in the two subgroups)
//
// nbrPermutations = number of different permutations
// nbrElement = number of element in the permutations
// nbrElementPerColor = number of element in the first group
// permutations1 = array where will be stored the permutations of the first group
// permutations2 = array where will be stored the permutations of the second group

inline void EvaluatePermutationsOfSubGroupsSymmetric(unsigned long nbrPermutations, int nbrElement, int nbrElementPerColor, unsigned long * permutations1, unsigned long * permutations2)
{
  unsigned long MinValue = (0x1ul << nbrElementPerColor) - 0x1ul;
  unsigned long MaxValue = MinValue << (nbrElement - nbrElementPerColor);
  unsigned long* TmpArrayPerm = new unsigned long [nbrElement];
  unsigned long Mask = (0x1ul << nbrElement) - 1;
  nbrPermutations = 0;
  for (; MinValue <= MaxValue; ++MinValue)
    {
      int Count = 0;
      int Pos = 0;
      while ((Pos < nbrElement) && (Count <= nbrElementPerColor))
	{
	  if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
	    ++Count;
	  ++Pos;
	}
      if (Count == nbrElementPerColor && (MinValue < ((~MinValue) & Mask  )) )
	{
	  int Pos1 = 0;
	  int Pos2 = nbrElementPerColor;
	  for (Pos = 0; Pos < nbrElement; ++Pos)
	    {
	      if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
		{
		  TmpArrayPerm[Pos1] = (unsigned long) Pos;
		  ++Pos1;
		}
	      else
		{
		  TmpArrayPerm[Pos2] = (unsigned long) Pos;
		  ++Pos2;
		}
	    }
	  unsigned long TmpPerm2 = 0ul;
	  unsigned long TmpPerm3 = 0ul;
	  for (int i = 0; i < nbrElementPerColor; ++i)
	    {
	      TmpPerm2 |= TmpArrayPerm[i] << (i * 5);
	      TmpPerm3 |= TmpArrayPerm[i + nbrElementPerColor] << (i *5);
	    }
	  permutations1[nbrPermutations] = TmpPerm2;
	  permutations2[nbrPermutations] = TmpPerm3;	      
	  ++nbrPermutations;
	}
    }
  delete[] TmpArrayPerm;
  return;
}

// evaluate all the different permutations of a given monomial and keep track of the signs
//
// nbrParticles = number of different permutations
// monomial = array storing the monomial
// nbrPermutations = reference to the number of permutations
// monomialPermutatiosn = reference to array where permutations will be stored
// permutationSigns = reference to array where permutation signs will be stored

inline void EvaluateMonomialPermutations(int nbrParticles, unsigned long* monomial, int& nbrPermutations, unsigned long**& monomialPermutations, double*& permutationSigns)
{
  FactorialCoefficient *FactCoef = new FactorialCoefficient();
  FactCoef->FactorialMultiply(nbrParticles);
  monomialPermutations = new unsigned long*[FactCoef->GetIntegerValue()];  
  permutationSigns = new double[FactCoef->GetIntegerValue()];
  delete FactCoef;
  unsigned long* TmpPermutation = new unsigned long[nbrParticles];
  unsigned long* TmpPermutation2 = new unsigned long[nbrParticles];
  unsigned long* TmpPermutation3 = new unsigned long[nbrParticles];
  memcpy(TmpPermutation, monomial, sizeof(unsigned long) * nbrParticles);
  int NumPermutations;
  
  nbrPermutations = 0;
  do 
    {
      monomialPermutations[nbrPermutations] = new unsigned long[nbrParticles];      
      memcpy(monomialPermutations[nbrPermutations], TmpPermutation, sizeof(unsigned long) * nbrParticles);
      
      memcpy(TmpPermutation2, monomial, sizeof(unsigned long) * nbrParticles);
      memcpy(TmpPermutation3, monomialPermutations[nbrPermutations], sizeof(unsigned long) * nbrParticles);
      NumPermutations = 0;
      SortArrayDownOrdering(TmpPermutation3,  TmpPermutation2, nbrParticles, NumPermutations);
      if ( (NumPermutations & 0x1) == 1 )
	{
	  permutationSigns[nbrPermutations] = -1.0;
	}
      else
	{
	  permutationSigns[nbrPermutations] = 1.0;
	}
      nbrPermutations++;
    }
  while (std::prev_permutation(TmpPermutation, TmpPermutation + nbrParticles));
  delete [] TmpPermutation;
  delete [] TmpPermutation2;
  delete [] TmpPermutation3;
}


#endif
