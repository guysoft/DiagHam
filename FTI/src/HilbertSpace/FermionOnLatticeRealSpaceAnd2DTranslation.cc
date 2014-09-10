////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 20/08/2014                      //
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


#include "config.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

FermionOnLatticeRealSpaceAnd2DTranslation::FermionOnLatticeRealSpaceAnd2DTranslation ()
{
  this->NbrFermions = 0;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = 0;
  this->NbrSite = 0;
  this->MomentumModulo = 0;
  this->XMomentum = 0; 
  this->YMomentum = 0;
  this->StateShift = 2 * (this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateMaxMomentum = 0;  
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// xMomentum = momentum sector in the x direction
// xTranslation = translation that has to be applied on the site index to connect two sites with a translation in the x direction
// yMomentum = momentum sector in the y direction
// yPeriodicity = periodicity in the y direction with respect to site numbering 
// memory = amount of memory granted for precalculations

FermionOnLatticeRealSpaceAnd2DTranslation::FermionOnLatticeRealSpaceAnd2DTranslation (int nbrFermions, int nbrSite, int xMomentum, int xTranslation,
										      int yMomentum, int yPeriodicity, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MomentumModulo = this->NbrSite / xTranslation;
  cout << "MomentumModulo=" << MomentumModulo<<endl;

  this->StateXShift = 1;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->MaxXMomentum = this->NbrSite / xTranslation;
  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = xTranslation;
  this->ComplementaryStateXShift = this->MaxMomentum - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;
  cout << "this->MaxXMomentum=" << this->MaxXMomentum << endl;
  cout << "this->XMomentum=" << this->XMomentum << endl;
  cout << "this->StateXShift=" << this->StateXShift << endl;
  cout << "this->ComplementaryStateXShift=" << this->ComplementaryStateXShift << endl;
  cout << "this->XMomentumMask=" << hex << this->XMomentumMask << dec << endl;

  this->MaxYMomentum = yPeriodicity;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = this->NbrSite / this->MaxXMomentum;
  this->StateYShift = (this->NbrYMomentumBlocks / this->MaxYMomentum);
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  cout << "this->NbrYMomentumBlocks=" << this->NbrYMomentumBlocks << endl;
  cout << "this->StateYShift=" << this->StateYShift << endl;
  cout << "this->YMomentumBlockSize=" << this->YMomentumBlockSize << endl;
  cout << "this->ComplementaryStateYShift=" << this->ComplementaryStateYShift << endl;
  cout << "this->YMomentumMask=" << hex << this->YMomentumMask << dec << endl;
  cout << "this->YMomentumBlockMask=" << hex << this->YMomentumBlockMask << dec << endl;
  this->NbrFermionsParity = (~((unsigned long) this->NbrFermions)) & 0x1ul;


  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->GenerateSignLookUpTable();
      this->LargeHilbertSpaceDimension  = this->GenerateStates();
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->HilbertSpaceDimension > 0)
	{
	  this->StateMaxMomentum = new int [this->LargeHilbertSpaceDimension];  
	  int CurrentMaxMomentum = this->MaxMomentum;
	  while (((this->StateDescription[0] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
	    --CurrentMaxMomentum;
	  this->StateMaxMomentum[0] = CurrentMaxMomentum;
	  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      while (((this->StateDescription[i] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
		--CurrentMaxMomentum;
	      this->StateMaxMomentum[i] = CurrentMaxMomentum;
	    }
	  this->GenerateLookUpTable(memory);
	}
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->NbrMomentum * sizeof(int);
      UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeRealSpaceAnd2DTranslation::FermionOnLatticeRealSpaceAnd2DTranslation(const FermionOnLatticeRealSpaceAnd2DTranslation& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->NbrSite = fermions.NbrSite;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnLatticeRealSpaceAnd2DTranslation::~FermionOnLatticeRealSpaceAnd2DTranslation ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeRealSpaceAnd2DTranslation& FermionOnLatticeRealSpaceAnd2DTranslation::operator = (const FermionOnLatticeRealSpaceAnd2DTranslation& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->NbrSite = fermions.NbrSite;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnLatticeRealSpaceAnd2DTranslation::Clone()
{
  return new FermionOnLatticeRealSpaceAnd2DTranslation(*this);
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long FermionOnLatticeRealSpaceAnd2DTranslation::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(this->NbrFermions, this->NbrSite - 1, 0l);
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY) == this->StateDescription[i]))
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = 0x0ul;
	    }
	}
      else
	{
	  this->StateDescription[i] = 0x0ul;
	}
    }
  cout << "new dim = " << TmpLargeHilbertSpaceDimension << endl;
  unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  this->ReorderingSign = new unsigned long [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != 0x0ul)
	{
	  TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	  this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = this->FindOrbitSize(this->StateDescription[i]);
	  unsigned long& TmpSign = this->ReorderingSign[TmpLargeHilbertSpaceDimension];
	  TmpSign = 0x0ul;	  
	  int Index = 1;
	  unsigned long TmpState =  this->StateDescription[i];
	  for (int m = 0; m < this->MaxYMomentum; ++m)
	    {
	      unsigned long TmpState2 = TmpState;
	      for (int n = 1; n < this->MaxXMomentum; ++n)
		{
		  TmpSign |= (this->GetSignAndApplySingleXTranslation(TmpState2) << Index) ^ ((TmpSign & (0x1ul << (Index - 1))) << 1);
		  ++Index;
		}
	      TmpSign |= this->GetSignAndApplySingleYTranslation(TmpState) << Index;
	      ++Index;
	    }
	  ++TmpLargeHilbertSpaceDimension;
	}
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  return TmpLargeHilbertSpaceDimension;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeRealSpaceAnd2DTranslation::RawGenerateStates(int nbrFermions, int currentSite, long pos)
{
  if (nbrFermions == 0)
    {
      this->StateDescription[pos] = 0x0ul;	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentSite; j >= 0; --j)
	{
	  this->StateDescription[pos] = 0x1ul << j;
	  ++pos;
	}
      return pos;
    }
  long TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, pos);
  unsigned long Mask = 0x1ul << currentSite;
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->RawGenerateStates(nbrFermions, currentSite - 1, pos);
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// return value = Hilbert space dimension

long FermionOnLatticeRealSpaceAnd2DTranslation::EvaluateHilbertSpaceDimension(int nbrFermions)
{
  BinomialCoefficients binomials(this->NbrSite);
  long dimension = binomials(this->NbrSite, this->NbrFermions);
  return dimension;
}

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int FermionOnLatticeRealSpaceAnd2DTranslation::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != this->MaxMomentum)
    return -1;
  unsigned long TmpState = 0x0ul;
  int TmpNbrParticles = 0;
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      if (TmpDescription[i][0] == '1')
	{
	  TmpState |= 0x1ul << i;
	  ++TmpNbrParticles;	  
	}
      else
	{
	  if (TmpDescription[i][0] != '0')
	    {
	      return -1;
	    }
	}
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if (TmpNbrParticles != this->NbrFermions)
    return -1;
  int TmpLzMax = this->MaxMomentum;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}

// find state index
//
// stateDescription = unsigned longeger describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnLatticeRealSpaceAnd2DTranslation::FindStateIndex(unsigned long stateDescription, int maxMomentum)
{
  if ((stateDescription < this->StateDescription[0]) || (stateDescription > this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;
  long PosMax = stateDescription >> this->LookUpTableShift[maxMomentum];
  long PosMin = this->LookUpTable[maxMomentum][PosMax];
  PosMax = this->LookUpTable[maxMomentum][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMin != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMax;
}

// apply a^+_m_sigma a_n_sigma operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator including the orbital and the spin index
// n = index of the annihilation operator including the orbital and the spin index
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int FermionOnLatticeRealSpaceAnd2DTranslation::AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long State = this->StateDescription[index];
  if ((State & (0x1ul << n)) == 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if ((State & (0x1ul << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
  coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
  coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
  State |= (0x1ul << m);
  return this->SymmetrizeAdAdResult(State, coefficient, nbrTranslationX, nbrTranslationY);
}
  
// apply a^+_m1_sigma a^+_m2_sigma operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int FermionOnLatticeRealSpaceAnd2DTranslation::AdAd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
  TmpState |= (0x1ul << m2);
  
  coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}
