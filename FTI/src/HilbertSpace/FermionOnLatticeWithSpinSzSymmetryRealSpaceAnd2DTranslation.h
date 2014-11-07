////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//       in real space with translation invariance in two directions and      //
//                               Sz<->-Sz symmetry                            //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 03/11/2014                      //
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


#ifndef FERMIONONLATTICEWITHSPINSZSYMMETRYREALSPACEAND2DTRANSLATION_H
#define FERMIONONLATTICEWITHSPINSZSYMMETRYREALSPACEAND2DTRANSLATION_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"

#include <iostream>


#define FERMION_LATTICE_REALSPACE_SU2_SZ_MASK 0x5555555555555555ul


class FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation : public FermionOnLatticeWithSpinRealSpaceAnd2DTranslation
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;

 protected:

  // sign of the parity sector for the Sz<->-Sz symmetry
  double SzParitySign;

  // sign of the parity sector for the Sz<->-Sz symmetry, 0 for 1, 1 for -1
  unsigned long SzParity;

 public:

  // default constructor
  // 
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (int nbrFermions, int nbrSite, bool minusSzParity, int xMomentum, int maxXMomentum,
							       int yMomentum, int maxYMomentum, unsigned long memory = 10000000);
  
  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the value of Sz
  // nbrSite = number of sites
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (int nbrFermions, int totalSpin, int nbrSite, bool minusSzParity, int xMomentum, int maxXMomentum,
							       int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation(const FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation& operator = (const FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();
  
 protected:

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(unsigned long& state, double& coefficient, 
				   int& nbrTranslationX, int& nbrTranslationY);

  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = reference on the number of translations to apply in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to apply in the y direction to the resulting state to obtain the return orbit describing state
  // nbrSpinFlip = reference on the number of full spin flip to apply to the resulting state to obtain the return orbit describing state
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY, int& nbrSpinFlip);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(unsigned long stateDescription);

  // Apply the Sz operator to flip all the spins
  //
  // stateDescription = reference on state description
  // stateDescription = state that has to be converted to its canonical expression
  virtual void ApplySzSymmetry (unsigned long& stateDescription);
  
  // get the fermonic sign when performing a flip all the spins, and apply the flip sign 
  //
  // stateDescription = reference on state description
  // return value = 0 if the sign is +1, 1 if the sign is -1
  virtual unsigned long GetSignAndApplySzSymmetry (unsigned long& stateDescription);

};


// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state  

inline int FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::SymmetrizeAdAdResult(unsigned long& state, double& coefficient, 
											     int& nbrTranslationX, int& nbrTranslationY)
{
  int NbrSpinFlip;
  state = this->FindCanonicalForm(state, nbrTranslationX, nbrTranslationY, NbrSpinFlip);
  int TmpMaxMomentum = 2 * this->NbrSite + 1;
  while ((state >> TmpMaxMomentum) == 0x0ul)
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslationX = (this->MaxXMomentum - nbrTranslationX) % this->MaxXMomentum;
      nbrTranslationY = (this->MaxYMomentum - nbrTranslationY) % this->MaxYMomentum;
      coefficient *= 1.0 - (2.0 * ((double) (((this->ReorderingSign[TmpIndex] >> ((((nbrTranslationY * this->MaxXMomentum) + nbrTranslationX) * 2) + NbrSpinFlip)) & 0x1u)))); 
    }
  return TmpIndex;
}


// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations toapply  in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to apply in the y direction to the resulting state to obtain the return orbit describing state
// nbrSpinFlip = reference on the number of full spin flip to apply to the resulting state to obtain the return orbit describing state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY,
												    int& nbrSpinFlip)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;  
  unsigned long TmpStateDescription;  
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  nbrSpinFlip = 0;
  TmpStateDescription = stateDescription;
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslationX = n;	      
	  nbrTranslationY = 0;	      
	}
    }
  for (int m = 1; m < this->MaxYMomentum; ++m)
    {
      this->ApplySingleYTranslation(stateDescription);      
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslationX = 0;	      
	  nbrTranslationY = m;	      
	}
      TmpStateDescription = stateDescription;
      for (int n = 1; n < this->MaxXMomentum; ++n)
	{
	  this->ApplySingleXTranslation(TmpStateDescription);      
	  if (TmpStateDescription < CanonicalState)
	    {
	      CanonicalState = TmpStateDescription;
	      nbrTranslationX = n;	      
	      nbrTranslationY = m;	      
	    }
	}
    }
  stateDescription = stateDescriptionReference;
  this->ApplySzSymmetry(stateDescription);
  TmpStateDescription = stateDescription;
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      nbrTranslationX = 0;	      
      nbrTranslationY = 0;	      
      nbrSpinFlip = 1;
    } 
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslationX = n;	      
	  nbrTranslationY = 0;	      
	  nbrSpinFlip = 1;
	}
    }
  for (int m = 1; m < this->MaxYMomentum; ++m)
    {
      this->ApplySingleYTranslation(stateDescription);      
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslationX = 0;	      
	  nbrTranslationY = m;	      
	  nbrSpinFlip = 1;
	}
      TmpStateDescription = stateDescription;
      for (int n = 1; n < this->MaxXMomentum; ++n)
	{
	  this->ApplySingleXTranslation(TmpStateDescription);      
	  if (TmpStateDescription < CanonicalState)
	    {
	      CanonicalState = TmpStateDescription;
	      nbrTranslationX = n;	      
	      nbrTranslationY = m;	      
	      nbrSpinFlip = 1;
	    }
	}
    }  
  return CanonicalState;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::TestMomentumConstraint(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  unsigned long TmpStateDescription2 = stateDescription;
  unsigned long TmpStateDescription3 = stateDescription;
  int XSize = 1;
  unsigned long TmpSign = this->GetSignAndApplySingleXTranslation(TmpStateDescription);   
  unsigned long TmpSign2 = 0x0ul;
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      TmpSign ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
    }
  if ((((this->XMomentum * XSize) + ((((int) TmpSign) * this->MaxXMomentum) >> 1)) % this->MaxXMomentum) != 0)
    return false;
  int YSize = this->MaxYMomentum;
  int TmpXSize = 0;
  TmpSign = 0x0ul;
  TmpStateDescription2 = stateDescription;
  for (int m = 1; m < YSize; ++m)
    {
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
      TmpSign2 = TmpSign;
      TmpStateDescription = TmpStateDescription2;
      TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  TmpSign2 ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
      else
	{
	  TmpXSize = 0;
	}
    } 
  if (YSize == this->MaxYMomentum)
    {
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
      TmpSign2 = TmpSign;
    }
  if ((((2 * this->YMomentum * YSize * this->MaxXMomentum)
	+ (2 * this->XMomentum * TmpXSize * this->MaxYMomentum)
	+ (((int) TmpSign2) * this->MaxXMomentum * this->MaxYMomentum)) % (2 * this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;

  stateDescription = TmpStateDescription3;
  TmpStateDescription = stateDescription;
  TmpSign = this->GetSignAndApplySzSymmetry(TmpStateDescription) ^ this->SzParity;
  if (stateDescription == TmpStateDescription)
    {
      if (TmpSign != 0x0ul)
	{
	  return false;
	}
      else
	{
	  return true;
	}
    }
  stateDescription = TmpStateDescription;
  TmpStateDescription = stateDescription;
  TmpStateDescription2 = TmpStateDescription3;
  for (int m = 1; m < YSize; ++m)
    {
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
      TmpSign2 = TmpSign;
      TmpStateDescription = TmpStateDescription2;
      TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  TmpSign2 ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
      else
	{
	  TmpXSize = 0;
	}
    } 
  if (YSize == this->MaxYMomentum)
    {
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
      TmpSign2 = TmpSign;
    }
  if ((((2 * this->YMomentum * YSize * this->MaxXMomentum)
	+ (2 * this->XMomentum * TmpXSize * this->MaxYMomentum)
	+ (((int) TmpSign2) * this->MaxXMomentum * this->MaxYMomentum)) % (2 * this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;
  return true;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::FindOrbitSize(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  unsigned long TmpStateDescription2 = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  int YSize = this->MaxYMomentum;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(stateDescription); 
      TmpStateDescription = TmpStateDescription2;
      int TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
    }
  
  unsigned long TmpStateDescription3 = TmpStateDescription2;
  this->ApplySzSymmetry(TmpStateDescription3);
  if (TmpStateDescription3 == TmpStateDescription2)
    return (XSize * YSize);  

  XSize = 1;
  TmpStateDescription = TmpStateDescription2;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((TmpStateDescription3 != TmpStateDescription) && (TmpStateDescription2 != TmpStateDescription))
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  
  stateDescription = TmpStateDescription3;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(stateDescription); 
      TmpStateDescription = TmpStateDescription2;
      int TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
    }
  return (2 * XSize * YSize);
}

// Apply the Sz operator to flip all the spins
//
// stateDescription = reference on state description

inline void FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::ApplySzSymmetry (unsigned long& stateDescription)
{
  
  unsigned long TmpState = stateDescription;
  stateDescription = ((TmpState >> 1) ^ TmpState) & FERMION_LATTICE_REALSPACE_SU2_SZ_MASK;
  stateDescription |= stateDescription << 1;
  stateDescription ^= TmpState; 
}

// get the fermonic sign when performing a flip all the spins, and apply the flip sign 
//
// stateDescription = reference on state description
// return value = 0 if the sign is +1, 1 if the sign is -1

inline unsigned long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::GetSignAndApplySzSymmetry (unsigned long& stateDescription)
{
  unsigned long TmpState = stateDescription;
  stateDescription = ((TmpState >> 1) ^ TmpState) & FERMION_LATTICE_REALSPACE_SU2_SZ_MASK;
  stateDescription |= stateDescription << 1;
  stateDescription ^= TmpState; 
  // compute the parity of the pair number
  TmpState &= (TmpState >> 1);
  TmpState &= FERMION_LATTICE_REALSPACE_SU2_SZ_MASK;
#ifdef __64_BITS__
  TmpState ^= (TmpState >> 32);
#endif
  TmpState ^= (TmpState >> 16);
  TmpState ^= (TmpState >> 8);
  TmpState ^= (TmpState >> 4);
  TmpState ^= (TmpState >> 2);
  return (TmpState & 0x01ul);
}

#endif


