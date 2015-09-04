////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of bosons hardcore on lattice                 //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 11/09/2014                      //
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


#ifndef BOSONONLATTICEGUTZWILLERPROJECTIONREALSPACEAND2DTRANSLATION_H
#define BOSONONLATTICEGUTZWILLERPROJECTIONREALSPACEAND2DTRANSLATION_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"

#include <iostream>



class BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation : public  FermionOnLatticeRealSpaceAnd2DTranslation
{

 protected:

   // number of bosons
   int NbrBosons;	

   // number of bosons + 1
   int IncNbrBosons;

 public:

  // default constructor		
  // 
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation (int nbrBosons, int nbrSite, int xMomentum, int maxXMomentum,
							       int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation (const BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation& bosons);

  // destructor
  //
  ~BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation& operator = (const BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // apply a^+_m_u a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);


 protected:

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(unsigned long& state, double& coefficient, 
				   int& nbrTranslationX, int& nbrTranslationY);

  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalFormAndTestMomentumConstraint(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);


};


// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::SymmetrizeAdAdResult(unsigned long& state, double& coefficient, 
										   int& nbrTranslationX, int& nbrTranslationY)
{
  this->FindCanonicalFormAndTestMomentumConstraint(state, nbrTranslationX, nbrTranslationY);
  if (nbrTranslationX < 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpMaxMomentum = this->NbrSite;
  while ((state >> TmpMaxMomentum) == 0x0ul)
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
    }
  return TmpIndex;
}



// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::FindCanonicalFormAndTestMomentumConstraint(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;  
  unsigned long TmpStateDescription;  
  nbrTranslationX = 0;
  nbrTranslationY = 0;
   for (int m = 0; (m < this->MaxYMomentum) && (stateDescriptionReference != stateDescription) ; ++m)
    {
      TmpStateDescription = stateDescription;
      for (int n = 0; (n < this->MaxXMomentum) && (TmpStateDescription != stateDescription) ; ++n)
	{
	  if (TmpStateDescription < CanonicalState)
	    {
	      CanonicalState = TmpStateDescription;
	      nbrTranslationX = n;	      
	      nbrTranslationY = m;	      
	    }
	  this->ApplySingleXTranslation(TmpStateDescription);      
	}
      this->ApplySingleYTranslation(stateDescription);
    }
  return CanonicalState;  
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::TestMomentumConstraint(unsigned long stateDescription)
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
  if (((this->XMomentum * XSize) % this->MaxXMomentum) != 0)
    return false;
  int YSize = this->MaxYMomentum;
  int TmpXSize = 0;
  TmpStateDescription2 = stateDescription;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(TmpStateDescription2);	
//      cout << hex << stateDescription << " " << TmpStateDescription2 << " " << dec <<endl;
      TmpStateDescription = TmpStateDescription2;
      TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);
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
    
  if ((((this->YMomentum * YSize * this->MaxXMomentum)
	+ (this->XMomentum * TmpXSize * this->MaxYMomentum)) % (this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;
  return true;
}


#endif


