////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of fermions on lattice with spin  and Gutzwiller          //
//   projection in real space with translation invariance in one direction    //
//                          and the Sz<->-Sz symmetry                         //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 30/07/2014                      //
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


#ifndef FERMIONONLATTICEWITHSPINSZSYMMETRYANDGUTZWILLERPROJECTIONREALSPACEAND1DTRANSLATION_H
#define FERMIONONLATTICEWITHSPINSZSYMMETRYANDGUTZWILLERPROJECTIONREALSPACEAND1DTRANSLATION_H


#include "config.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation.h"

#include <iostream>



class FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation : public FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;

 protected:


 public:

  // default constructor
  // 
  FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // momentum = momentum sector
  // periodicity = periodicity with respect to site numbering 
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation (int nbrFermions, int nbrSite, int momentum, int periodicity, bool minusSzParity, 
										      unsigned long memory = 10000000);
  
  // basic constructor when Sz is  conserved
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the total value of Sz
  // nbrSite = number of sites
  // momentum = momentum sector
  // periodicity = periodicity with respect to site numbering 
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation (int nbrFermions, int totalSpin, int nbrSite, int momentum, int periodicity, bool minusSzParity, 
										      unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation(const FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation& operator = (const FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions);
  
  // evaluate Hilbert space dimension when Sz is preserved
  //
  // nbrFermions = number of fermions
  // nbrSpinUp = number of fermions with spin up
  // return value = Hilbert space dimension

  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, long pos);

  // generate all states corresponding to the constraints, knowing the number of holes
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // nbrHoles = number of unoccupied sites
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStatesWithHoleCounting(int nbrFermions, int currentSite, int nbrHoles, long pos);

};


#endif


