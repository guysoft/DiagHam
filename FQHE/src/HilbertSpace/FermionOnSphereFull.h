////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//                                                                            //
//                        last modification : 06/07/2006                      //
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


#ifndef FERMIONONSPHEREFULL_H
#define FERMIONONSPHEREFULL_H


#include "config.h"
#include "HilbertSpace/FermionOnSphere.h"

#include <iostream>



class FermionOnSphereFull :  public FermionOnSphere
{

 protected:

 public:

  // default constructor
  //
  FermionOnSphereFull();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = reference on twice the momentum total value, totalLz will be recomputed from referenceState and stored in totalLz
  // lzMax = maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  // referenceState = array that describes the reference state to start from
  FermionOnSphereFull (int nbrFermions, int lzMax);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereFull(const FermionOnSphereFull& fermions);

  // destructor
  //
  virtual ~FermionOnSphereFull ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereFull& operator = (const FermionOnSphereFull& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // generate all states (i.e. all possible skew symmetric polynomials)
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int lzMax, int currentLzMax, long pos);

};


#endif


