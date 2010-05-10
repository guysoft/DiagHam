////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                  class of fermions on sphere including three               //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 20/04/2010                      //
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


#ifndef FERMIONONSPHERETHREELANDAULEVELS_H
#define FERMIONONSPHERETHREELANDAULEVELS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include <iostream>


class FermionOnSphere;


class FermionOnSphereThreeLandauLevels : public FermionOnSphereWithSU3Spin
{

 protected:


 public:

  // default constructor
  // 
  FermionOnSphereThreeLandauLevels ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  FermionOnSphereThreeLandauLevels (int nbrFermions, int totalLz, int lzMax, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereThreeLandauLevels(const FermionOnSphereThreeLandauLevels& fermions);

  // destructor
  //
  virtual ~FermionOnSphereThreeLandauLevels ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereThreeLandauLevels& operator = (const FermionOnSphereThreeLandauLevels& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

   // compute the projection of the product of a bosonic state and a fermionic state
  //
  // bosonState = real vector where the bosonic state is stored
  // fermionState = real vector where the fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // finalStates = array where the obtained states are stored in their fermionic representation
  // weigth = array where the coefficients for each obtained state are stored
  // bosonSpace = pointer to the bosonic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  virtual void BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, unsigned long* finalStates, double* weigth, 
					      BosonOnSphereShort* bosonSpace,FermionOnSphere * finalSpace, int firstComponent,int nbrComponent);

  // compute the product and the projection of a Slater determinant and a monomial 
  // 
  // slater = array where the slater is stored in its monomial representation
  // monomial = array where the monomial is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // finalStates = array where the obtained states are stored in their fermionic representation
  // weigth = array where the coefficients for each obtained state are stored
  // finalSpace = pointer to the final HilbertSpace
  // return value = number of different obtained states
  virtual unsigned int MonomialsTimesSlaterProjection(unsigned long* slater, unsigned long* monomial, unsigned long* variable, int nbrVariable, unsigned long*& finalStates, 
						      double*& weigth, FermionOnSphere* finalSpace);
  
 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // nbrFluxQuanta = number of flux quanta
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // nbrFluxQuanta = number of flux quanta
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz, long pos);

  // convert a fermionic state to its monomial variable representation
  //
  // state = initial fermionic state in its fermionic representation
  // slater = reference on the array where the monomial representation has to be stored
  // nbrVariable = number of fermions in the second Landau level
  // variable = reference on the array where the indice of fermions in the second Landau and the third level has to be stored
  virtual void ConvertToMonomialVariable(unsigned long state, unsigned long*& slater, int& nbrVariable, unsigned long*& variable);

};

// convert a fermionique state to its monomial variable representation
//
// state = initial fermionic state in its fermionic representation
// slater = reference on the array where the monomial representation has to be stored
// nbrVariable = number of fermions in the second Landau level
// variable = reference on the array where the indice of fermions in the second Landau and the third level has to be stored

inline void FermionOnSphereThreeLandauLevels::ConvertToMonomialVariable(unsigned long state, unsigned long*& slater, int& nbrVariable, unsigned long*& variable)
{
  unsigned long Tmp;
  int TmpPos=0;
  for(int i=this->LzMax;i>=0;i--)
    {
      Tmp=((state >> (3*i)) & ((unsigned long) 0x7));
      switch (Tmp)
	{
	case 0x1l: 
	  {
	    slater[TmpPos]=i;
	    TmpPos++;
	    break;
	  }
	case 0x2l:
	  {
	    slater[TmpPos]=i;
	    variable[nbrVariable]=2*TmpPos;
	    nbrVariable++;
	    TmpPos++;
	    break;
	  }
	case 0x3l:
	  {
	    slater[TmpPos]=i;
	    variable[nbrVariable]=2*TmpPos;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos]=i;
	    TmpPos++;
	    break;
	  }
	case 0x4l: 
	  {
	    slater[TmpPos]=i;
	    variable[nbrVariable]=2*TmpPos+0x01ul;
	    nbrVariable++;
	    TmpPos++;
	    break;
	  }
	case 0x5l:
	  {
	    slater[TmpPos]=i;
	    variable[nbrVariable]=2*TmpPos+0x01ul;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos]=i;
	    TmpPos++;
	    break;
	  }
	case 0x6l:
	  {
	    slater[TmpPos]=i;
	    variable[nbrVariable]=2*TmpPos+0x01ul;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos]=i;
	    variable[nbrVariable]=2*TmpPos;
	    nbrVariable++;
	    TmpPos++;
	    break;
	  }
	case 0x7l:
	  {
	    slater[TmpPos]=i;
	    variable[nbrVariable]=2*TmpPos+0x01ul;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos]=i;
	    variable[nbrVariable]=2*TmpPos;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos]=i;
	    TmpPos++;
	    break;
	  }
	  break;
	default : 
	  {
	    break;
	  }
	}
    }
}


#endif


