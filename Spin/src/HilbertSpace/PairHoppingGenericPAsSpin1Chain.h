////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of pair hopping with generic p                  //
//                     Hilbert space written as spin 1 chain                  //
//                             for more than 32 spins                         //
//                                                                            //
//                        last modification : 17/10/2019                      //
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


#ifndef PAIRHOPPINGGENERICPSPIN1CHAINLONG_H
#define PAIRHOPPINGGENERICPSPIN1CHAINLONG_H


#include "config.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainLong.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class ComplexMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class PairHoppingGenericPAsSpin1ChainLong : public  PairHoppingP2AsSpin1ChainLong
{

  friend class Spin1ChainWithTranslationsLong;
  friend class PairHoppingP1AsSpin1ChainWithTranslationsLong;
  friend class PairHoppingP2AsSpin1ChainWithTranslationsLong;

 protected:

  // number of spin 1 per unit cell
  int UnitCellSize;

  // mask to isolate the first unit cell
  ULONGLONG FirstUnitCellMask;

  // number of possible unit cells satisfying the pair hopping constraints
  int NbrPossibleUnitCells;
  // possible unit cells satisfying the pair hopping constraints
  ULONGLONG* PossibleUnitCells;
  // number of plusses for each possible unit cells satisfying the pair hopping constraints
  int* PossibleUnitCellsNbrPlus;
  
  // number of possible unit cells satisfying the pair hopping constraints for each possible number of minuses
  int* NbrPossibleUnitCellsPerNbrMinus;
  // possible unit cells satisfying the pair hopping constraints for each possible number of minuses
  ULONGLONG** PossibleUnitCellsPerNbrMinus;
  // number of plusses for each possible unit cells satisfying the pair hopping constraints for each possible number of minuses
  int** PossibleUnitCellsNbrPlusPerNbrMinus;
  
public:

  // default constructor
  //
  PairHoppingGenericPAsSpin1ChainLong ();

  // constructor for complete Hilbert space
  //
  // chainLength = number of spin / group of 2p+1 orbitals
  // pValue = p value
  // periodicBoundaryConditions = true if the system uses periodic boundary conditions
  // memorySize = memory size in bytes allowed for look-up table
  // useForEntanglementMatrix = true if the hilbert space has to be generated for the entanglement matrix calculations
  PairHoppingGenericPAsSpin1ChainLong (int chainLength, int pValue, bool periodicBoundaryConditions = false, int memorySize = -1, bool useForEntanglementMatrix = false);

  // constructor for the Hilbert space, fixing the number of pluses in the rightmost unit cell and the number of minuses in the leftmost unit cell
  //
  // chainLength = number of spin / group of 2p+1 orbitals
  // nbrPlusRight = number of pluses in the rightmost unit cell
  // nbrMinusLeft = number of minuses in the lefttmost unit cell
  PairHoppingGenericPAsSpin1ChainLong (int chainLength, int nbrPlusRight, int nbrMinusLeft, int memorySize = -1);
  
  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  PairHoppingGenericPAsSpin1ChainLong (const PairHoppingGenericPAsSpin1ChainLong& chain);

  // destructor
  //
  virtual ~PairHoppingGenericPAsSpin1ChainLong ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  PairHoppingGenericPAsSpin1ChainLong& operator = (const PairHoppingGenericPAsSpin1ChainLong& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // apply the swap operator within the unit cell
  //
  // unitCellCoordinate = coordinate of the unit cell
  // siteIndex = index of the leftmost site within the unit cell
  // state = index of the state on which the operator has to be applied
  // return value = index of resulting state 
  virtual int SwapOperator (int unitCellCoordinate, int siteIndex, int state);

  // apply the operator coupling neighboring unit cells
  //
  // leftUnitCellCoordinate = coordinate of the left unit cell
  // rightUnitCellCoordinate = coordinate of the right unit cell
  // state = index of the state on which the operator has to be applied
  // return value = index of resulting state 
  virtual int PlusMinusOperator (int leftUnitCellCoordinate, int rightUnitCellCoordinate, int state);  
  
 protected:

  // evaluate Hilbert space dimension
  //
  // previousSpin = value of the previous spin (0 for -1, 1 for 0 and 2 for +1)
  // sitePosition = site on chain where spin has to be changed
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int previousSpin, int sitePosition);

  // generate all states with neither constraint from boundary conditions nor discrete symmtry constraint
  //
  // statePosition = position for the new states
  // previousSpin = value of the previous spin (0 for -1, 1 for 0 and 2 for +1)
  // sitePosition = site on chain where spin has to be changed
  // return value = number of generated states
  virtual long RawGenerateStates(long statePosition, int previousSpin, int sitePosition);

  // generate all states with constraints 
  //
  // return value = number of generated states
  virtual long GenerateStates();
  
};

#endif


