////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Cecile Repellin                 //
//                                                                            //
//                                                                            //
//                          class of bosons on CP2                            //
//                                                                            //
//                                                                            //
//                        last modification : 08/01/2013                      //
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


#ifndef BOSONONCP2_H
#define BOSONONCP2_H

#include "config.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include <iostream>



class BosonOnCP2 : public BosonOnSphereShort
{

  friend class FQHESquareLatticeSymmetrizeU1U1StateOperation;
  
 protected:

  // Number of flux quanta
  int NbrFluxQuanta;
  // total value of Tz
  int TotalTz;
  // total value of Y
  int TotalY;
  // total value of r
  int TotalR;
  // total value of s
  int TotalS;
  // array that gives the value of tz for one particle corresponding to the linearized index
  int* quantumNumberTz;
  // array that gives the value of y for one particle corresponding to the linearized index
  int* quantumNumberY;
  // array that gives the value of r for one particle corresponding to the linearized index
  int* quantumNumberR;
  // array that gives the value of s for one particle corresponding to the linearized index
  int* quantumNumberS;
  

 public:

  // default constructor
  // 
  BosonOnCP2 ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta = number of flux quanta (p)
  // totalJz = total value of jz
  // totalKz = total value of kz
  // memory = amount of memory granted for precalculations
  BosonOnCP2 (int nbrBosons, int nbrFluxQuanta, int totalTz, int totalY, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnCP2(const BosonOnCP2& bosons);

  // destructor
  //
  ~BosonOnCP2 ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnCP2& operator = (const BosonOnCP2& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);
  
//   // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given (Jz, Kz) sector.
//   // 
//   // nbrBosonSector = number of particles that belong to the subsytem 
//   // jzSector = Jz sector in which the density matrix has to be evaluated 
//   // kzSector = Kz sector in which the density matrix has to be evaluated 
//   // groundState = reference on the total system ground state
//   // architecture = pointer to the architecture to use parallelized algorithm 
//   // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
//   RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition(int nbrBosonSector, int tzSector, int ysector,  RealVector& groundState, AbstractArchitecture* architecture = 0);
// 
//   // core part of the evaluation density matrix particle partition calculation
//   // 
//   // minIndex = first index to consider in source Hilbert space
//   // nbrIndex = number of indices to consider in source Hilbert space
//   // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
//   // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
//   // groundState = reference on the total system ground state
//   // densityMatrix = reference on the density matrix where result has to stored
//   // return value = number of components that have been added to the density matrix
//   virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
// 								  RealVector& groundState, RealSymmetricMatrix* densityMatrix);
  
  // get the quantum numbers j, jz, kz of a one particle state 
  //
  //quantumNumberJ = array that gives the quantum number j for a single particle state
  //quantumNumberJz = array that gives the quantum number j for a single particle stateDescription
  //quantumNumberKz = array that gives the quantum number j for a single particle state
  inline void GetQuantumNumbersFromLinearizedIndex(int* quantumNumberTz, int* quantumNumberY, int* quantumNumberR, int* quantumNumberS)
  {
    for (int r = 0; r <= this->NbrFluxQuanta; ++r)
	{
	  for (int s = 0; s <= this->NbrFluxQuanta - r ; ++s)
	  {
	    int index = (this->NbrFluxQuanta + 1)*r - (r - 1)*r/2 + s;
	    int t = this->NbrFluxQuanta - r - s;
	    quantumNumberR[index] = r;
	    quantumNumberS[index] = s;
	    quantumNumberTz[index] = r - s;
	    quantumNumberY[index] = r + s - 2*t;
	    
	  }
	}
  }

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentJ = current value of j for a single particle
  // currentJz = current value of jz for a single particle
  // currentKz = current value of kz for a single particle
  // currentTotalJz = current total value of Jz
  // currentTotalKz = current total value of Kz
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentR,  int currentS, int currentTotalR, int currentTotalS);

  // generate all states corresponding to the constraints
  // 
  // stateDescription = array that gives each state description
  // nbrBosons = number of bosons
  // currentJ = current value of j for a single particle
  // currentJz = current value of jz for a single particle
  // currentKz = current value of kz for a single particle
  // currentTotalJz = current total value of Jz
  // currentTotalKz = current total value of Kz
  // currentFermionicPosition = current fermionic position within the state description
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(unsigned long* stateDescription, int nbrBosons, int currentR, int currentS, int currentTotalR, int currentTotalS, int currentFermionicPosition, long pos);

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual bool HasPauliExclusions(int index, int pauliK, int pauliR);
  
};


#endif


