////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of fermions on square lattice                 //
//                                  in momentum space                         //
//                                                                            //
//                        last modification : 16/02/2011                      //
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


#ifndef FERMIONONSQUARELATTICEMOMENTUMSPACE_H
#define FERMIONONSQUARELATTICEMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSphere.h"

#include <iostream>



class FermionOnSquareLatticeMomentumSpace : public FermionOnSphere
{

 protected:

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;
  // momentum along the x direction
  int KxMomentum;
  // momentum along the y direction
  int KyMomentum;

 public:

  // default constructor
  // 
  FermionOnSquareLatticeMomentumSpace ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeMomentumSpace(const FermionOnSquareLatticeMomentumSpace& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeMomentumSpace& operator = (const FermionOnSquareLatticeMomentumSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given momentum sector and fixed number of particles
  // 
  // subsytemSizeX = number of sites along the x direction that belong to the subsytem
  // subsytemSizeY = number of sites along the y direction that belong to the subsytem
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
  // kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem
  virtual HermitianMatrix EvaluatePartialDensityMatrixMomentumSpace (int subsytemSizeX, int subsytemSizeY, int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given momentum sector and fixed number of particles
  // 
  // subsytemSizeX = number of sites along the x direction that belong to the subsytem
  // subsytemSizeY = number of sites along the y direction that belong to the subsytem
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
  // kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrixMomentumSpace (int subsytemSizeX, int subsytemSizeY, int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
  // kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, bool removeBinomialCoefficient = false);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos);

};


#endif


