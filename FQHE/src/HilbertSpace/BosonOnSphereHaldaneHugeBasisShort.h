////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on sphere  using the Haldane basis           //
//                            for system size such that                       //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                          work for huge Hilbert space                       //
//                                                                            //
//                        last modification : 08/02/2009                      //
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


#ifndef BOSONONSPHEREHALDANEHUGEBASISSHORT_H
#define BOSONONSPHEREHALDANEHUGEBASISSHORT_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"

#include <iostream>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class BosonOnSphereHaldaneHugeBasisShort :  public BosonOnSphereShort
{

 protected:

  // the fermionic huge Hilbert space associated to the bosonic one
  FermionOnSphereHaldaneHugeBasis* FermionHugeBasis;

 public:

  // default constructor
  //
  BosonOnSphereHaldaneHugeBasisShort ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  // maxFileSize = maximum file size (in MBytes)
  // referenceState = array that describes the reference state to start from
  // memory = amount of memory granted for precalculations
  BosonOnSphereHaldaneHugeBasisShort (int nbrBosons, int totalLz, int lzMax, unsigned long maxFileSize, int* referenceState, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  BosonOnSphereHaldaneHugeBasisShort (char* fileName, unsigned long memory);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereHaldaneHugeBasisShort(const BosonOnSphereHaldaneHugeBasisShort& bosons);

  // destructor
  //
  virtual ~BosonOnSphereHaldaneHugeBasisShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereHaldaneHugeBasisShort& operator = (const BosonOnSphereHaldaneHugeBasisShort& bosons);

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // create the Jack polynomial decomposition corresponding to the root partition
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized hugebasis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized hugebasis
  RealVector& GenerateJackPolynomial(RealVector& jack, double alpha);

  // create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized hugebasis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized hugebasis
  RealVector& GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha);

};

#endif
