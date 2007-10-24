////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    set of functions used to managed files related to QHE on sphere         //
//                                                                            //
//                        last modification : 06/06/2006                      //
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


#ifndef QHEONSPHEREFILETOOLS_H
#define QHEONSPHEREFILETOOLS_H

#include "config.h"


// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured
bool QHEOnSphereFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& lzMax, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured
bool QHEOnSphereFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, bool& statistics);

// try to guess system information from file name for system suth an SU(2) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// sz = reference to twice the z projection of the total spin (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured
bool QHEOnSphereWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, int& sz, bool& statistics);

// try to guess system information from file name for system suth an SU(4) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// sz = reference to twice the z projection of the total spin (grab it only if initial value is 0)
// iz = reference to twice the z projection of the total isospin (grab it only if initial value is 0)
// ez = reference to twice the entanglement (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnSphereWithSU4SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, int& sz, int& iz, int& ez, bool& statistics);


#endif
