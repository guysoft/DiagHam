////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             set of functions used to managed files related to FQHE         //
//                          on square lattice disk                            //
//                                                                            //
//                        last modification : 01/03/2011                      //
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


#ifndef FQHEONSQUARELATTICEFILETOOLS_H
#define FQHEONSQUARELATTICEFILETOOLS_H

#include "config.h"


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnSquareLatticeFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, bool& statistics);

// try to guess system information from file name 
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// mass = reference to the mass term
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, double& mass, bool& statistics);

// try to guess system information from file name for a cubic lattice
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// nbrSiteZ = reference to the number sites along the z direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnCubicLatticeFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& nbrSiteZ, bool& statistics);

// try to guess system information from file name for a cubic lattice
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// nbrSiteZ = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumY = reference to the momentum along the y direction
// momentumZ = reference to the momentum along the z direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnCubicLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& nbrSiteZ, int& momentumX, int& momentumY, int& momentumZ, bool& statistics);


#endif
