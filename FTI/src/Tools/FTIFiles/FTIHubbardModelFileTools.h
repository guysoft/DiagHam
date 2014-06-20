////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     set of functions used to managed files related to Hubbard models       //
//                                                                            //
//                        last modification : 19/06/2014                      //
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


#ifndef FTIHUBBARDMODELFILETOOLS_H
#define FTIHUBBARDMODELFILETOOLS_H

#include "config.h"


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference to the number sites
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// return value = true if no error occured
bool FTIHubbardModelFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSites, bool& statistics);

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference to the number sites
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference to flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, bool& statistics, bool& gutzwiller);


#endif
