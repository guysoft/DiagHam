////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          set of functions used to handle torus pseudopotential files       //
//                                                                            //
//                        last modification : 17/04/2012                      //
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


#ifndef FQHETORUSPSEUDOPOTENTIALTOOLS_H
#define FQHETORUSPSEUDOPOTENTIALTOOLS_H

#include "config.h"



// get pseudopototentials for particles on torus with SU(2) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// return value = true if no error occured
bool FQHETorusSU2GetPseudopotentials (char* fileName, int* nbrPseudoPotentials, double** pseudoPotentials);


#endif

