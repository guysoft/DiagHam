////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     set of functions used to managed files related to spin system          //
//                                                                            //
//                        last modification : 22/06/2010                      //
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


#ifndef SPINFILETOOLS_H
#define SPINFILETOOLS_H

#include "config.h"


// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins 
// return value = true if no error occured
bool SpinFindSystemInfoFromFileName(char* filename, int& nbrSpins);

// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins 
// spin = reference to twice the spin value per site
// return value = true if no error occured
bool SpinFindSystemInfoFromFileName(char* filename, int& nbrSpins, int& spin);

// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// return value = true if no error occured
bool SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin);

// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// momentum = reference on the momentum
// return value = true if no error occured
bool SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& momentum);

// try to guess system information from file name for a 2d spin system with translations
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// xMomentum = reference on the momentum along the x direction
// xPeriodicity = reference on the number of sites along the x direction
// yMomentum = reference on the momentum along the y direction
// yPeriodicity = reference on the number of sites along the y direction
// return value = true if no error occured
bool SpinWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& xMomentum, int& xPeriodicity,
							   int& yMomentum, int& yPeriodicity);

// try to guess system information from file name for a 2d spin system with translations
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// xMomentum = reference on the momentum along the x direction
// xPeriodicity = reference on the number of sites along the x direction
// yMomentum = reference on the momentum along the y direction
// yPeriodicity = reference on the number of sites along the y direction
// inversion =  reference on the inversion parity
// return value = true if no error occured
bool SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& xMomentum, int& xPeriodicity,
							   int& yMomentum, int& yPeriodicity, int& inversion);

// try to guess system information from file name for a 2d spin system with translations
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// xMomentum = reference on the momentum along the x direction
// xPeriodicity = reference on the number of sites along the x direction
// yMomentum = reference on the momentum along the y direction
// yPeriodicity = reference on the number of sites along the y direction
// return value = true if no error occured
bool SpinWith2DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& spin, int& xMomentum, int& xPeriodicity,
							   int& yMomentum, int& yPeriodicity);

// try to guess system information from file name for a 2d spin system with translations
//
// filename = file name
// nbrSpins = reference to the number of spins 
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// xMomentum = reference on the momentum along the x direction
// xPeriodicity = reference on the number of sites along the x direction
// yMomentum = reference on the momentum along the y direction
// yPeriodicity = reference on the number of sites along the y direction
// inversion =  reference on the inversion parity
// return value = true if no error occured
bool SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& spin, int& xMomentum, int& xPeriodicity,
							   int& yMomentum, int& yPeriodicity, int& inversion);

// try to guess system information from file name
//
// filename = file name
// nbrSites = reference to the number of sites
// qValue = reference to Zn charge
// return value = true if no error occured
bool PottsFindSystemInfoFromVectorFileName(char* filename, int& nbrSites, int& qValue);

#endif
