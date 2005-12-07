////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class for absorption spectra of quantum well            //
//                                in magnetic field                           //
//                                                                            //
//                        last modification : 06/12/2005                      //
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


#ifndef QUANTUMWELLBFIELDABSORPTIONSPECTRA_H
#define QUANTUMWELLBFIELDABSORPTIONSPECTRA_H

#include "config.h"

#include "Tools/QuantumDot/Spectra/Spectra.h"

class QuantumWellBFieldAbsorptionSpectra : public Spectra
{
 public:

  // constructor from a set of energy files. Each peak is assimilated to a Lorentzian function.
  //
  // fileNumber=  number of files
  // files = name of files
  // stateNumber = integer array containing number of states in each file
  // gamma = lorentzian broadening parameter
  // beta = beta factor (inverse of the temperature in energy unit 1/kT) 
  // eMin = photon minimum energy
  // eMax = photon maximum energy
  // deltaE = photon energy step
  QuantumWellBFieldAbsorptionSpectra(int fileNumber, char** files, int* stateNumber, double gamma, double beta, double eMin, double eMax, double deltaE);

};
