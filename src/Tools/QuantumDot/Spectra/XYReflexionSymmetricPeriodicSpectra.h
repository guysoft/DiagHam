////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//       class for periodic average spectra with XY reflexion symmetry        //
//                                                                            //
//                        last modification : 04/04/2004                      //
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



#ifndef XYREFLEXIONSYMMETRICPERIODICSPECTRA_H
#define XYREFLEXIONSYMMETRICPERIODICSPECTRA_H

#include "config.h"

#include "HilbertSpace/QuantumDotHilbertSpace/XYReflexionSymmetricPeriodic3DOneParticle.h"


class XYReflexionSymmetricPeriodicSpectra
{
 protected:
  // wave function basis dimension in the x direction
  int NbrStateX;
  int LowerImpulsionX;
  
  // wave function basis dimension in the y direction
  int NbrStateY;
  int LowerImpulsionY;
  
  // wave function basis dimension in the z direction
  int NbrStateZ;
  int LowerImpulsionZ;

  double*** RealCoefficients;
  double*** ImaginaryCoefficients;

 public:

  // constructor from a Hilbert space and a file
  //
  // space = Hilbert space describing the particle
  // fileName = name of the state file
  XYReflexionSymmetricPeriodicSpectra(XYReflexionSymmetricPeriodic3DOneParticle* space, char* fileName);

  // get the value of impulsion operators with another wavefunction <this|p|another>
  //
  // space = Hilbert space describing the other particle
  // fileName = the file to stock the other function
  // sizeX, sizeY, sizeZ = size of sample in X, Y and Z directions
  // impulsionX, impulsionY, impulsionZ = reference to the return values
  void GetImpulsion(XYReflexionSymmetricPeriodic3DOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &realImpulsionX, double &imaginaryImpulsionX, double &realImpulsionY, double &imaginaryImpulsionY, double &realImpulsionZ, double &imaginaryImpulsionZ);
  
};


#endif
