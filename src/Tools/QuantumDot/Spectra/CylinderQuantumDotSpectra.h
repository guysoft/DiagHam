////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//        class for periodic average spectra with Fourier-Bessel basis        //
//                                                                            //
//                        last modification : 07/05/2004                      //
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


#ifndef CYLINDERQUANTUMDOTSPECTRA_H
#define CYLINDERQUANTUMDOTSPECTRA_H


#include "config.h"

#include "HilbertSpace/QuantumDotHilbertSpace/VerticalPeriodicParticleInMagneticField.h"

class CylinderQuantumDotSpectra
{
 protected:

  // quantum number of kinetic momentum in Z direction
  int NumberM;
  
  // wave function basis dimension in the plane
  int NbrStateR;
  
  // wave function basis dimension in the z direction
  int NbrStateZ;
  int LowerImpulsionZ;
  
  // magnetic field in Z direction
  double Bz;

  double** RealCoefficients;
  double** ImaginaryCoefficients;

 public:

  // constructor from a Hilbert space and a file
  //
  // space = Hilbert space describing the particle
  // fileName = name of the state file
  // bz = magnetic field in Z direction
  CylinderQuantumDotSpectra(VerticalPeriodicParticleInMagneticField* space, char* fileName, double bz);

  // get the value of impulsion operators with another wavefunction <this|p|another>
  //
  // space = Hilbert space describing the other particle
  // fileName = the file to stock the other function
  // sizeZ = size of sample in Z direction
  // sizeR = size of the super-cylinder in plane
  // impulsionX, impulsionY, impulsionZ = reference to the return values
  void GetImpulsion(VerticalPeriodicParticleInMagneticField* space, char* fileName, double sizeZ, double sizeR, double &realImpulsionX, double &imaginaryImpulsionX, double &realImpulsionY, double &imaginaryImpulsionY, double &realImpulsionZ, double &imaginaryImpulsionZ);
  
};

#endif
