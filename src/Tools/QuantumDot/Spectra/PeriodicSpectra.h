////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                      class for periodic average spectra                    //
//                                                                            //
//                        last modification : 12/20/2003                      //
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



#ifndef PERIODICSPECTRA_H
#define PERIODICSPECTRA_H

#include "config.h"

#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"


class PeriodicSpectra
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

  // wave function overlaps on a cell in the x direction 
  double* RealWaveFunctionOverlapX;
  double* ImaginaryWaveFunctionOverlapX;
  double* RealSquareOverlapX;
  double* ImaginarySquareOverlapX;

  // wave function overlaps on a cell in the y direction 
  double* RealWaveFunctionOverlapY;
  double* ImaginaryWaveFunctionOverlapY;
  double* RealSquareOverlapY;
  double* ImaginarySquareOverlapY;

  // wave function overlaps on a cell in the z direction 
  double* RealWaveFunctionOverlapZ;
  double* ImaginaryWaveFunctionOverlapZ;
  double* RealSquareOverlapZ;
  double* ImaginarySquareOverlapZ;
  
  double*** RealCoefficients;
  double*** ImaginaryCoefficients;

 public:

  // constructor from a Hilbert space and a file
  //
  // space: Hilbert space describing the particle
  // fileName: name of the state file
  PeriodicSpectra(Periodic3DOneParticle* space, char* fileName);

  // get mean value in X direction
  //
  // squareX: reference to the mean square value in X direction
  // return: position in 1.0 scale
  double GetMeanValueX(double& squareX);

  // get mean value in Y direction
  //
  // squareY: reference to the mean square value in Y direction
  // return: position in 1.0 scale
  double GetMeanValueY(double& squareY);

  // get mean value in Z direction
  //
  // squareZ: reference to the mean square value in Z direction
  // return: position in 1.0 scale 
  double GetMeanValueZ(double& squareZ);

  // get the wave function value of a state at a given point
  //
  // x, y, z : the position of the point
  // SizeX, SizeY, SizeZ : the 3D-sizes of the sample
  // Real, Imaginary : references to the real and imaginary components of the wave function
  void WaveFunctionValue(double x, double SizeX, double y, double SizeY, double z, double SizeZ, double& Real, double& Imaginary);
};


#endif
