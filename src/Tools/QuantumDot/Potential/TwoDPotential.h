////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                            class for 2D potential                          //
//                                                                            //
//                        last modification : 09/15/2003                      //
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


#ifndef TWODPOTENTIAL_H
#define TWODPOTENTIAL_H

#include "Tools/QuantumDot/Potential/Potential.h"

#include "config.h"

using std::ostream;

class TwoDPotential : public Potential
{
 protected:

  // the potential in 2D, the first index is for y and the second for x
  double** Potential;

 public:

  // default constructor
  //
  TwoDPotential();

  // construct a frame potential
  //
  // x, y, z: three dimensions in respective direction
  TwoDPotential(int x, int y, int z);

  // destructor
  //
  ~TwoDPotential();

  // get reference of a given potential
  //
  // x, y = x and y coordinate
  // return : potential value of the cell
  double& operator() (int x, int y);

  // read 2D potential from a file
  //
  // FileName = name of the data file
  // convention: x columns and y lines 
  void ReadTwoDPotential(char* FileName);
  
  // read 2D potential from a file when the number of lines and columns are known
  //
  // FileName = nam of the data file
  // convention: x columns and y lines
  void ReadTwoDPotential(char* FileName, int NX, int NY, int NZ);
  
  // print the potential in a output
  //
  // output = output (file, screen, ...)
  // return: the ostream, which contains the Potential array
  // convention: there are x columns and y lines
  ostream& PrintPotential(ostream& output);

  // construct a well potential with uniform probability
  //
  // proportion = proportion of the alloy
  // offset = offset of the alloy and the host meterials, withour virtual crystal approximation
  // weight = weight of potential in each mono-layer
  // scratch = true if constructed from scratch, false if from existing diagram
  void UniformWell(double proportion, double offset, double* weight, bool scratch);

  // construct a dot potential with uniform probability
  void UniformPyramidDot(double proportion);
};

// get reference of a given potential
//
// x = x coordinate, y = y coordinate, all are 0-based
// return : potential value of the cell

inline double& TwoDPotential::operator() (int x, int y)
{
  return (this->Potential)[y][x];
}

#endif
