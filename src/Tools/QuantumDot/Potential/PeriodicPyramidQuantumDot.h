////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                    class for periodic pyramid quantum dot                  //
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


#ifndef PERIODICPYRAMIDQUANTUMDOTPOTENTIAL_H
#define PERIODICPYRAMIDQUANTUMDOTPOTENTIAL_H

#include "Tools/QuantumDot/Potential/ThreeDPotential.h"

#include "config.h"

class PeriodicPyramidQuantumDotPotential : public ThreeDPotential
{

 protected:

  // width of wetting layer
  int WettingWidth;

  // base radius
  int BaseRadius;

  // top radius
  int TopRadius;

 

 public:

  // constructor
  PeriodicPyramidQuantumDotPotential(int NbrCellX, int NbrCellY, int NbrCellZ, double Lz, int u, int a, int rb, int rt, int w, double offset, double concentration, double piezofield, bool scratch, char* logfile);

  // destructor
  ~PeriodicPyramidQuantumDotPotential();

  bool GeneratePotential(double C, double F, double Lz, double offset, bool scratch, char* filename);

  // determine if a cell is in the dot (wetting layer is included)
  // x: x coordinate of the cell
  // y: y coordinate of the cell
  // z: z coordinate of the cell
  //
  // return : true if the cell is in the dot, false otherwise
  bool InTheDot(int x, int y, int z);
};


#endif
