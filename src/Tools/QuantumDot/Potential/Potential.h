////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                           base class for potential                         //
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


#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Color/PicRGB.h"
#include "BitmapPicture/AbstractBitmapPicture.h"

#include "config.h"

#include <iostream>
#include <fstream>

using std::ostream;
using std::ifstream;

class Potential
{
 protected:

  // number of points in x direction
  int NumberX;

  // number of points in y direction
  int NumberY;

  // number of points in z direction
  int NumberZ;

  // type of atoms, true: alloy atom, false: host atom
  // order: zyx
  bool*** Alloy;


 public:

  // concentration
  double Concentration;

  // default constructor
  //
  Potential();

  // constructor from a 3D dimension
  // 
  // x, y, z: three dimensions in respective direction
  Potential(int x, int y, int z);

  // virtual destructor
  //
  virtual ~Potential();
  
  // print the potential dimension
  //
  // output = an ostream to display
  // display Nx, Ny and Nz in an output
  virtual ostream& PrintDimension(ostream& output);

  // print the diagram of atomic distribution in an output
  //
  // output = an ostream to display
  // display Alloy arry as z blocks, each block contains x columns and y lines
  virtual ostream& PrintDiagram(ostream& output);
  
  // read the atom diagram from an input, which has z blocks each block contains x columns and y lines
  //
  // fileName = name of the storage file
  virtual void ReadDiagram(char* fileName);

  // generate an arbitrary distribution of atoms with a given proportion
  //
  // proportion = nominal proportion
  virtual void ArbitraryDistribution(double proportion);

  // calculate the mean first neighbors in a potential diagram
  //
  // Total: reference to the total number of InN 
  // return = the mean number of first neighbors (6)
  double MeanFirstNeighbors(double& Total);
  
  // save the diagram to a bitmap file
  //
  // u: the number of monolayers omitted under the considered structure
  // a: the number of monolayers omitted above the considered structure
  // startX: point to start in diagram in X direction
  // endX: point to end in diagram in X direction
  // startY: point to start in diagram in Y direction
  // endY: point to end in diagram in Y direction
  // choice: choice of orientation: 1 for XY, 2 for XZ, 3 for YZ
  // sizeX: the size in X direction of each cell in pixel
  // sizeY: the size in Y direction of each cell in pixel  
  // InN:   color (in RGB definition) of InN cell
  // GaN:   color (in RGB definition) of GaN cell
  // background : color (in RGB definition) of background
  // NbrX: number of cell displayed in X direction
  // fileName: name of the file to store the picture

  bool Potential::SaveBmpPicture(int under, int above, int startX, int endX, int startY, int endY, int choice, int sizeX, int sizeY, PicRGB& InN, PicRGB& GaN, PicRGB& background, int NbrX, char* fileName);
};

// fill a cell with a given color
//
// startX: X coordination to start
// sizeX:  size in X direction
// startY: Y coordination to start
// sizeY:  size in Y direction
// Col: color in RGB definition
// picture: pointer to picture which is filled
void CellFill(int startX, int sizeX, int startY, int sizeY, PicRGB& Col, AbstractBitmapPicture* picture);

#endif
