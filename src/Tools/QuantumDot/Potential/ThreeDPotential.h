////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                            class for 3D potential                          //
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


#ifndef THREEDPOTENTIAL_H
#define THREEDPOTENTIAL_H

#include "Tools/QuantumDot/Potential/Potential.h"

#include "config.h"

#include <iostream>
#include <fstream>

using std::ifstream;

class ThreeDPotential : public Potential
{ 
  friend class QuantumDots3DHamiltonian;
  friend class PeriodicQuantumDots3DHamiltonian;
  friend class NewPeriodicQuantumDots3DHamiltonian;

 protected:

  // the potential in 3D, the first index is for z, the second for y and the third for x
  double*** Potential;

  // number of mono-layers under the considered structure
  int under;

  // size 
  double* UnderSize;

  // potential along growth direction, under the considered structure
  double* UnderPotential;

  // number of mono-layers above the considered structure
  int above;

  // size
  double* AboveSize;
  
  // potential along the growth direction, above the considered structure
  double* AbovePotential;

 public:

  // default constructor
  //
  ThreeDPotential();

  // construct a frame potential
  //
  // x, y, z: three dimensions in respective direction
  ThreeDPotential(int x, int y, int z);

  // construct a frame potential with more geometric parameter
  //
  // u and a: number of mono-layers under and above the structure respectively
  ThreeDPotential(int x, int y, int z, int u, int a);

  // destructor
  //
  ~ThreeDPotential();

  // get reference of a given potential
  //
  // x, y, z : respective coordinate
  // return : potential value of the cell
  double& operator() (int x, int y, int z);

  // construct a well with arbitrary distribution of alloy in the well
  //
  // down_field = electric field under the considered structure
  // field = electric field in the considered structure
  // up_field = electric field above the considered structure
  // cell = the height of a monolayer
  // proportion = proportion of alloy in the considered structure
  // offset = offset between host material and alloy
  // scratch = true if construct from scratch
  // fileName = name of the file to store parameters
  // 0 reference: the monolayer just before the downside of the well
  void UniformWell(double down_field, double field, double upfield, double cell, double proportion, double offset, bool scratch, char* fileName);

  // construct a well with a pocket alloy upon the well
  //
  // height = height of the pocket
  // width = width of the pocket
  // down_field = electric field under the considered structure
  // field = electric field in the considered structure
  // up_field = electric field above the considered structure
  // cell = the height of a monolayer
  // offset = offset between host material and alloy
  // scratch = true if construct from scratch
  // fileName = name of the file to store parameters
  // 0 reference: the monolayer just before the downside of the well
  // ATTENTION TO Resize of Potential
  void FluctuatedWell(int height, int width, double down_field, double field, double up_field, double cell, double offset, char* fileName);

  // potential 0 reference: upside monolayer
  // ATTENTION TO Resize of Potential, UnderSize, UnderPotential, AboveSize & AbovePotential
  void IEFFluctuatedWell(double p, double down_field, double field, double up_field, double cell, double offset, char* fileName);

  // construct a well with segragation distribution of alloy
  // down_field = electric field under the considered structure
  // field = electric field in the considered structure
  // up_field = electric field above the considered structure
  // cell = the height of a monolayer
  // proportion1 = probability without InN neighborhood
  // proportion2 = probability with InN neighborhood
  // offset = offset between host material and alloy
  // scratch = true if construct from scratch
  void SegregationWell(double down_field, double field, double up_field, double cell, double proportion1, double proportion2, double offset, bool scratch, char* fileName);
  
  // m, n, h : three coordinations of the considered cell
  // return true if there is any InN neighborhood
  bool Neighborhood(int m, int n, int h);

  // m, n, h : three coordinations of the considered cell
  // return true if there is any InN neighborhood
  bool NeighborhoodBis(int m, int n, int h);

  // construct a truncated hexagonal pyramid dot used for nitride materials simulation with uniform distribution
  //
  // proportion = proportion of alloy materials
  // Rb = base radius, Rt = top radius
  // w = thickness of wetting layer
  // down_field, wetting_field, dot_field, up_field = the electric field in respective part
  // offset = offset bulk-alloy material
  // c: height of a monolayer
  // scratch = true if constructed from nothing, else from existing Diagram
  // fileName = file to store parameters
  void ArbitraryPyramidDot(double proportion, int Rb, int Rt, int w, double down_field, double wetting_field, double dot_field, double up_field, double offset, double c, bool scratch, char* fileName);

  // construct a truncated hexagonal pyramid dot used for nitride materials simulation with uniform distribution
  //
  // proportion1 = probability without InN neighborhood
  // proportion2 = probability with InN neighborhood 
  // Rb = base radius, Rt = top radius
  // w = thickness of wetting layer
  // down_field, wetting_field, dot_field, up_field = the electric field in respective part
  // offset = offset bulk-alloy material
  // c: height of a monolayer
  // scratch = true if constructed from nothing, else from existing Diagram
  // fileName = file to store parameters
  void SegregationPyramidDot(double proportion1, double proportion2, int Rb, int Rt, int w, double down_field, double wetting_field, double dot_field, double up_field, double offset, double c, bool scratch, char* fileName);

  // print the 3-D potential in an output
  //
  // output = output (file, screen, ...)
  // output has z block, each block has y lines and x columns
  ostream& PrintPotential(ostream& output);

  // print the 3-D potential with potential array along growth direction in an output
  //  
  // output = output (file, screen, ...)
  // output has z block, each block has y lines and x columns
  ostream& PrintPotentialWithField(ostream& output);

  // construct a potential from a file
  //
  // fileName = name of the data file
  // convention: z blocks, each block has y lines and x columns
  // order of potential: zyx
  void ReadPotential(char* fileName);

  // construct a potential from a file with electric field
  //
  // fileName = name of the storage file
  // convention: z blocks, each block has y lines and x columns
  // order of potential: zyx
  void ReadPotentialWithField(char* fileName);

};

// get reference of a given potential
//
// i : y coordinate, j : x coordinate
// return : potential value of the cell

inline double& ThreeDPotential::operator() (int x, int y, int z)
{
  return (this->Potential)[z][y][x];
}

#endif
