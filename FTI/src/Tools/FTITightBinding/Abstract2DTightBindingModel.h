////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 2D tight binding model                  //
//                                                                            //
//                        last modification : 01/05/2012                      //
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


#ifndef ABSTRACT2DTIGHTBINDINGMODEL_H
#define ABSTRACT2DTIGHTBINDINGMODEL_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract1DTightBindingModel.h"


class Abstract2DTightBindingModel : public Abstract1DTightBindingModel
{

 protected:

   // number of sites in the y direction
  int NbrSiteY;

  // numerical factor for momentum along y
  double KyFactor;

  // boundary condition twisting angle along y
  double GammaY;

 public:

  // default constructor
  //
  Abstract2DTightBindingModel();

  // destructor
  //
  ~Abstract2DTightBindingModel();

  // get the linearized momentum index
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // return value = inearized momentum index
  virtual int GetLinearizedMomentumIndex(int kx, int ky);

  // get momentum value from a linearized momentum index
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // return value = inearized momentum index
  virtual void GetLinearizedMomentumIndex(int index, int& kx, int& ky);

  // get the linearized momentum index, without assuming k to be in the first BZ
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // return value = inearized momentum index
  virtual int GetLinearizedMomentumIndexSafe(int kx, int ky);

  // get momentum value from a linearized momentum index, without assuming k to be in the first BZ
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // return value = inearized momentum index
  virtual void GetLinearizedMomentumIndexSafe(int index, int& kx, int& ky);

  // get the number of sites in the y direction
  //
  // return value = number of sites in the y direction
  int GetNbrSiteY();

  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

  // write the full band structure information in an ASCII file
  //
  // fileName = name of the output file 
  // return value = true if no error occured  
  virtual bool WriteBandStructureASCII(char* fileName);

  // compute the Chern number of a given band
  //
  // band = band index
  // return value = Chern number
  virtual double ComputeChernNumber(int band);

};

// get the linearized momentum index
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// return value = linearized momentum index

inline int Abstract2DTightBindingModel::GetLinearizedMomentumIndex(int kx, int ky)
{
  return ((kx * this->NbrSiteY) + ky);
}

// get momentum value from a linearized momentum index
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// return value = inearized momentum index

inline void Abstract2DTightBindingModel::GetLinearizedMomentumIndex(int index, int& kx, int& ky)
{
  kx = index / this->NbrSiteY;
  ky = index % this->NbrSiteY;
}

// get the linearized momentum index, without assuming k to be in the first BZ
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// return value = linearized momentum index

inline int Abstract2DTightBindingModel::GetLinearizedMomentumIndexSafe(int kx, int ky)
{
  while (kx < 0)
      kx += this->NbrSiteX;
  kx %= this->NbrSiteX;
  while (ky < 0)
      ky += this->NbrSiteY;
  ky %= this->NbrSiteY;
  return ((kx * this->NbrSiteY) + ky);
}

// get momentum value from a linearized momentum index, without assuming k to be in the first BZ
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// return value = inearized momentum index

inline void Abstract2DTightBindingModel::GetLinearizedMomentumIndexSafe(int index, int& kx, int& ky)
{
  int n = this->NbrSiteX * this->NbrSiteY;
  while (index < 0)
      index += n;
  index %= n;
  kx = index / this->NbrSiteY;
  ky = index % this->NbrSiteY;
}

// get the number of sites in the y direction
//
// return value = number of sites in the y direction

inline int Abstract2DTightBindingModel::GetNbrSiteY()
{
  return this->NbrSiteY;
}

#endif
