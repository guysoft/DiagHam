////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 3D tight binding model                  //
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


#ifndef ABSTRACT3DTIGHTBINDINGMODEL_H
#define ABSTRACT3DTIGHTBINDINGMODEL_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


using std::cout;
using std::endl;


class Abstract3DTightBindingModel : public Abstract2DTightBindingModel
{

 protected:

   // number of sites in the z direction
  int NbrSiteZ;
  // number of sites in the direction perpendicular to X
  int NbrSiteYZ;

  // numerical factor for momentum along z
  double KzFactor;

  // boundary condition twisting angle along z
  double GammaZ;

  // embedding of sublattices relative to the unit cell reference point along z
  RealVector EmbeddingZ;

 public:

  // default constructor
  //
  Abstract3DTightBindingModel();

  // destructor
  //
  ~Abstract3DTightBindingModel();

  // get the linearized momentum index
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // ky = momentum along the z direction
  // return value = inearized momentum index
  virtual int GetLinearizedMomentumIndex(int kx, int ky, int kz);

  // get momentum value from a linearized momentum index
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // kz = reference on the momentum along the z direction
  // return value = inearized momentum index
  void GetLinearizedMomentumIndex(int index, int& kx, int& ky, int& kz);

  // get the linearized momentum index, without assuming k to be in the first BZ
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // ky = momentum along the z direction
  // return value = inearized momentum index
  virtual int GetLinearizedMomentumIndexSafe(int kx, int ky, int kz);

  // get momentum value from a linearized momentum index, without assuming k to be in the first BZ
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // kz = reference on the momentum along the z direction
  // return value = inearized momentum index
  void GetLinearizedMomentumIndexSafe(int index, int& kx, int& ky, int& kz);

  // get the number of sites in the z direction
  //
  // return value = number of sites in the z direction
  int GetNbrSiteZ();

  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

 protected:

  // write an header that describes the tight binding model
  // 
  // output = reference on the output stream
  // return value  = reference on the output stream
  virtual ofstream& WriteHeader(ofstream& output);

};

// get the linearized momentum index
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// kz = momentum along the z direction
// return value = linearized momentum index

inline int Abstract3DTightBindingModel::GetLinearizedMomentumIndex(int kx, int ky, int kz)
{
  return (((kx * this->NbrSiteY) + ky) * this->NbrSiteZ + kz);
}

// get momentum value from a linearized momentum index
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// kz = reference on the momentum along the z direction
// return value = inearized momentum index

inline void Abstract3DTightBindingModel::GetLinearizedMomentumIndex(int index, int& kx, int& ky, int& kz)
{
  kx = index / this->NbrSiteYZ;
  ky = index % this->NbrSiteYZ;
  kz = ky % this->NbrSiteZ;
  ky /= this->NbrSiteZ;
}

// get the linearized momentum index, without assuming k to be in the first BZ
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// kz = momentum along the z direction
// return value = linearized momentum index

inline int Abstract3DTightBindingModel::GetLinearizedMomentumIndexSafe(int kx, int ky, int kz)
{
  while (kx < 0)
      kx += this->NbrSiteX;
  kx %= this->NbrSiteX;
  while (ky < 0)
      ky += this->NbrSiteY;
  ky %= this->NbrSiteY;
  while (kz < 0)
      kz += this->NbrSiteZ;
  kz %= this->NbrSiteZ;
  return (((kx * this->NbrSiteY) + ky) * this->NbrSiteZ + kz);
}

// get momentum value from a linearized momentum index, without assuming k to be in the first BZ
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// kz = reference on the momentum along the z direction
// return value = inearized momentum index

inline void Abstract3DTightBindingModel::GetLinearizedMomentumIndexSafe(int index, int& kx, int& ky, int& kz)
{
  int n = this->NbrSiteX * this->NbrSiteYZ;
  while (index < 0)
      index += n;
  index %= n;
  kx = index / this->NbrSiteYZ;
  ky = index % this->NbrSiteYZ;
  kz = ky % this->NbrSiteZ;
  ky /= this->NbrSiteZ;
}

// get the number of sites in the z direction
//
// return value = number of sites in the z direction

inline int Abstract3DTightBindingModel::GetNbrSiteZ()
{
  return this->NbrSiteZ;
}

#endif
