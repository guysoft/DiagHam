////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 1D tight binding model                  //
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


#ifndef ABSTRACT1DTIGHTBINDINGMODEL_H
#define ABSTRACT1DTIGHTBINDINGMODEL_H


#include "config.h"
#include "Tools/FTITightBinding/AbstractTightBindingModel.h"


class Abstract1DTightBindingModel : public AbstractTightBindingModel
{

 protected:

   // number of sites in the x direction
  int NbrSiteX;

  // numerical factor for momentum along x
  double KxFactor;

  // boundary condition twisting angle along x
  double GammaX;

  // embedding of sublattices relative to the unit cell reference point along x
  RealVector EmbeddingX;

  // One body eigenstate basis associated to each point of the band structure, the array index corresponds to the linearized momentum
  ComplexMatrix* OneBodyBasis;

  // energy spectrum of the band structure, first index is the band index, second index is the linearized momentum
  double** EnergyBandStructure;

 public:

  // default constructor
  //
  Abstract1DTightBindingModel();

  // destructor
  //
  ~Abstract1DTightBindingModel();

  // get the energy at a given momentum of the band structure
  //
  // bandIndex = index of the band to consider
  // momentumIndex = linearized momentum
  // return value = energy
  virtual double GetEnergy(int bandIndex, int momentumIndex);

  // ask if the one body transformation matrices are available
  //
  // return value = true if the one body transformation matrices are available
  virtual bool HaveOneBodyBasis();

  // get  the one body transformation matrix corresponding to a given momentum of the band structure
  //
  // momentumIndex = linearized momentum
  // return value = reference on the one body transformation matrix
  virtual ComplexMatrix& GetOneBodyMatrix(int momentumIndex);
  
  // get the number of sites in the x direction
  //
  // return value = number of sites in the x direction
  int GetNbrSiteX();

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

// get the energy at a given momentum of the band structure
//
// bandIndex = index of the band to consider
// momentumIndex = linearized momentum
// return value = energy

inline double Abstract1DTightBindingModel::GetEnergy(int bandIndex, int momentumIndex)
{
  return this->EnergyBandStructure[bandIndex][momentumIndex];
}

// get  the one body transformation matrix corresponding to a given momentum of the band structure
//
// momentumIndex = linearized momentum
// return value = reference on the one body transformation matrix

inline ComplexMatrix& Abstract1DTightBindingModel::GetOneBodyMatrix(int momentumIndex)
{
  return this->OneBodyBasis[momentumIndex];
}

// get the number of sites in the x direction
//
// return value = number of sites in the x direction

inline int Abstract1DTightBindingModel::GetNbrSiteX()
{
  return this->NbrSiteX;
}

// ask if the one body transformation matrices are available
//
// return value = true if the one body transformation matrices are available

inline bool Abstract1DTightBindingModel::HaveOneBodyBasis()
{
  return (this->OneBodyBasis != 0);
}

#endif
