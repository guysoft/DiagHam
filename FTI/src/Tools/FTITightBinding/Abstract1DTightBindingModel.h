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

  // get the linearized momentum index
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // return value = linearized momentum index
  virtual int GetLinearizedMomentumIndex(int kx, int ky);

  // get momentum value from a linearized momentum index
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // return value = linearized momentum index
  virtual void GetLinearizedMomentumIndex(int index, int& kx, int& ky);

  // get the energy at a given momentum of the band structure
  //
  // bandIndex = index of the band to consider
  // momentumIndex = linearized momentum
  // return value = energy
  virtual double GetEnergy(int bandIndex, int momentumIndex);

  // get all the energies of a band, sorted from the smallest to the largest
  //
  // energies = reference to the array where the energies will be stored (the allocation is done by the method)
  // momenta = reference to the array where the linearized momenta associated to each energy will be stored (the allocation is done by the method)
  // bandIndex = index of the band to consider
  virtual void GetEnergies(double*& energies, int*& momenta, int bandIndex);

  // get all the energies, sorted from the smallest to the largest
  //
  // energies = reference to the array where the energies will be stored (the allocation is done by the method)
  // momenta = reference to the array where the linearized momentum associated to each energy will be stored (the allocation is done by the method)
  // bandIndices = reference to the array where the band index associated to each energy will be stored (the allocation is done by the method)
  virtual void GetAllEnergies(double*& energies, int*& momenta, int*& bandIndices);

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

  // flatten the energy spectrum
  //
  // flatteningFactor = flattening factor applied to each band (the band is rescale with respect to its average value)
  // bandShift = shift each band by bandShift times the band index
  virtual void FlattenBands(double flatteningFactor, double bandShift);

  // evaluate the two point correlation function in a given region
  //
  // maxX = x coordinate of the region upper right corner 
  // maxY = y coordinate of the region upper right corner 
  // occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
  // nbrOccupiedMomenta = number of occupied momenta
  // bandIndex = index of the band to consider
  // return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)
  virtual HermitianMatrix EvaluateFullTwoPointCorrelationFunction(int maxX, int maxY, int* occupiedMomenta, int nbrOccupiedMomenta, int bandIndex);

  // evaluate the mixed two point correlation function in a given region, assuming translation invariance along one direction
  //
  // maxX = length along the borken translation direction of the region 
  // ky = momentum along the translation invariant direction
  // occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
  // bandIndices = array that gives the band index of each occupied state
  // nbrOccupiedMomenta = number of occupied momenta
  // return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)
  virtual HermitianMatrix EvaluateFullMixedTwoPointCorrelationFunctionWithK(int maxX, int ky, int* occupiedMomenta, int* bandIndices, int nbrOccupiedMomenta);

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

inline ComplexMatrix&  Abstract1DTightBindingModel::GetOneBodyMatrix(int momentumIndex)
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

// get the linearized momentum index
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// return value = linearized momentum index

inline int Abstract1DTightBindingModel::GetLinearizedMomentumIndex(int kx, int ky)
{
  return kx;
}

// get momentum value from a linearized momentum index
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// return value = linearized momentum index

inline void Abstract1DTightBindingModel::GetLinearizedMomentumIndex(int index, int& kx, int& ky)
{
  kx = index;
  ky = 0;
}

#endif
