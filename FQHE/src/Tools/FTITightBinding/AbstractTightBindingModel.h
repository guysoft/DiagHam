////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of abstract tight binding model                    //
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


#ifndef ABSTRACTTIGHTBINDINGMODEL_H
#define ABSTRACTTIGHTBINDINGMODEL_H


#include "config.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;


class AbstractTightBindingModel
{

  friend class FTITightBindingModelBandStructureOperation;

 protected:

  // number of bands (including any internal degree of freedom such as spin)
  int NbrBands;

  // number of states per band
  long NbrStatePerBand;

 public:

  // destructor
  //
  virtual ~AbstractTightBindingModel();

  // get the number of states per band
  //
  // return value =  number of states per band
  virtual long GetNbrStatePerBand();

  // get the energy at a given momentum of the band structure
  //
  // bandIndex = index of the band to consider
  // momentumIndex = linearized momentum
  // return value = energy
  virtual double GetEnergy(int bandIndex, int momentumIndex) = 0;

  // ask if the one body transformation matrices are available
  //
  // return value = true if the one body transformation matrices are available
  virtual bool HaveOneBodyBasis() = 0;

  // get  the one body transformation matrix corresponding to a given momentum of the band structure
  //
  // momentumIndex = linearized momentum
  // return value = reference on the one body transformation matrix
  virtual ComplexMatrix& GetOneBodyMatrix(int momentumIndex) = 0;

  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

  // write the full band structure information i a binary file
  //
  // fileName = name of the output file 
  // return value = true if no error occured  
  virtual bool WriteBandStructure(char* fileName);

  // compute the spread of a given band
  // 
  // band = band index
  // return value = band spread 
  virtual double ComputeBandSpread(int band);

  // compute the direct gap between two bands
  // 
  // band1 = index of the lower band 
  // band2 = index of the upper band (-1 if it has to be set to band1 + 1)
  // return value =  direct band gap
  virtual double ComputeDirectBandGap(int band1, int band2 = -1);

 protected:

  // write an ASCII header that describes the tight binding model
  // 
  // output = reference on the output stream
  // commentChar = optional ASCII character to add in front of each header line
  // return value  = reference on the output stream
  virtual ostream& WriteASCIIHeader(ostream& output, char commentChar = '\0');

  // compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void ComputeBandStructure(long minStateIndex = 0l, long nbrStates = 0l) = 0;

};

// get the number of states per band
//
// return value =  number of states per band

inline long AbstractTightBindingModel::GetNbrStatePerBand()
{
  return this->NbrStatePerBand;
}

#endif
