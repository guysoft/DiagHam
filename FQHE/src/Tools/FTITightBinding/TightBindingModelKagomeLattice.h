////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of tight binding model for the Kagome lattice          //
//                                                                            //
//                        last modification : 07/05/2012                      //
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


#ifndef TIGHTBINDINGMODELKAGOMELATTICE_H
#define TIGHTBINDINGMODELKAGOMELATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelKagomeLattice : public Abstract2DTightBindingModel
{

 protected:

  // hopping amplitude between neareast neighbor sites
  double NNHopping;
  // hopping amplitude between next neareast neighbor sites
  double NextNNHopping;
  // spin orbit coupling to neareast neighbor sites
  double NNSpinOrbit;
  // spin orbit coupling to next neareast neighbor sites
  double NextNNSpinOrbit;
  
  // four times the sublattice staggered chemical potential 
  double MuS;

 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t1 = real part of the hopping amplitude between neareast neighbor sites
  // t2 = real part of the hopping amplitude between next neareast neighbor sites
  // lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
  // lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
  // mus = sublattice chemical potential on A1 sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelKagomeLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double lambda1, double lambda2, double mus, 
				 double gammaX, double gammaY, bool storeOneBodyMatrices = true);

  // destructor
  //
  ~TightBindingModelKagomeLattice();

 protected :

  // compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void ComputeBandStructure(long minStateIndex = 0l, long nbrStates = 0l);

};


#endif
