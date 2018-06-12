////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                         Class author : Cecile Repellin                     //
//                                                                            //
//           class of tight binding model for trilayer grahene-hBN            //
//                with a relative twist angle and Moire coupling              //
//                      with valley and spin conservation                     //
//               (must be TR conjugated to obtain other valley index)         //
//                   last modification : 04/06/2018                           //
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


#ifndef TIGHTBINDINGMODELTRILAYERGRAPHENETWISTBN_H
#define TIGHTBINDINGMODELTRILAYERGRAPHENETWISTBN_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModelBilayerTwistBilayer.h"


class TightBindingModelTrilayerGrapheneTwistBN : public TightBindingModelBilayerTwistBilayer
{

 protected:
     

 public:

// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t = tunneling amplitude within one layer
// t3 = trigonal warping
// g1 = tunneling amplitude within one layer
// tM = tunneling amplitude between layers
// uVoltage = amplitude of voltage
// theta = twisting angle     
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelTrilayerGrapheneTwistBN(int nbrSiteX, int nbrSiteY, int nbrPointsX, int nbrPointsY, double uVoltage, double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices);

  // destructor
  //
  ~TightBindingModelTrilayerGrapheneTwistBN();

  // compute the Bloch hamiltonian at a point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // ky = momentum along the x axis
  // return value = Bloch hamiltonian
  virtual HermitianMatrix ComputeBlochHamiltonian(double kx, double ky);
  
  
  // compute the form factor for the density operator 
  // 
  // kx = momentum along x of annihilation operator
  // ky = momentum along y of creation operator
  // qx = momentum transfer along x direction
  // qy = momentum transfer along y direction
  // valleyIndex = valley index of density operator
  virtual Complex ComputeDensityFormFactor(int kx, int ky, int qx, int qy, int valleyIndex);
  

 protected :

  // core part that computes the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif
