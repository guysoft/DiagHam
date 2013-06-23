////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                          class author Cecile Repellin                      //
//                                                                            //
//               class of tight binding model for the Kagome lattice          //
//                      with tilted boundary conditions                       //
//                        last modification : 14/05/2013                      //
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


#ifndef TIGHTBINDINGMODELKAGOMELATTICETILTED_H
#define TIGHTBINDINGMODELKAGOMELATTICETILTED_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelKagomeLatticeTilted : public Abstract2DTightBindingModel
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
  
  //first coordinate of the first spanning vector of the tilted lattice
  int Nx1;
  //second coordinate of the first spanning vector of the tilted lattice
  int Ny1;
  //first coordinate of the second spanning vector of the tilted lattice
  int Nx2;
  //second coordinate of the second spanning vector of the tilted lattice
  int Ny2;
  //array of projected momenta
  double** ProjectedMomenta;
  //second coordinate in momentum space of the second spanning vector of the reciprocal lattice
  int Offset;
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
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelKagomeLatticeTilted(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double lambda1, double lambda2, double mus, 
				 double gammaX, double gammaY, 
				 AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

  // destructor
  //
  ~TightBindingModelKagomeLatticeTilted();

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);
  
  //computes all the values of the momentum projected and stores them in a double array
  //
  virtual void ComputeAllProjectedMomenta();
  
  //Computes value of projected momentum along the lattice directions
  //
  //kx = first coordinate of the given point in the Brillouin zone
  //ky = second coordinate of the given point in the Brillouin zone
  //latticeComponent = index of the lattice vector along which the projection is to be performed
  //return value = projected momentum
  virtual double GetProjectedMomentum(int kx, int ky, int latticeComponent);

};

  inline double TightBindingModelKagomeLatticeTilted::GetProjectedMomentum(int kx, int ky, int latticeComponent)
  {
    return this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky)][latticeComponent];
  }
#endif