////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of tight binding model for the 2D atomic limit         //
//                                                                            //
//                        last modification : 03/07/2014                      //
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


#ifndef TIGHTBINDINGMODEL2DATOMICLIMITLATTICE_H
#define TIGHTBINDINGMODEL2DATOMICLIMITLATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModel2DAtomicLimitLattice : public Abstract2DTightBindingModel
{

 protected:

  // nbrSitesUnitCell = number of sites per unit cell
  int NbrSitesUnitCell;

  // on-site chemical potentials
  double* ChemicalPotentials;


 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSitesUnitCell = number of sites per unit cell
  // chemicalPotentials = on site chemical potentials
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModel2DAtomicLimitLattice(int nbrSiteX, int nbrSiteY,
					int nbrSitesUnitCell, double* chemicalPotentials,
					double gammaX, double gammaY, 
					AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

  // constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // nbrSitesUnitCell = number of sites per unit cell
  // NbrLowChemicalPotentials = number of sites with chemical potential -1, the other sites having chemical potential +1
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModel2DAtomicLimitLattice(int nbrSiteX, int nbrSiteY,
					int nbrSitesUnitCell, int nbrLowChemicalPotentials,
					double gammaX, double gammaY, 
					AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

  // destructor
  //
  ~TightBindingModel2DAtomicLimitLattice();

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif
