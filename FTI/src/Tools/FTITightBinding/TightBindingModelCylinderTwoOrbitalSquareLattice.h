////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//  class of tight binding model for the two-orbital model on square lattice  //
//          with periodic boundary conditions only along the x direction      //
//                                                                            //
//                        last modification : 27/10/2013                      //
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


#ifndef TIGHTBINDINGMODELCYLINDERTWOORBITALSQUARELATTICE_H
#define TIGHTBINDINGMODELCYLINDERTWOORBITALSQUARELATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract1DTightBindingModel.h"


class TightBindingModelCylinderTwoOrbitalSquareLattice : public Abstract1DTightBindingModel
{

 protected:
  
   // number of sites in the y direction
  int NbrSiteY;

  // imag part of the inter-orbital hopping amplitude between nearest neighbors along the x direction
  double NNHoppingInterX;
  // the inter-orbital hopping amplitude between nearest neighbors along the y direction
  double NNHoppingInterY;
  // the intra-orbital hopping amplitude between nearest neighbors
  double NNHoppingIntra;
  // folding factor for the momenta along sigma_x and sigma_y
  double FoldingFactor;  
  // four times the sublattice staggered chemical potential 
  double MuS;

public:

    // default constructor
    //
    // nbrSiteX = number of sites in the x direction
    // nbrSiteY = number of sites in the y direction
    // t1 = imag part of the inter-orbital hopping amplitude between nearest neighbors along the x direction
    // t2 = the inter-orbital hopping amplitude between nearest neighbors along the y direction
    // t3 = the intra-orbital hopping amplitude between nearest neighbors
    // foldingFactor = folding factor for the momenta along sigma_x and sigma_y
    // mus = sublattice chemical potential on A sites
    // gammaX = boundary condition twisting angle along x
    // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
    TightBindingModelCylinderTwoOrbitalSquareLattice(int nbrSiteX, int nbrSiteY, int t1, int t2, int t3, int foldingFactor, 
						     double mus, double gammaX, 
						     AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

    // destructor
    //
    ~TightBindingModelCylinderTwoOrbitalSquareLattice();

    // compute the many-body real space entanglement spectrum of a full band
    // 
    // outputFile = name of the output file where the spectrum has to be stored
    // minEnergy = lowest energy of the full band
    // maxEnergy = highest energy of the full band
    // nbrSiteYA = number of site to keep for the A part along the y direction    
    virtual void ComputeManyBodyRealSpaceEntanglementSpectrum(char* outputFile, double minEnergy, double maxEnergy, int nbrSiteYA);

protected :

    // compute the band structure
    //
    // minStateIndex = minimum index of the state to compute
    // nbrStates = number of states to compute
    virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif

