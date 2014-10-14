////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of bosons on lattice                          //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 11/09/2014                      //
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


#ifndef BOSONONLATTICEGUTZWILLERPROJECTIONONEORBITALPERSITEREALSPACEAND2DTRANSLATION_H
#define BOSONONLATTICEGUTZWILLERPROJECTIONONEORBITALPERSITEREALSPACEAND2DTRANSLATION_H

#include "config.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include <iostream>



class BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation : public  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation
{

 protected:
   
   int Lx;
   int Ly;
	
 public:

  // default constructor		
  // 
  BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation ();

  // basic constructor
  // 
  // nbrBosons = number of fermions
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation (int nbrBosons, int lx, int ly, int xMomentum, int maxXMomentum,
						 int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation (const  BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation& bosons);
  
  // destructor
  //
  ~BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation & operator = (const BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation & bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // compute sum of positions in the x and y direction for lattice class
  //
  // index = index of the state in the basis whose position sums are to be computed
  // positionX = reference on the sum of positions in the x direction
  // positionY = reference on the sum of positions in the y direction
  virtual void GetPositionSum(int index,int & positionX, int & positionY);
  virtual void GetPositionSum(unsigned long * monomial, int & positionX, int & positionY);

  virtual void GetCompositeFermionWavefunction(ComplexVector & trialState, ComplexMatrix & jastrowEigenVecs,ComplexMatrix & cFEigenVecs, double phaseTranslationX);
};

// compute sum of positions in the x and y direction for lattice class
//
// index = index of the state in the basis whose position sums are to be computed
// positionX = reference on the sum of positions in the x direction
// positionY = reference on the sum of positions in the y direction

inline void BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::GetPositionSum(int index, int & positionX, int & positionY)
{
  positionX =0;
  positionY =0;
  int Nx = this->Lx / this->MaxXMomentum;
  int Ny = this->Ly / this->MaxYMomentum;
  
  for (int i = 0; i < this->NbrSite; i++)
    { 	
         if (this->StateDescription[index]>> i & 1ul == 1 )
{
	 positionX +=  (i/this->StateXShift) *  Nx;
         positionY += ((i%this->StateXShift)/this->StateYShift) * Ny;
         positionX += (i % Nx);
         positionY += ((i%this->StateYShift)/ Nx) ;
}
    }
}

// compute sum of positions in the x and y direction for lattice class
//
// index = index of the state in the basis whose position sums are to be computed
// positionX = reference on the sum of positions in the x direction
// positionY = reference on the sum of positions in the y direction

inline void BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::GetPositionSum(unsigned long * monomial, int & positionX, int & positionY)
{
  positionX =0;
  positionY =0;
  int Nx = this->Lx / this->MaxXMomentum;
  int Ny = this->Ly / this->MaxYMomentum;
  
  for (int i = 0; i < this->NbrBosons; i++)
    { 	
	 positionX +=  ( monomial[i]/this->StateXShift) *  Nx;
         positionY += ((monomial[i]%this->StateXShift)/this->StateYShift) * Ny;
         positionX += (monomial[i] % Nx);
         positionY += ((monomial[i]%this->StateYShift)/ Nx) ;
    }
}

#endif
