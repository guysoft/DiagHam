////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                   Class author: Cecile Repellin                            //
//                                                                            //
//                                                                            //
//         class of particle on lattice projected density operator            //
//                                                                            //
//                        last modification : 17/10/2013                      //
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



#ifndef PARTICLEONLATTICEPROJECTEDDENSITYOPERATOR_H
#define PARTICLEONLATTICEPROJECTEDDENSITYOPERATOR_H


#include "config.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"

#include <iostream>



class ParticleOnLatticeProjectedDensityOperator : public AbstractOperator
{

protected:
  // hilbert space associated to the particles of the original space
  ParticleOnSphere* ParticleSource;
  
  // hilbert space associated to the particles of the destination space
  ParticleOnSphere* ParticleDestination;
  
  // tight binding model
  Abstract2DTightBindingModel* TightBindingModel;
  
  //total momentum in x-direction of the density operator
  int Kx;
  //total momentum in y-direction of the density operator
  int Ky;
  
  
  
public:
  
  // default constructor
  //
  // particleSource = Hilbert space associated to the original system
  // particleDestination = Hilbert space associated to the destination system
  // tightBindingModel = tight binding model
  // kx = total momentum in x-direction of the density operator
  // ky = total momentum in y-direction of the density operator
  ParticleOnLatticeProjectedDensityOperator(ParticleOnSphere* particleSource, ParticleOnSphere* particleDestination, Abstract2DTightBindingModel* tightBindingModel, int kx, int ky);
  
  // destructor
  //
  ~ParticleOnLatticeProjectedDensityOperator();
  
   // clone operator (without duplicating datas)
  //
  // return value = pointer to cloned operator
  AbstractOperator* Clone();
  
  
  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);
  
  
  // get Hilbert space on which operator acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();
  
  // return dimension of Hilbert space where operator acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();
  
  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored

  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								    int firstComponent, int nbrComponent);
  
};

#endif