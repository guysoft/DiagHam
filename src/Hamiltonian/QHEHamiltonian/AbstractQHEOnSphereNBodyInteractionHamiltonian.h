////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract quantum Hall hamiltonian associated          //
//            to particles on a sphere with n-body interaction terms          //
//                                                                            //
//                        last modification : 22/09/2004                      //
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


#ifndef ABSTRACTQHEONSPHERENBODYINTERACTIONHAMILTONIAN_H
#define ABSTRACTQHEONSPHERENBODYINTERACTIONHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/QHEHilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/QHEHamiltonian/AbstractQHEOnSphereHamiltonian.h"

#include <iostream>


using std::ostream;


class AbstractQHEOnSphereNBodyInteractionHamiltonian : public AbstractQHEOnSphereHamiltonian
{

 protected:
  
  // maximum number of particles that can interact together
  int MaxNBody;
  // indicates which n-body interaction terms are present in the Hamiltonian
  bool* NBodyFlags;

  // array containing all interaction factors per N body interaction
  double** NBodyInteractionFactors;
  // number of interaction factors per N body interaction
  long* NBodyNbrInteractionFactors;
  // arrays for creation operartor indices attached to each interaction factor and per N body interaction
  int*** NBodyMValue;
  // arrays for annihilation operartor indices attached to each interaction factor and per N body interaction
  int*** NBodyNValue;
  
 public:

  // destructor
  //
  virtual ~AbstractQHEOnSphereNBodyInteractionHamiltonian() = 0;

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					  int firstComponent, int nbrComponent);

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication();

  // enable fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int lastComponent);

};

#endif
