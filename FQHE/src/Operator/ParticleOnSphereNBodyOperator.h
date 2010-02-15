////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of particle on sphere n-body operator               //
//                                                                            //
//                        last modification : 13/12/2004                      //
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


#ifndef PARTICLEONSPHERENBODYOPERATOR_H
#define PARTICLEONSPHERENBODYOPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnSphere.h"


class ParticleOnSphereNBodyOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  ParticleOnSphere* Particle;
  
  // array containing the indices of the creation operators from left to right (a+_0.....a+_(n-1))
  int* CreationIndices;
  // array containing the indices of the annihilation operators from left to right (a_0.....a_(n-1))
  int* AnnihilationIndices;
  // number of creation (or annihilation) operators
  int NbrNBody;
  
 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // creationIndices = array containing the indices of the creation operators from left to right (a+_0.....a+_(n-1))
  // annihilationIndices = rray containing the indices of the annihilation operators from left to right (a_0.....a_(n-1))
  // nbrNBody = number of creation (or annihilation) operators
  ParticleOnSphereNBodyOperator(ParticleOnSphere* particle, int* creationIndices,
				int* annihilationIndices, int nbrNBody);

  // copy constructor
  //
  // oper = reference on the operator to copy
  ParticleOnSphereNBodyOperator(const ParticleOnSphereNBodyOperator& oper);

  // destructor
  //
  ~ParticleOnSphereNBodyOperator();
  
  // clone operator without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractOperator* Clone ();

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
  
  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  Complex PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent);

  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
			       int firstComponent, int nbrComponent);
  
};

#endif
