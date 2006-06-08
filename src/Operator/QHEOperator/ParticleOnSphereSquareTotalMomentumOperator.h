////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of particle on sphere sqare of total momentum operator         //
//                                                                            //
//                        last modification : 06/03/2005                      //
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


#ifndef PARTICLEONSPHERESQUARETOTALMOMENTUMOPERATOR_H
#define PARTICLEONSPHERESQUARETOTALMOMENTUMOPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/QHEHilbertSpace/ParticleOnSphere.h"
#include "Matrix/RealMatrix.h"


class ParticleOnSphereSquareTotalMomentumOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  ParticleOnSphere* Particle;

  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a particle
  int LzMax;

  // matrix where all coefficents that come from the L+L- terms are stored
  RealMatrix Coefficients;
  // shift due to the Lz^2 term
  double Shift;

 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a fermion
  ParticleOnSphereSquareTotalMomentumOperator(ParticleOnSphere* particle, int totalLz, int lzMax);

  // copy constructor
  //
  // oper = reference on the operator to copy
  ParticleOnSphereSquareTotalMomentumOperator(const ParticleOnSphereSquareTotalMomentumOperator& oper);

  // destructor
  //
  ~ParticleOnSphereSquareTotalMomentumOperator();
  
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
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);
   
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
