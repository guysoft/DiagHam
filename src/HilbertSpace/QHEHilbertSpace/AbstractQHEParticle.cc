////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  abstract class of particle in QHE systems                 //
//                                                                            //
//                        last modification : 10/10/2004                      //
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


#include "config.h"
#include "HilbertSpace/QHEHilbertSpace/AbstractQHEParticle.h"


// virtual destructor
//

AbstractQHEParticle::~AbstractQHEParticle ()
{
}

// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex AbstractQHEParticle::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// return value = wave function evaluated at the given location

Complex AbstractQHEParticle::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, 
						     this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location
Complex AbstractQHEParticle::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
						int firstComponent, int nbrComponent)
{
  return Complex(0.0, 0.0);
}


// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex AbstractQHEParticle::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, 
								 int nextCoordinates, int firstComponent, 
								 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void AbstractQHEParticle::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
                                    
