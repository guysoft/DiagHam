////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of U(1) wave function obtained form symmetrization of a   //
//                        SU(K) wave function on sphere                       //
//                                                                            //
//                        last modification : 17/03/2008                      //
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


#ifndef FQHESPHERESYMMETRIZEDSUKTOU1WAVEFUNCTION_H
#define FQHESPHERESYMMETRIZEDSUKTOU1WAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "GeneralTools/GarbageFlag.h"


class FQHESphereSymmetrizedSUKToU1WaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles
  int NbrParticles;

  // number of particles per color
  int NbrParticlesPerColor;

  // K index (i.e. SU(K))
  int KValue;

  // true if the final state should be a fermionic state
  bool FermionFlag;  

  // internal storage of spinor coordinates
  Complex* SpinorUCoordinates;
  Complex* SpinorVCoordinates;

  // pointer to the base SU(K) wave function
  Abstract1DComplexFunctionOnSphere* SUKWaveFunction;
  
  // array containing description of each permutation that appears in the calculation symmetrization process
  unsigned long** Permutations;
  // number of permutations that appears in the symmetrization process
  unsigned long NbrPermutations;

  unsigned long* ColorPermutations;
  // garable flag associated to the Permutations array
  GarbageFlag Flag;

 public:

  // constructor
  //
  // nbrParticles = number of particles
  // kValue = number of particle per cluster
  // sUKWavefunction = pointer to the base SU(K) wave function
  // fermionFlag = true if the final state should be a fermionic state
  FQHESphereSymmetrizedSUKToU1WaveFunction(int nbrParticles, int kValue, Abstract1DComplexFunctionOnSphere* sUKWavefunction, bool fermionFlag = false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  FQHESphereSymmetrizedSUKToU1WaveFunction(const FQHESphereSymmetrizedSUKToU1WaveFunction& function);

  // destructor
  //
  ~FQHESphereSymmetrizedSUKToU1WaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

  // evaluate function at a given point
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  // return value = function value at (uv)
  Complex CalculateFromSpinorVariables(ComplexVector& uv);

 private:

  // evaluate all permutations requested to symmetrize the SU(K) state
  //
  void EvaluatePermutations();

};

#endif
