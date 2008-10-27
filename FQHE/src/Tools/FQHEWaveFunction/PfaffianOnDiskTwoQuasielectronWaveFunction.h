////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of Pfaffian wave function with two quasielectrons on disk      //
//                                                                            //
//                        last modification : 23/10/2008                      //
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


#ifndef PFAFFIANONDISKTWOQUASIELECTRONWAVEFUNCTION_H
#define PFAFFIANONDISKTWOQUASIELECTRONWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "GeneralTools/GarbageFlag.h"


class PfaffianOnDiskTwoQuasielectronWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // number of particles
  int NbrParticles;

  // position of the first quasielectron
  Complex ZElectron1;
  // position of the second quasielectron
  Complex ZElectron2;
  // gaussian weight associated to the quasielectrons
  double GaussianWeight;

  // Flag for bosons/fermions
  bool FermionFlag;

  // temporary array where the Pfaffian has to be stored
  Complex** TmpPfaffian;
  // temporary array where indices are stored
  int* TmpIndexArray;
  // temporary array used to store gaussian weights
  Complex* TmpGaussianWeights;

  // array containing description of each permutation that appears in the calculation symmetrization process
  unsigned long* Permutations;
  // number of permutations that appears in the symmetrization process
  unsigned long NbrPermutations;
  // garable flag associated to the Permutations array
  GarbageFlag Flag;


 public:

  // constructor
  //
  // nbrParticles = number of particles
  // zElectron1 = position of the first quasielectron
  // zElectron2 = position of the second quasielectron (spherical coordinates, theta angle)
  // fermions = flag indicating whether to calculate bosonic or fermionic pfaffian
  PfaffianOnDiskTwoQuasielectronWaveFunction(int nbrParticles, Complex zElectron1, Complex zElectron2, bool fermions=false);

  // constructor from data file 
  //
  // filename = pointer to the file name that described the symmetrization procedure
  // zElectron1 = position of the first quasielectron
  // zElectron2 = position of the second quasielectron (spherical coordinates, theta angle)
  // fermions = flag indicating whether to calculate bosonic or fermionic pfaffian
  PfaffianOnDiskTwoQuasielectronWaveFunction(char* filename, Complex zElectron1, Complex zElectron2, bool fermions=false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  PfaffianOnDiskTwoQuasielectronWaveFunction(const PfaffianOnDiskTwoQuasielectronWaveFunction& function);

  // destructor
  //
   ~PfaffianOnDiskTwoQuasielectronWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

  // write all permutations requested to symmetrize the state to data file 
  //
  // filename = pointer to the file name that described the symmetrization procedure
  // return value = true if no error occured
  virtual bool WritePermutations(char* filename);

 protected:


  // get all permutations requested to symmetrize the state from data file 
  //
  // filename = pointer to the file name that described the symmetrization procedure
  // return value = true if no error occured
  virtual bool ReadPermutations(char* filename);

  // evaluate all permutations requested to symmetrize the state
  //
  virtual void EvaluatePermutations();


};

#endif
