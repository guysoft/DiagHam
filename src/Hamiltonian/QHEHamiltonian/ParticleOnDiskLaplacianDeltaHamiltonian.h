////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated to particles on a disk with         //
//                          laplacian delta interaction                       //
//                                                                            //
//                        last modification : 05/02/2004                      //
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


#ifndef PARTICLEONDISKDELTAHAMILTONIAN_H
#define PARTICLEONDISKDELTAHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/QHEHilbertSpace/ParticleOnDisk.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/AbstractQHEOnDiskHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnDiskLaplacianDeltaHamiltonian : public AbstractQHEOnDiskHamiltonian
{

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  ParticleOnDiskLaplacianDeltaHamiltonian(ParticleOnDisk* particles, int nbrParticles);

  // destructor
  //
  ~ParticleOnDiskLaplacianDeltaHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);

  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, ParticleOnDiskLaplacianDeltaHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnDiskLaplacianDeltaHamiltonian& H);

 private:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

};

#endif
