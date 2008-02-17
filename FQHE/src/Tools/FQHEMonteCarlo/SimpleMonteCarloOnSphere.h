////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//      class for a basic Monte Carlo algorith for particles on a sphere      //
//                                                                            //
//                        last modification : 23/01/2008                      //
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


#ifndef SIMPLEMONTECARLOONSPHERE_H
#define SIMPLEMONTECARLOONSPHERE_H


#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "AbstractMCSamplingFunction.h"
#include "Tools/FQHEWaveFunction/ParticleOnSphereCollection.h"

class OptionManager;

class SimpleMonteCarloOnSphere
{
 protected:

  // number of particles
  int NbrParticles;

  // wavefunction to simulate
  Abstract1DComplexFunction *WaveFunction;

  // sampling function
  AbstractMCSamplingFunction *SamplingFunction;  
  
  // class holding the particle coordinates
  ParticleOnSphereCollection *System;

  // pointer to the option manager
  OptionManager* Options;
  
  
 public:

  // default constructor
  SimpleMonteCarloOnSphere();

  // set up for basic monte-carlo scheme
  // nbrParticles = number of particles in system
  // waveFunction = wavefunction to be simulated
  // samplingFunction = function to be used to generate samples
  SimpleMonteCarloOnSphere(int nbrParticles, Abstract1DComplexFunction *waveFunction,
			   AbstractMCSamplingFunction *samplingFunction);
  
  // destructor
  ~SimpleMonteCarloOnSphere();

  // add an option group containing all options related to the wave functions
  //
  // manager = pointer to the option manager
  void AddOptionGroup(OptionManager* manager);

  
  
};

#endif // SIMPLEMONTECARLOONSPHERE_H
