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


#include "SimpleMonteCarloOnSphere.h"

#include "Options/Options.h"

#include <iostream>
using std::cout;
using std::endl;

// default constructor
SimpleMonteCarloOnSphere::SimpleMonteCarloOnSphere()
{
  this->NbrParticles=0;
}

// set up for basic monte-carlo scheme
// nbrParticles = number of particles in system
// waveFunction = wavefunction to be simulated
// samplingFunction = function to be used to generate samples
SimpleMonteCarloOnSphere::SimpleMonteCarloOnSphere(int nbrParticles, Abstract1DComplexFunction *waveFunction,
			 AbstractMCSamplingFunction *samplingFunction)
{
  this->NbrParticles=nbrParticles;
  this->WaveFunction=waveFunction; // is returned normalized by wavefunction handler...
  this->SamplingFunction=samplingFunction;
  this->System=new ParticleOnSphereCollection(this->NbrParticles);
  this->SamplingFunction->RegisterSystem(System);
  this->SamplingFunction->AdaptAverageMCNorm(); // this also relaxes the particle positions in System
}
  
// destructor
SimpleMonteCarloOnSphere::~SimpleMonteCarloOnSphere()
{
  if (this->NbrParticles!=0)
    delete this->System;  
}


// add an option group containing all options related to the wave functions
//
// manager = pointer to the option manager
void SimpleMonteCarloOnSphere::AddOptionGroup(OptionManager* manager)
{
 this->Options = manager;
 OptionGroup* MCGroup  = new OptionGroup ("Monte-Carlo options");
  (*(this->Options)) += MCGroup;
  
  (*MCGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MCGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MCGroup) += new SingleIntegerOption  ('\n', "randomSeed", "random seed for internal random number generator", -1);
    
}
