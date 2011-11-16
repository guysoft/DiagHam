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


#ifndef MRBLOCKSAMPLINGFUNCTION_H
#define MRBLOCKSAMPLINGFUNCTION_H

#define TESTING_MRBLOCKSAMPLINGFUNCTION_H

#include "AbstractMCBlockSamplingFunction.h"

class MRBlockSamplingFunction : public AbstractMCBlockSamplingFunction
{
 protected:
  // total number of particles
  int NbrParticles;
  // number of particles in each layer
  int NbrParticlesPerBlock;

  // pointers to spinor coordinates (external)
  Complex *SpinorUCoordinates;
  Complex *SpinorVCoordinates;

  // for access to ParticleCollection
  Complex LastU;
  Complex LastV;

  // norm for total function value
  double ElementNorm;

  // number of distinct blocks
  int NbrBlocks;
  // array for temporary storage of indices
  int **BlockPermutations;  

  // weights associated with each of the blocks
  double *BlockWeights;

  // table for JastrowElements
  Complex *JastrowElements;

  // critical distance for evaluation
  double SqrCriticalDistance;

  // flags for critical distances
  unsigned long *CriticalParticles;

  bool CriticalFlag;

  Complex JastrowFactor;
  
 public:
  // constructor
  // nbrParticlesPerBlock = N/2
  // sqrCriticalDistance = tuning parameter for recalculation of blocks
  MRBlockSamplingFunction(int nbrParticlesPerBlock, double sqrCriticalDistance = 1e-20);
  // virtual destructor
  virtual ~MRBlockSamplingFunction();

  // register basic system of particles
  virtual void RegisterSystem(AbstractParticleCollection *system);

  // method for ratio of probabilities with respect to the last configuration
  // allows for more rapid calculation due to cancellation of factors
  virtual double GetTransitionRatio();

  // get the full function value for a system of particles
  virtual Complex GetFunctionValue();

    // get the estimate of the sampling function value calculated over the sampling block only, for the given system of particles
  virtual Complex GetSamplingBlockValue();

  // get the current phase of the sampling function
  Complex GetSamplingBlockPhase();

  // call this method to scale the sampling function (needed to normalize the function)
  // scale = total scaling factor
  virtual void ScaleByFactor(double scale);

   // accessor routine for NbrBlocks
  virtual int GetNbrBlocks() {return this->NbrBlocks;}

    // query weights of individual blocks
  // weights = array reference to values of weighting factors
  virtual void GetBlockWeights(double *weights);

  // get the Monte Carlo amplitude for the requested block with nbrPermute particles exchanged
  // nbrBlock = nbrPermute = number of particles to exchange between blocks
  // amplitude = reference of return argument
  virtual void GetBlockAmplitude(int nbrBlock, Complex &amplitude);

  // get the Monte Carlo amplitude for the requested block with nbrPermute particles exchanged
  // amplitudes = pointer to array of return arguments
  virtual void GetAllBlockAmplitudes(Complex *amplitudes);

  // query permutation of particles applied to given block
  // nbrBlock = block index
  // permutations = pointer to array to be filled with return argument
  virtual void GetBlockPermutations(int nbrBlock, int* permutations);

  // query permutation of particles applied to given block
  // nbrBlock = block index
  // permutations = pointer to array to be filled with return argument
  virtual void GetAllBlockPermutations(int** permutations);

 protected:

  // get linearized index for symmetric array
  inline int SymIndex(int i,int j);

  // precalculate Jastrow factor elements
  void EvaluateTable();
  
  
};

#ifdef TESTING_MRBLOCKSAMPLINGFUNCTION_H
#include <iostream>
using std::cout;
using std::endl;
#endif

inline int MRBlockSamplingFunction::SymIndex(int i, int j)
{
#ifdef TESTING_MRBLOCKSAMPLINGFUNCTION_H
  if (i<j)
#endif
    return (this->NbrParticles*i-i*(i-1)/2+j-i);
#ifdef TESTING_MRBLOCKSAMPLINGFUNCTION_H
  else
    {
      cout << "Attention, need sign IN SymIndex!"<<endl;
      //return (this->NbrParticles*j-j*(j-1)/2+i-j);
      return 0;
    }
#endif
}

#endif 
