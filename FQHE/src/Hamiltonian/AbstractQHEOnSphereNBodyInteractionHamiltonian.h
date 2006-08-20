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
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"

#include <iostream>


using std::ostream;


class AbstractQHEOnSphereNBodyInteractionHamiltonian : public AbstractQHEOnSphereHamiltonian
{

 protected:
  
  // maximum number of particles that can interact together
  int MaxNBody;
  // indicates which n-body interaction terms are present in the Hamiltonian
  bool* NBodyFlags;

  // maximum value that the sum of indices can reach per n-body interaction
  int* MaxSumIndices;
  // minimum value that the sum of indices can reach per n-body interaction
  int* MinSumIndices;
  // array containing all interaction factors per N body interaction (first index for interaction type and second for index sum)
  double*** NBodyInteractionFactors;
  // array containing the number of index group per interaction (first index)and per index sum (second indices) for the creation (or annhilation) operators
  int** NbrSortedIndicesPerSum;
  // array containing the index group per interaction (first index)and per index sum (second indices) for the creation (or annhilation) operators
  int*** SortedIndicesPerSum;
  // sign in front of each n-body interaction term (including the one coming from the statistics)
  double* NBodySign;

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

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual  RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
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

  // enable fast multiplication algorithm
  //
  // fileName = prefix of the name of the file where temporary matrix elements will be stored
  virtual void EnableFastMultiplicationWithDiskStorage(char* fileName);

  // get all indices needed to characterize a completly skew symmetric tensor, ordered by the sum of the indices
  //
  // nbrValues = number of different values an index can have
  // nbrIndices = number of indices 
  // nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
  // sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
  // return value = total number of index groups
  virtual long GetAllSkewSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum);

  // get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
  //
  // nbrValues = number of different values an index can have
  // nbrIndices = number of indices 
  // nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
  // sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
  // sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
  //                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
  // return value = total number of index groups
  virtual long GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum,
				       double**& sortedIndicesPerSumSymmetryFactor);

};

#endif
