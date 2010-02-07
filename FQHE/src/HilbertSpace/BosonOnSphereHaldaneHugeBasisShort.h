////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on sphere  using the Haldane basis           //
//                            for system size such that                       //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                          work for huge Hilbert space                       //
//                                                                            //
//                        last modification : 08/02/2009                      //
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


#ifndef BOSONONSPHEREHALDANEHUGEBASISSHORT_H
#define BOSONONSPHEREHALDANEHUGEBASISSHORT_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"

#include <iostream>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class AbstractArchitecture;


class BosonOnSphereHaldaneHugeBasisShort :  public BosonOnSphereShort
{

 protected:

  // the fermionic huge Hilbert space associated to the bosonic one
  FermionOnSphereHaldaneHugeBasis* FermionHugeBasis;

  // temporary array to store monomial representation
  unsigned long* TemporaryMonomial;
  unsigned long* TemporaryMonomial2;

  friend class FQHESphereJackGeneratorOperation;

 public:

  // default constructor
  //
  BosonOnSphereHaldaneHugeBasisShort ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  // maxFileSize = maximum file size (in MBytes)
  // referenceState = array that describes the reference state to start from
  // memory = amount of memory granted for precalculations
  // symmetricFlag = indicate if a symmetric basis has to be used (only available if totalLz = 0)
  // fullDimension = provide the full (i.e. without squeezing) Hilbert space dimension (0 if it has to be computed)
  BosonOnSphereHaldaneHugeBasisShort (int nbrBosons, int totalLz, int lzMax, unsigned long maxFileSize, int* referenceState, unsigned long memory = 10000000, bool symmetricFlag = false, long fullDimension = 0l);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memoryHilbert = amount of memory granted to store the Hilbert space (in Mbytes)
  BosonOnSphereHaldaneHugeBasisShort (char* fileName, unsigned long memoryHilbert);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereHaldaneHugeBasisShort(const BosonOnSphereHaldaneHugeBasisShort& bosons);

  // destructor
  //
  virtual ~BosonOnSphereHaldaneHugeBasisShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereHaldaneHugeBasisShort& operator = (const BosonOnSphereHaldaneHugeBasisShort& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // check if disk is used to store the Hilbert space
  //
  // return value = true if disk storage is used
  bool CheckDiskStorage();

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState);

  // create the Jack polynomial decomposition corresponding to the root partition
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized hugebasis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // partialSave = save partial results in a given vector file
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized hugebasis
  RealVector& GenerateJackPolynomial(RealVector& jack, double alpha, long minIndex = 0l, long maxIndex = 0l, char* partialSave = 0);

  // create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized hugebasis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // partialSave = save partial results in a given vector file
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized hugebasis
  RealVector& GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha, long minIndex = 0l, long maxIndex = 0l, char* partialSave = 0);

  // create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry and using sparse storage
  //
  // alpha = value of the Jack polynomial alpha coefficient
  // architecture = architecture to use for precalculation
  // partialSave = save partial results in a given vector file
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // memory = amount of memory (in bytes) allowed for temporary vector storage (0 if 
  // memoryBlock = amount of memory (in bytes) allowed for precomputing state indices  (0 if the whole vector has to be stored in memory)
  virtual void GenerateSymmetrizedJackPolynomialSparse(double alpha, AbstractArchitecture* architecture, char* partialSave, long minIndex = 0l, long maxIndex = 0l, long memory = 0l, long memoryBlock = 0l);

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has to be normalized to 1
  // symmetryFactor = if true also remove the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // convert a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // symmetryFactor = if true also add the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // fuse two states which belong to different Hilbert spaces 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // leftVector = reference on the vector whose Hilbert space will be fuse to the left
  // rightVector = reference on the vector whose Hilbert space will be fuse to the right
  // padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
  // leftSpace = point to the Hilbert space that will be fuse to the left
  // rightSpace = point to the Hilbert space that will be fuse to the right
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // return value = reference on the fused state
  RealVector& FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
			  ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
			  bool symmetrizedFlag);

  // use product rule to produce part of the components of a system from a smaller one
  //
  // outputVector = reference on the vector which will contain the product rule state  (without zeroing components which do not occur in the fusion)
  // inputVector = reference on the vector associated to the smaller system
  // inputSpace = pointer to the Hilbert space of the smaller system
  // commonPattern = array describing the shared leftmost pattern between the n-body states in both the smaller and larger system sizes
  // commonPatterSize = number of elements in the commonPattern array
  // addedPattern = array describing the pattern that has to be inserted to go from the smaller system to the larger one
  // addedPatterSize = number of elements in the addedPattern array
  // coefficient = multiplicqtive fqctor to go fron the component of the smaller system to the larger one
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // return value = reference on the product rule state
  RealVector& ProductRules (RealVector& outputVector, RealVector& inputVector, ParticleOnSphere* inputSpace, 
			    int* commonPattern, int commonPatterSize, int* addedPattern, int addedPatterSize,
			    double coefficient, bool symmetrizedFlag);

  // compute part of the Jack polynomial square normalization in a given range of indices
  //
  // state = reference on the unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  virtual double JackSqrNormalization (RealVector& outputVector, long minIndex = 0l, long nbrComponents = 0l);

 protected :
  
  // core part of the Jack generator using the Lz<->-Lz symmetry and the factorized algorithm
  //
  // invAlpha = inverse of the Jack polynomial alpha coefficient
  // maxRoot = root partition (in fermionic binary representation)
  // partialSave = save partial results in a given vector file
  // minIndex = start computing the Jack polynomial from the minIndex-th component
  // maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
  // indexArray = array where state indices are stored
  // stateArray = array use to store computed state description
  // componentArray = array where computed component numerical factors are stored
  // nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
  // rhoArray = rho factor associated to each state
  virtual void GenerateSymmetrizedJackPolynomialFactorizedCore(double invAlpha, unsigned long maxRoot, long minIndex, long maxIndex, unsigned long** stateArray, double** componentArray, long** indexArray, int* nbrComputedComponents, double* rhoArray);

};

// check if disk is used to store the Hilbert space
//
// return value = true if disk storage is used

inline  bool BosonOnSphereHaldaneHugeBasisShort::CheckDiskStorage()
{
  return this->FermionHugeBasis->CheckDiskStorage();
}
#endif
