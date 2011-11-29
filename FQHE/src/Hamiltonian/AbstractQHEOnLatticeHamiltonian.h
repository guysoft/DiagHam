////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                 class of quatum Hall hamiltonian associated                //
//                 to particles on a lattice in magnetic field                //
//                                                                            //
//                        last modification : 12/02/2008                      //
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


#ifndef ABSTRACTQHEONLATTICEHAMILTONIAN_H
#define ABSTRACTQHEONLATTICEHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"
#include "GeneralTools/RealUniqueArray.h"
#include "GeneralTools/ComplexUniqueArray.h"

#include <iostream>
#include <cstdlib>
#include <climits>


using std::ostream;


class AbstractArchitecture;

using std::cout;
using std::endl;

// threshold before something is defined different from zero
#define LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD 1.71245451764e-13


class AbstractQHEOnLatticeHamiltonian : public AbstractQHEHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:
  
  // Hilbert space associated to the system
  ParticleOnLattice* Particles;

  // number of particles
  int NbrParticles;

  // lattice dimension in the x-direction
  int Lx;
  // lattice dimension in the y-direction
  int Ly;
  // number of sublattice sites per unit cell
  int SubLattices;
  // ky-symmetry-flag
  bool HaveKySymmetry;
  // maximum ky-momentum
  int KyMax;  
  // number of sites
  int NbrSites;  
  // number of unit cells
  int NbrCells;
  // number of flux through lattice
  int NbrFluxQuanta;
  // flux through a unit cell
  double FluxDensity;

  // array containing kinetic energy (hopping) terms
  Complex* HoppingTerms;
  // number of hopping terms
  int NbrHoppingTerms;
  // arrays for indices attached to hopping terms (i nitial and f inal)
  int* KineticQi;
  int* KineticQf;

  // fields for general four-fermion interactions:
  // array containing interaction factors
  Complex* InteractionFactors;
  // number of interaction factors
  int NbrInteractionFactors;
  // arrays for quantum numbers attached to each interaction factor
  int* Q1Value;
  int* Q2Value;
  int* Q3Value;
  int* Q4Value;
  
  // alternative method used to store indices attached to each interaction factor
  // Q1Value and Q2Value still store the NbrQ12Indices different combinations of Q1/Q2
  // InteractionFactors still hold the matrix elements in a completely linearized array
  int NbrQ12Indices;
  // number of different combinations of Q3/Q4 for each pair of Q1/Q2
  int* NbrQ34Values;
  // and the corresponding values of Q3/Q4 themselves
  int** Q3PerQ12;
  int** Q4PerQ12;

  // separate out completely diagonal contact-interactions:
  double* DiagonalInteractionFactors;
  // number of interaction factors (usually = number of sites)
  int NbrDiagonalInteractionFactors;
  // quantum number attached to each "self-interaction" factor
  int* DiagonalQValues;

  // flag for storing the result of test for complex matrix elements
  bool HaveComplexMatrixElements;
  // flag whether test has been carried out
  bool HaveTestedForComplexMatrixElement;

  // shift to apply to go from precalculation index to the corresponding index in the HilbertSpace
  int PrecalculationShift;
  
  // flag for implementation of hermitian symmetry
  bool HermitianSymmetryFlag;
  
  // amount of memory (in bytes) that can be used to store precalculated matrix elements
  long AllowedMemory;
  // flag indicating whether fast multiplication data was loaded from a file
  bool LoadedPrecalculation;
  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // step between each precalculated index (main part: start at 0, FastMultiplicationStep, 2*FastMultiplicationStep, ...)
  int FastMultiplicationStep;
  // number of non-null real terms in the hamiltonian for each state (typically a small number)
  unsigned short* NbrRealInteractionPerComponent;
  // number of non-null complex terms in the hamiltonian for each state
  unsigned short* NbrComplexInteractionPerComponent;
  // index of the state obtained for each term of the hamiltonian when applying on a given state
  // holding indices for both real (1st) and complex matrix elements
  int** InteractionPerComponentIndex;
  // index of real (first in enumeration) and complex (following) multiplicative coefficients obtained for each term of the hamiltonian when applying on a given state and with a given destination state
  unsigned short** InteractionPerComponentCoefficientIndex;
  

  // array storing a single copy of each real matrix element value
  RealUniqueArray RealInteractionCoefficients;

  // array storing a single copy of each real matrix element value
  ComplexUniqueArray ComplexInteractionCoefficients;

  // number of tasks for load balancing
  int NbrBalancedTasks;
  // load balancing array for parallelisation, indicating starting indices
  long *LoadBalancingArray;
  // cumulative count of non-zero matrix elements

  // flag to indicate if a hamiltonian is temporary stored on disk
  bool DiskStorageFlag;
  // name of the file that contains hamiltonian matrix elements
  char* DiskStorageFileName;
  // index of the first row that appears in the on-disk hamiltonian
  int DiskStorageStart;
  // maximum number of non-null terms in the hamiltonian for each state
  int MaxNbrInteractionPerComponent;
  // size of the in-memory temporary buffer
  long BufferSize;

  // shift to apply to the Hamiltonian diagonal elements
  double HamiltonianShift;



 public:

  // default constructor
  //
  AbstractQHEOnLatticeHamiltonian();

  // destructor
  //
  virtual ~AbstractQHEOnLatticeHamiltonian() = 0;

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // set flux density in units of flux quanta through the lattice
  //
  // nbrFluxQuanta = flux quantua piercing the lattice
  virtual void SetNbrFluxQuanta(int nbrFluxQuanta);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  virtual AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  virtual int GetHilbertSpaceDimension ();
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  virtual void ShiftHamiltonian (double shift);
  
  // ask if Hamiltonian implements methods using hermitian symmetry 
  //
  virtual bool IsHermitian();

  // ask if Hamiltonian implements methods applying the conjugate of the Hamiltonian
  //
  virtual bool IsConjugate();

  // ask if Hamiltonian has any complex matrix elements
  //
  virtual bool IsComplex();

  // symmetrize interaction factors to enable hermitian matrix multiplication
  // return = true upon success
  virtual bool HermitianSymmetrizeInteractionFactors();

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied<
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						     int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied<
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& ConjugateLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
  virtual ComplexVector* ConjugateLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							      int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied<
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							      int firstComponent, int nbrComponent);


    // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied<
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
  virtual RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						     int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied<
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
  virtual RealVector* ConjugateLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							      int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied<
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
  virtual RealVector* HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							      int firstComponent, int nbrComponent);

 
  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  virtual List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  virtual List<Matrix*> RightInteractionOperators();  

  // save precalculations in a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be stored
  // return value = true if no error occurs
  virtual bool SavePrecalculation (char* fileName);

 protected:
 
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								int firstComponent, int nbrComponent);
	
	// multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
	virtual ComplexVector& LowLevelAddMultiplyFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int lastComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
							int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								int firstComponent, int nbrComponent);
	
	// multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
	virtual ComplexVector* LowLevelMultipleAddMultiplyFastMultiply(ComplexVector* vSources, ComplexVector * vDestinations,  int nbrVectors, int firstComponent, int lastComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& ConjugateLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									 int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& ConjugateLowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
								 int firstComponent, int nbrComponent);
	
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& ConjugateLowLevelAddMultiplyFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int lastComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										 int firstComponent, int nbrComponent);
	
  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored  
  virtual ComplexVector* ConjugateLowLevelMultipleAddMultiplyFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									  int firstComponent, int lastComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* ConjugateLowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									 int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									 int firstComponent, int nbrComponent);
	
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiplyFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								  int firstComponent, int lastComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
								 int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent);
	
	  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiplyFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									 int firstComponent, int nbrComponent);

    // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
								int firstComponent, int nbrComponent);
	
	// multiply a vector by the current hamiltonian for a given range of indices 
	// and add result to another vector, low level function (no architecture optimization)
	// using partial fast multiply option
	//
	// vSource = vector to be multiplied
	// vDestination = vector at which result has to be added
	// firstComponent = index of the first component to evaluate
	// nbrComponent = number of components to evaluate
	// return value = reference on vector where result has been stored
	virtual RealVector& LowLevelAddMultiplyFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int lastComponent);
	
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
							int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent);
	
	// multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
	virtual RealVector* LowLevelMultipleAddMultiplyFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int lastComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
								int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& ConjugateLowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
									 int firstComponent, int nbrComponent);
	
	// multiply a vector by the current hamiltonian for a given range of indices 
	// and add result to another vector, low level function (no architecture optimization)
	// using partial fast multiply option
	//
	// vSource = vector to be multiplied
	// vDestination = vector at which result has to be added
	// firstComponent = index of the first component to evaluate
	// nbrComp	onent = number of components to evaluate
	// return value = reference on vector where result has been stored
	virtual RealVector& ConjugateLowLevelAddMultiplyFastMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int lastComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& ConjugateLowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
								 int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										 int firstComponent, int nbrComponent);
	
	// multiply a set of vectors by the current hamiltonian for a given range of indices 
	// and add result to another et of vectors, low level function (no architecture optimization)
	// using partial fast multiply option
	//
	// vSources = array of vectors to be multiplied
	// vDestinations = array of vectors at which result has to be added
	// nbrVectors = number of vectors that have to be evaluated together
	// firstComponent = index of the first component to evaluate
	// nbrComponent = number of components to evaluate
	// return value = pointer to the array of vectors where result has been stored
	virtual RealVector* ConjugateLowLevelMultipleAddMultiplyFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
																																			 int firstComponent, int lastComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* ConjugateLowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									 int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
									 int firstComponent, int nbrComponent);
	
	
	// multiply a vector by the current hamiltonian for a given range of indices 
	// and add result to another vector, low level function (no architecture optimization)
	// using partial fast multiply option
	//
	// vSource = vector to be multiplied
	// vDestination = vector at which result has to be added
	// firstComponent = index of the first component to evaluate
	// nbrComponent = number of components to evaluate
	// return value = reference on vector where result has been stored
	virtual RealVector& HermitianLowLevelAddMultiplyFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int lastComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
								 int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										 int firstComponent, int nbrComponent);
	
	 // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
	virtual RealVector* HermitianLowLevelMultipleAddMultiplyFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										   int firstComponent, int lastComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									 int firstComponent, int nbrComponent);

  
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

  // get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
  // nbrThreads = number of threads requested
  // segmentIndices = array returning the reference to an array of the first index of each of the segments
  //
  virtual bool GetLoadBalancing(int nbrTasks, long* &segmentIndices);

  // test the amount of memory needed for fast multiplication algorithm
  //
  // allowedMemory = amount of memory that cam be allocated for fast multiplication
  // return value = amount of memory needed
  virtual long FastMultiplicationMemory(long allowedMemory);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = number of components that have to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int nbrComponent);

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication();

  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int nbrComponent);
  
  // enable fast multiplication algorithm using on disk cache 
  //
  // fileName = prefix of the name of the file where temporary matrix elements will be stored
  virtual void EnableFastMultiplicationWithDiskStorage(char* fileName);

  // load precalculations from a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be read
  // return value = true if no error occurs
  virtual bool LoadPrecalculation (char* fileName);

  // core part of the FastMultiplication method involving 1- and 2-body terms and diagonal elements
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientIndexArray = array of the indices of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnLattice* particles, int index, 
					   int* indexArray, unsigned short* coefficientIndexArray,int& positionR, int & positionC);
	
	  // core part of the FastMultiplication method involving 1- and 2-body terms and diagonal elements
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientIndexArray = array of the indices of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnLattice* particles, int index, 
					   int* indexArray, unsigned short* coefficientIndexArray, int & positionR, int & positionC);
	
	// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector& vSource, ComplexVector& vDestination);
	
	// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

	// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector& vSource, RealVector& vDestination);
	
	// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors);

	// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector& vSource, ComplexVector& vDestination);
	
	// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

virtual void EvaluateMNOneBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

	// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector& vSource, RealVector& vDestination);
	
	// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

virtual void EvaluateMNOneBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors);
	
virtual void HermitianEvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector& vSource, ComplexVector& vDestination);

virtual void HermitianEvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector& vSource, RealVector& vDestination);

virtual void HermitianEvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector * vSources, ComplexVector* vDestinations, int nbrVectors);

virtual void HermitianEvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors);


	
	// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

	virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);
	
	virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector* vSources,  ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);
	
	virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector& vSource, RealVector& vDestination);
	
	virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector* vSources,  RealVector* vDestinations, int nbrVectors, double* tmpCoefficients);
	
	// core part of the AddMultiply method involving the two-body interaction
	// 
	// particles = pointer to the Hilbert space
	// index = index of the component on which the Hamiltonian has to act on
	// vSource = vector to be multiplied
	// vDestination = vector at which result has to be added
	
	virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);
	
	virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector* vSources,  ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);
	
	virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector& vSource, RealVector& vDestination);
	
	virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector* vSources,  RealVector* vDestinations, int nbrVectors, double* tmpCoefficients);

	
	virtual void EvaluateMNTwoBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector&  vSources,  ComplexVector&  vDestinations);
	
	
	virtual void EvaluateMNTwoBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector* vSources,  ComplexVector* vDestinations, int nbrVectors, Complex * tmpCoefficients);
	
	
	virtual void EvaluateMNTwoBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector&  vSources,  RealVector& vDestinations);
	
	
	virtual void EvaluateMNTwoBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector* vSources,  RealVector* vDestinations, int nbrVectors, double* tmpCoefficients);
	
	
	// core part of the PartialFastMultiplicationMemory method involving one-body term
	// 
	// particles = pointer to the Hilbert space
	// firstComponent = index of the first component that has to be precalcualted
	// lastComponent  = index of the last component that has to be precalcualted
	// memory = reference on the amount of memory required for precalculations
	
	virtual void EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, long& memory);
	
	
	// core part of the PartialFastMultiplicationMemory method involving two-body term
	// 
	// particles = pointer to the Hilbert space
	// firstComponent = index of the first component that has to be precalcualted
	// lastComponent  = index of the last component that has to be precalcualted
	// memory = reference on the amount of memory required for precalculations
	
	virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, long& memory);
	
};


// core part of the FastMultiplication method involving 1- and 2-body terms and diagonal elements
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientIndexArray = array of the indices of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnLattice* particles, int index, 
										 int* indexArray, unsigned short* coefficientIndexArray, int& positionR, int & positionC)
{
  int qi, qf;
  int Index2;
  int tmpElementPos;
  double Coefficient;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  // deal with kinetic energy terms first!             
  for (int j = 0; j < NbrHoppingTerms; ++j)
    {
      qi = this->KineticQi[j];
      qf = this->KineticQf[j];
      // considering: this->HoppingTerms[j			
      Index2 = particles->AdA(index, qf, qi, Coefficient);	  
      if (Index2 < Dim)
	{
	  //cout << "Element ("<<qi<<"->"<<qf<<"): "<<Coefficient<<endl;
	  if (fabs(this->HoppingTerms[j].Im)<LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD) // real element
	    {
	      indexArray[positionR] = Index2;
	      tmpElementPos = RealInteractionCoefficients.InsertElement(Coefficient*this->HoppingTerms[j].Re);
	      if (tmpElementPos > USHRT_MAX )
		{
		  cout << "Error: too many different real matrix elements for fast storage"<<endl;
		  exit(1);
		}
	      coefficientIndexArray[positionR] = (unsigned short) tmpElementPos;
	      ++positionR;
	    }
	  else
	    {
	      indexArray[positionC] = Index2;
	      tmpElementPos = ComplexInteractionCoefficients.InsertElement
		(Coefficient*this->HoppingTerms[j]);
	      if (tmpElementPos > USHRT_MAX )
		{
		  cout << "Error: too many different complex matrix elements for fast storage"<<endl;
		  exit(1);
		}
	      coefficientIndexArray[positionC] = (unsigned short) tmpElementPos;
	      ++positionC;
	    }
	  //cout << "Hopping connecting :"<<Index2<<", "<<i<<": "<<Coefficient*this->HoppingTerms[j]<<endl;
	}
    }
}


// core part of the FastMultiplication method involving 2-body terms and diagonal elements
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientIndexArray = array of the indices of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnLattice* particles, int index, 
										 int* indexArray, unsigned short* coefficientIndexArray, int& positionR, int & positionC)
{
  int Index2;
  int tmpElementPos;
  double Coefficient;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  //       cout << "Interactions for component "<<i<<":  real "<< this->NbrRealInteractionPerComponent[position]<<" complex "
  // 	   << this->NbrComplexInteractionPerComponent[position] <<endl;    


  // four-fermion interactions:
  if (this->NbrQ12Indices == 0) // full storage
    { 	  
      for (int j = 0; j < NbrInteractionFactors; ++j) 
	{
	  int q1 = this->Q1Value[j];
	  int q2 = this->Q2Value[j];
	  int q3 = this->Q3Value[j];
	  int q4 = this->Q4Value[j];	       
	  Index2 = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
	  if (Index2 < Dim)
	    {
	      if (fabs(this->InteractionFactors[j].Im)<LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD) // real element
		{
		  indexArray[positionR] = Index2;
		  tmpElementPos = RealInteractionCoefficients.InsertElement(Coefficient*this->InteractionFactors[j].Re);
		  if (tmpElementPos > USHRT_MAX )
		    {
		      cout << "Error: too many different real matrix elements for fast storage"<<endl;
		      exit(1);
		    }
		  coefficientIndexArray[positionR] = (unsigned short) tmpElementPos;
		  ++positionR;
		}
	      else
		{
		  indexArray[positionC] = Index2;
		  tmpElementPos = ComplexInteractionCoefficients.InsertElement
		    (Coefficient*this->InteractionFactors[j]);
		  if (tmpElementPos > USHRT_MAX )
		    {
		      cout << "Error: too many different complex matrix elements for fast storage"<<endl;
		      exit(1);
		    }
		  coefficientIndexArray[positionC] = (unsigned short) tmpElementPos;
		  ++positionC;
		}
	      //cout << "4b - connecting :"<<Index2<<", "<<i<<": "<<Coefficient*this->InteractionFactors[j]<< " (q's=["<<q1<<","<<q2<<","<<q3<<","<<q4<<"])"<<endl;
	    }
	}
    }
  else // intelligent storage
    {
      double Coefficient2;
      int ProcessedNbrInteractionFactors = 0;
      int TmpNbrQ34Values;
      int* TmpQ3Values;
      int* TmpQ4Values;
      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
	{
	  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
	  if (Coefficient != 0.0)
	    {
	      TmpNbrQ34Values = this->NbrQ34Values[i12];
	      TmpQ3Values = this->Q3PerQ12[i12];
	      TmpQ4Values = this->Q4PerQ12[i12];
	      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
		{
		  Index2 = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
		  if (Index2 < Dim)
		    {
		      if (fabs(this->InteractionFactors[ProcessedNbrInteractionFactors].Im)<LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD) 
			{
			  indexArray[positionR] = Index2;
			  tmpElementPos = RealInteractionCoefficients.InsertElement
			    (Coefficient*Coefficient2*this->InteractionFactors[ProcessedNbrInteractionFactors].Re);
			  if (tmpElementPos > USHRT_MAX )
			    {
			      cout << "Error: too many different real matrix elements for fast storage"<<endl;
			      exit(1);
			    }
			  coefficientIndexArray[positionR] = (unsigned short) tmpElementPos;
			  ++positionR;
			}
		      else
			{
			  indexArray[positionC] = Index2;
			  tmpElementPos = ComplexInteractionCoefficients.InsertElement
			    (Coefficient*Coefficient2*this->InteractionFactors[ProcessedNbrInteractionFactors]);
			  if (tmpElementPos > USHRT_MAX )
			    {
			      cout << "Error: too many different complex matrix elements for fast storage"<<endl;
			      exit(1);
			    }
			  coefficientIndexArray[positionC] = (unsigned short) tmpElementPos;
			  ++positionC;
			}
		      //cout << "4b - connecting :"<<Index2<<", "<<i<<": "<<Coefficient<<"*"<<Coefficient2<<"*"<<this->InteractionFactors[ProcessedNbrInteractionFactors]<< " (q's=["<<this->Q1Value[i12]<<", "<<this->Q2Value[i12]<<", "<<TmpQ3Values[i34]<<", "<<TmpQ4Values[i34]<<"])"<<endl;
		    }
		  ++ProcessedNbrInteractionFactors;
		}
	    }
	  else
	    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
	}
    }
      
  // separated diagonal terms as these will be the general rule for contact interactions
  if (NbrDiagonalInteractionFactors>0)
    {
      // need additional symmetry factor of 1/2 in hermitian mode, as diagonal elements will not be treated separately if stored in memory!
      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors,
						 DiagonalInteractionFactors, DiagonalQValues);
      if (this->IsHermitian())
	Coefficient *= 0.5;
      if (fabs(Coefficient)>LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD)
	{
	  indexArray[positionR] = index;
	  tmpElementPos = RealInteractionCoefficients.InsertElement(Coefficient);
	  if (tmpElementPos > USHRT_MAX )
	    {
	      cout << "Error: too many different real matrix elements for fast storage"<<endl;
	      exit(1);
	    }
	  coefficientIndexArray[positionR] = (unsigned short) tmpElementPos;
	  ++positionR;
	  //cout << "diag - connecting :"<<i<<", "<<i<<": "<<Coefficient<<endl;
	}	   
    }
}


// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector& vSource, ComplexVector& vDestination)
{
	cout <<"Use new code "<<endl;
	if (this->NbrHoppingTerms>0)
	{
 int Dim = particles->GetHilbertSpaceDimension();  
 double Coefficient;
	// deal with kinetic energy terms first!      
	int qi;
	int qf;
	int Index;
	double TmpInteractionRe,TmpInteractionIm;
		  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
				TmpInteractionRe = this->HoppingTerms[j].Re;
	      TmpInteractionIm = this->HoppingTerms[j].Im;
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
		      vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteractionRe = this->HoppingTerms[ReducedNbrHoppingTerms].Re;
	  TmpInteractionIm = this->HoppingTerms[ReducedNbrHoppingTerms].Im;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
		  vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
		}
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    } 
	}
}

// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNOneBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector& vSource, ComplexVector& vDestination)
{
	cout <<"Use new code "<<endl;
	if (this->NbrHoppingTerms>0)
	{
 int Dim = particles->GetHilbertSpaceDimension();  
 double Coefficient;
	// deal with kinetic energy terms first!      
	int qi;
	int qf;
	int Index;
Complex TmpInteraction;
		  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
				TmpInteraction = this->HoppingTerms[j];
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
					vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms];
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		{
		  vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
		}
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    } 
	}
}

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
	cout <<"Use new code "<<endl;
	int Dim = particles->GetHilbertSpaceDimension(); 
	double Coefficient;
	int Index;
	
     if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
		 Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
					vDestination[Index] += Coefficient * this->InteractionFactors[j] * vSource[index];
		    }
		}
	}
      else // intelligent storage
	{
	  double Coefficient2, TmpRe, TmpIm;
		double TmpInteractionRe,TmpInteractionIm;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpRe = vSource[index].Re*Coefficient;
		      TmpIm = vSource[index].Im*Coefficient;
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpInteractionRe = this->InteractionFactors[ProcessedNbrInteractionFactors].Re;
			      TmpInteractionIm = this->InteractionFactors[ProcessedNbrInteractionFactors].Im;
			      vDestination.Re(Index) += Coefficient2 * (TmpRe*TmpInteractionRe-TmpIm*TmpInteractionIm);
			      vDestination.Im(Index) += Coefficient2 * (TmpRe*TmpInteractionIm+TmpIm*TmpInteractionRe);
			    }
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	    }

      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors, DiagonalInteractionFactors, DiagonalQValues);
	      vDestination[index] +=  Coefficient * vSource[index];
	}
}


// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNTwoBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
	cout <<"Use new code "<<endl;
	int Dim = particles->GetHilbertSpaceDimension(); 
	double Coefficient;
	int Index;
	
     if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
		 Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
					vDestination[index] += Conj(Coefficient * this->InteractionFactors[j]) * vSource[Index];
		    }
		}
	}
      else // intelligent storage
	{
	  double Coefficient2;
		Complex TmpSum;
	  int ProcessedNbrInteractionFactors= 0;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
						TmpSum += (Coefficient*Coefficient2)*Conj(this->InteractionFactors[ProcessedNbrInteractionFactors])*vSource[Index];	
			    }
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
		vDestination[index]+=TmpSum;
	    }

      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors, DiagonalInteractionFactors, DiagonalQValues);
				if (Coefficient != 0.0)
					vDestination[index] +=  Coefficient * vSource[index];
	}
}

// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
cout <<"Use new code "<<endl;
    if (this->NbrHoppingTerms>0)
	{
		int Dim = particles->GetHilbertSpaceDimension(); 
		int Index;
		double Coefficient;
		Complex TmpInteraction;
		Complex TmpCoefficient;
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j];
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms];
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      //cout << "element "<<qi<<"->"<<qf<<" on "<<i<<": "; 
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      //cout << "target: "<<Index<<endl;
	      if (Index < Dim)
		{
		  TmpCoefficient = Coefficient * TmpInteraction;
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += this->HamiltonianShift * vSources[l][i];
	    }
	}
}

// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNOneBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
cout <<"Use new code "<<endl;
    if (this->NbrHoppingTerms>0)
	{
		int Dim = particles->GetHilbertSpaceDimension(); 
		int Index;
		double Coefficient;
		Complex TmpInteraction;
		Complex TmpCoefficient;
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = Conj(this->HoppingTerms[j]);
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = Conj(this->HoppingTerms[ReducedNbrHoppingTerms]);
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      //cout << "element "<<qi<<"->"<<qf<<" on "<<i<<": "; 
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      //cout << "target: "<<Index<<endl;
	      if (Index < Dim)
		{
		  TmpCoefficient = Coefficient * TmpInteraction;
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += this->HamiltonianShift * vSources[l][i];
	    }
	}
}

// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector* vSources,  ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
	cout <<"Use new code "<<endl;
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
	Complex TmpInteraction;
	Complex TmpCoefficient;
	if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteraction = this->InteractionFactors[j];
		  Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][index];
		    }
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  Complex* TmpCoefficients = new Complex[nbrVectors];
	      ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      for (int l = 0; l < nbrVectors; ++l)
			tmpCoefficients[l] = vSources[l][index]*Coefficient;
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpCoefficient = this->InteractionFactors[ProcessedNbrInteractionFactors] * Coefficient2;
			      for (int l = 0; l < nbrVectors; ++l)
				vDestinations[l][Index] += TmpCoefficient * TmpCoefficients[l];
			    }
			  ++ProcessedNbrInteractionFactors;
			}
				}
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (this->NbrDiagonalInteractionFactors > 0)
	{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors , DiagonalInteractionFactors, DiagonalQValues);
	      for (int l = 0; l < nbrVectors; ++l)
					vDestinations[l][index] +=  Coefficient * vSources[l][index];
	}
      
}

// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNTwoBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector* vSources,  ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
	cout <<"Use new code "<<endl;
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
	Complex TmpInteraction;
	Complex TmpCoefficient;
	if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteraction = Conj(this->InteractionFactors[j]);
		  Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][index] += TmpCoefficient * vSources[l][Index];
		    }
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  int ProcessedNbrInteractionFactors= 0;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
		for (int l = 0; l < nbrVectors; ++l)
			tmpCoefficients[l] = 0.0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpCoefficient = Conj(this->InteractionFactors[ProcessedNbrInteractionFactors])*(Coefficient * Coefficient2);
			      for (int l = 0; l < nbrVectors; ++l)
				tmpCoefficients[l] += TmpCoefficient * vSources[l][Index];
			    }
			  ++ProcessedNbrInteractionFactors;
			}
				}
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
		for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][index] += tmpCoefficients[l];
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (this->NbrDiagonalInteractionFactors > 0)
	{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors , DiagonalInteractionFactors, DiagonalQValues);
	      for (int l = 0; l < nbrVectors; ++l)
					vDestinations[l][index] +=  Coefficient * vSources[l][index];
	}
      
}

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnLatticeHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{	
	cout <<"Use new code "<<endl;
	Complex TmpInteraction;
	int Index;
	int Dim = particles->GetHilbertSpaceDimension();
	double Coefficient;
	    // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteraction = this->InteractionFactors[j];
		  Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += Coefficient * TmpInteraction*vSource[index];
		      vDestination[index] += Coefficient * Conj(TmpInteraction)*vSource[Index];
		    }
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  Complex TmpSum;
	  Complex TmpCoefficient;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	      TmpSum=0.0;
	      ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      //cout << "Non-zero: a_"<<this->Q1Value[i12]<<" a_"<< this->Q2Value[i12]<<"| "<<i<<">  [following: "<<NbrQ34Values[i12]<<"]"<<endl;
		      TmpCoefficient = vSource[index];
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  //cout << "Testing a^+_"<<TmpQ3Values[i34]<<" a^+_"<< TmpQ4Values[i34]<<" ";
			  if (Index < Dim)
			    {
			      // cout << " -> "<<Index<<endl;
			      TmpInteraction = this->InteractionFactors[ProcessedNbrInteractionFactors];
			      Coefficient2*=Coefficient;
			      vDestination[Index] += Coefficient2 * TmpCoefficient * TmpInteraction;
			      //cout << "H["<<i<<", "<<Index<<"]="<<Coefficient2 * TmpInteraction<<endl;
			      TmpSum += Coefficient2*Conj(TmpInteraction)*vSource[Index];
			      //cout << "H["<<Index<<", "<<i<<"]="<<Coefficient2 * Conj(TmpInteraction)<<endl;
			    }
// 			  else
// 			    cout << " -> ZERO"<<endl;
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	      vDestination[index] += TmpSum;
	}	  
      // separated diagonal terms as these will be the general rule for contact interactions
      if (this->NbrDiagonalInteractionFactors > 0)
	{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors,
							DiagonalInteractionFactors, DiagonalQValues);
	      vDestination[index] +=  Coefficient * vSource[index];
	    }
}

inline void AbstractQHEOnLatticeHamiltonian::HermitianEvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector& vSource, ComplexVector& vDestination)
{
	cout <<"Use new code "<<endl;
	if (NbrHoppingTerms>0)
	{
			int Index;
	double Coefficient;
	Complex TmpInteraction;
	int Dim = particles->GetHilbertSpaceDimension();
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j];
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		      vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms];
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		{
		  vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		  vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
		}
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    }
	}
}

// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

inline void AbstractQHEOnLatticeHamiltonian::HermitianEvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
cout <<"Use new code "<<endl;
    if (this->NbrHoppingTerms>0)
	{
		int Dim = particles->GetHilbertSpaceDimension(); 
		int Index;
		double Coefficient;
		Complex TmpInteraction;
		Complex TmpCoefficient;
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j];
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
					      TmpCoefficient.Conjugate();
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms];
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      //cout << "element "<<qi<<"->"<<qf<<" on "<<i<<": "; 
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      //cout << "target: "<<Index<<endl;
	      if (Index < Dim)
		{
		  TmpCoefficient = Coefficient * TmpInteraction;
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
					  TmpCoefficient.Conjugate();
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += this->HamiltonianShift * vSources[l][i];
	    }
	}
}

// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void AbstractQHEOnLatticeHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, ComplexVector* vSources,  ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
	cout <<"Use new code "<<endl;
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
	Complex TmpInteraction;
	Complex TmpCoefficient;
	if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteraction = this->InteractionFactors[j];
		  Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][index];
					      TmpCoefficient.Conjugate();
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][index] += TmpCoefficient * vSources[l][Index];
		    }
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  int ProcessedNbrInteractionFactors = 0;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	      ProcessedNbrInteractionFactors = 0;
				for (int l = 0; l < nbrVectors; ++l)
		tmpCoefficients[l] = 0.0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpCoefficient = this->InteractionFactors[ProcessedNbrInteractionFactors] * Coefficient2 *Coefficient;
			      for (int l = 0; l < nbrVectors; ++l)
				vDestinations[l][Index] += TmpCoefficient * vSources[l][index];
				TmpCoefficient.Conjugate();
			      for (int l = 0; l < nbrVectors; ++l)
				tmpCoefficients[l] += TmpCoefficient * vSources[l][Index];
			    }
			  ++ProcessedNbrInteractionFactors;
			}
				}
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
		for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][index] += tmpCoefficients[l];
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (this->NbrDiagonalInteractionFactors > 0)
	{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors , DiagonalInteractionFactors, DiagonalQValues);
	      for (int l = 0; l < nbrVectors; ++l)
					vDestinations[l][index] +=  Coefficient * vSources[l][index];
	}
      
}


// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector& vSource, RealVector& vDestination)
{
	cout <<"Use new code "<<endl;
	int Dim = particles->GetHilbertSpaceDimension(); 
	double Coefficient;
	int Index;
	
     if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
		 Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
					vDestination[Index] += Coefficient * this->InteractionFactors[j].Re * vSource[index];
		    }
		}
	}
      else // intelligent storage
	{
	  double Coefficient2, TmpRe;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpRe = vSource[index]*Coefficient;
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
 vDestination[Index] += Coefficient2 * TmpRe * this->InteractionFactors[ProcessedNbrInteractionFactors].Re;
			    }
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	    }

      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors, DiagonalInteractionFactors, DiagonalQValues);
	      vDestination[index] +=  Coefficient * vSource[index];
	}
}

// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector& vSource, RealVector& vDestination)
{
	cout <<"Use new code "<<endl;
	if (this->NbrHoppingTerms>0)
	{
 int Dim = particles->GetHilbertSpaceDimension();  
 double Coefficient;
	// deal with kinetic energy terms first!      
	int qi;
	int qf;
	int Index;
	double TmpInteraction;
		  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
				TmpInteraction = this->HoppingTerms[j].Re;
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
					vDestination[Index] += Coefficient * TmpInteraction*vSource[i];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
		TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms].Re;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		{
			vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		}
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    } 
	}
}

// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
cout <<"Use new code "<<endl;
    if (this->NbrHoppingTerms>0)
	{
		int Dim = particles->GetHilbertSpaceDimension(); 
		int Index;
		double Coefficient;
  double TmpInteraction;
  double TmpCoefficient;
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j].Re;
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms].Re;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      //cout << "element "<<qi<<"->"<<qf<<" on "<<i<<": "; 
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      //cout << "target: "<<Index<<endl;
	      if (Index < Dim)
		{
		  TmpCoefficient = Coefficient * TmpInteraction;
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += this->HamiltonianShift * vSources[l][i];
	    }
	}
}

// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector* vSources,  RealVector* vDestinations, int nbrVectors, double* tmpCoefficients)
{
	cout <<"Use new code "<<endl;
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
	double TmpCoefficient;
	if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
		  Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * this->InteractionFactors[j].Re;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][index];
		    }
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	      ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      for (int l = 0; l < nbrVectors; ++l)
			tmpCoefficients[l] = vSources[l][index]*Coefficient;
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpCoefficient = this->InteractionFactors[ProcessedNbrInteractionFactors].Re  * Coefficient2;
			      for (int l = 0; l < nbrVectors; ++l)
				vDestinations[l][Index] += TmpCoefficient * tmpCoefficients[l];
			    }
			  ++ProcessedNbrInteractionFactors;
			}
				}
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (this->NbrDiagonalInteractionFactors > 0)
	{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors , DiagonalInteractionFactors, DiagonalQValues);
	      for (int l = 0; l < nbrVectors; ++l)
					vDestinations[l][index] +=  Coefficient * vSources[l][index];
	}
      
}

// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNOneBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector& vSource, RealVector& vDestination)
{
	cout <<"Use new code "<<endl;
	if (this->NbrHoppingTerms>0)
	{
 int Dim = particles->GetHilbertSpaceDimension();  
 double Coefficient;
	// deal with kinetic energy terms first!      
	int qi;
	int qf;
	int Index;
double TmpInteraction;
		  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
				TmpInteraction = this->HoppingTerms[j].Re;
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
					vDestination[i] += Coefficient * TmpInteraction *vSource[Index];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms].Re;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		{
		  vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
		}
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    } 
	}
}

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNTwoBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector& vSource, RealVector& vDestination)
{
	cout <<"Use new code "<<endl;
	int Dim = particles->GetHilbertSpaceDimension(); 
	double Coefficient;
	int Index;
	
     if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
		 Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
					vDestination[index] += Coefficient * this->InteractionFactors[j].Re * vSource[Index];
		    }
		}
	}
      else // intelligent storage
	{
	  double Coefficient2;
		double TmpSum = 0.0;
	  int ProcessedNbrInteractionFactors= 0;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
						TmpSum += (Coefficient*Coefficient2)*this->InteractionFactors[ProcessedNbrInteractionFactors].Re*vSource[Index];
			    }
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
		vDestination[index] += TmpSum;
	    }

      // separated diagonal terms as these will be the general rule for contact interactions
      if (this->NbrDiagonalInteractionFactors > 0)
	{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors, DiagonalInteractionFactors, DiagonalQValues);
				if (Coefficient != 0.0)
					vDestination[index] +=  Coefficient * vSource[index];
	}
}
// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNOneBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
cout <<"Use new code "<<endl;
    if (this->NbrHoppingTerms>0)
	{
		int Dim = particles->GetHilbertSpaceDimension(); 
		int Index;
		double Coefficient;
		double TmpInteraction;
		double TmpCoefficient;
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j].Re;
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms].Re;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		{
		  TmpCoefficient = Coefficient * TmpInteraction;
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += this->HamiltonianShift * vSources[l][i];
	    }
	}
}

// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNTwoBodyConjugateAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector* vSources,  RealVector* vDestinations, int nbrVectors, double* tmpCoefficients)
{
	cout <<"Use new code "<<endl;
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
	double TmpCoefficient;
	if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
		  Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * this->InteractionFactors[j].Re;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][index] += TmpCoefficient * vSources[l][Index];
		    }
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  int ProcessedNbrInteractionFactors= 0;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
		for (int l = 0; l < nbrVectors; ++l)
		tmpCoefficients[l] = 0.0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpCoefficient = this->InteractionFactors[ProcessedNbrInteractionFactors].Re * Coefficient * Coefficient2;
			      for (int l = 0; l < nbrVectors; ++l)
				tmpCoefficients[l]+= TmpCoefficient * vSources[l][Index];
			    }
			  ++ProcessedNbrInteractionFactors;
			}
				}
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
		for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][index] += tmpCoefficients[l];
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (this->NbrDiagonalInteractionFactors > 0)
	{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors , DiagonalInteractionFactors, DiagonalQValues);
	      for (int l = 0; l < nbrVectors; ++l)
					vDestinations[l][index] +=  Coefficient * vSources[l][index];
	}
      
}

inline void AbstractQHEOnLatticeHamiltonian::HermitianEvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector& vSource, RealVector& vDestination)
{
	cout <<"Use new code "<<endl;
	if (NbrHoppingTerms>0)
	{
			int Index;
	double Coefficient;
	double TmpInteraction;
	int Dim = particles->GetHilbertSpaceDimension();
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j].Re;
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		      vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms].Re;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		{
		  vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		  vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
		}
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    }
	}
}

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void AbstractQHEOnLatticeHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector& vSource, RealVector& vDestination)
{	
	cout <<"Use new code "<<endl;
	double TmpInteraction;
	int Index;
	int Dim = particles->GetHilbertSpaceDimension();
	double Coefficient;
	    // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteraction = this->InteractionFactors[j].Re;
		  Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += Coefficient * TmpInteraction * vSource[index];
		      vDestination[index] += Coefficient * TmpInteraction * vSource[Index];
		    }
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  double TmpSum = 0.0;
	  double TmpCoefficient;
	  int ProcessedNbrInteractionFactors= 0;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;

	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      //cout << "Non-zero: a_"<<this->Q1Value[i12]<<" a_"<< this->Q2Value[i12]<<"| "<<i<<">  [following: "<<NbrQ34Values[i12]<<"]"<<endl;
		      TmpCoefficient = vSource[index];
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  //cout << "Testing a^+_"<<TmpQ3Values[i34]<<" a^+_"<< TmpQ4Values[i34]<<" ";
			  if (Index < Dim)
			    {
			      // cout << " -> "<<Index<<endl;
			      TmpInteraction = this->InteractionFactors[ProcessedNbrInteractionFactors].Re;
			      Coefficient2*=Coefficient;
			      vDestination[Index] += Coefficient2 * TmpCoefficient * this->InteractionFactors[ProcessedNbrInteractionFactors].Re;
			      //cout << "H["<<i<<", "<<Index<<"]="<<Coefficient2 * TmpInteraction<<endl;
			      TmpSum += Coefficient2 * TmpInteraction * vSource[Index];
			      //cout << "H["<<Index<<", "<<i<<"]="<<Coefficient2 * Conj(TmpInteraction)<<endl;
			    }
// 			  else
// 			    cout << " -> ZERO"<<endl;
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	      vDestination[index] += TmpSum;
	}	  
      // separated diagonal terms as these will be the general rule for contact interactions
      if (this->NbrDiagonalInteractionFactors > 0)
			{
	      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors, DiagonalInteractionFactors, DiagonalQValues);
	      vDestination[index] +=  Coefficient * vSource[index];
	    }
}

// core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together

inline void AbstractQHEOnLatticeHamiltonian::HermitianEvaluateMNOneBodyAddMultiplyComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, int step, RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
cout <<"Use new code "<<endl;
    if (this->NbrHoppingTerms>0)
	{
		int Dim = particles->GetHilbertSpaceDimension(); 
		int Index;
		double Coefficient;
  double TmpInteraction;
  double TmpCoefficient;
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j].Re;
	      for (int i = firstComponent; i < lastComponent; i += step)
		{
		  Index = particles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
					{
			vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
			vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
					}
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms].Re;
	  for (int i = firstComponent; i < lastComponent; i += step)
	    {
	      //cout << "element "<<qi<<"->"<<qf<<" on "<<i<<": "; 
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      //cout << "target: "<<Index<<endl;
	      if (Index < Dim)
		{
		  TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
					{
			vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
			vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
					}
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += this->HamiltonianShift * vSources[l][i];
	    }
	}
}

// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

inline void AbstractQHEOnLatticeHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnLattice* particles, int index, RealVector* vSources,  RealVector* vDestinations, int nbrVectors, double* tmpCoefficients)
{
	cout <<"Use new code "<<endl;
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
	double TmpCoefficient;
	if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
		  Index = particles->AdAdAA(index, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * this->InteractionFactors[j].Re;
		      for (int l = 0; l < nbrVectors; ++l)
					{
			vDestinations[l][Index] += TmpCoefficient * vSources[l][index];
			vDestinations[l][index] += TmpCoefficient * vSources[l][Index];
					}
		    }
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  int ProcessedNbrInteractionFactors = 0;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
				for (int l = 0; l < nbrVectors; ++l)
		tmpCoefficients[l] = 0.0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = particles->AA(index, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpCoefficient = this->InteractionFactors[ProcessedNbrInteractionFactors].Re  * Coefficient2 *Coefficient;
			      for (int l = 0; l < nbrVectors; ++l)
						{
							vDestinations[l][Index] += TmpCoefficient * vSources[l][index];
							tmpCoefficients[l] += TmpCoefficient * vSources[l][Index];
						}
					}
			  ++ProcessedNbrInteractionFactors;
			}
				}
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
		for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][index] += tmpCoefficients[l];
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (this->NbrDiagonalInteractionFactors > 0)
	{
	    Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors , DiagonalInteractionFactors, DiagonalQValues);
	      for (int l = 0; l < nbrVectors; ++l)
					vDestinations[l][index] +=  Coefficient * vSources[l][index];
	}
      
}

// core part of the PartialFastMultiplicationMemory method involving two-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int Dim = particles->GetHilbertSpaceDimension();
if (this->NbrQ12Indices == 0) // full storage
    {
      for (int j = 0; j < NbrInteractionFactors; ++j) 
	{
	  int q1 = this->Q1Value[j];
	  int q2 = this->Q2Value[j];
	  int q3 = this->Q3Value[j];
	  int q4 = this->Q4Value[j];
	  if (fabs(this->InteractionFactors[j].Im) < LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD)
	    {
	      for (int i = firstComponent; i < lastComponent; ++i)
		{
		  Index = particles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];
		    }
		}
	    }
	  else
	    {
	      for (int i = firstComponent; i < lastComponent; ++i)
		{
		  Index = particles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      ++memory;
		      ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];
		    }
		}
	    }	 	  
	}
    }
  else // intelligent storage
    {
      double Coefficient2;
      int ProcessedNbrInteractionFactors;
      int TmpNbrQ34Values;
      int* TmpQ3Values;
      int* TmpQ4Values;
      for (int i = firstComponent; i < lastComponent; ++i)
	{	  
	  ProcessedNbrInteractionFactors = 0;
	  for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
	    {
	      Coefficient = particles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
	      if (Coefficient != 0.0)
		{
		  TmpNbrQ34Values = this->NbrQ34Values[i12];
		  TmpQ3Values = this->Q3PerQ12[i12];
		  TmpQ4Values = this->Q4PerQ12[i12];
		  for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
		    {
		      Index = particles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
		      if (Index < Dim)
			{
			  ++memory;
			  if (fabs(this->InteractionFactors[ProcessedNbrInteractionFactors].Im) < LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD)
			    ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];
			  else
			    ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];
			  // cout << "4b - connecting :"<<Index<<", "<<i<<": "<<Coefficient<<"*"<<Coefficient2<<"*"<<this->InteractionFactors[ProcessedNbrInteractionFactors]<< " (q's=["<<this->Q1Value[i12]<<", "<<this->Q2Value[i12]<<", "<<TmpQ3Values[i34]<<", "<<TmpQ4Values[i34]<<"])"<<endl;
			}
		      ++ProcessedNbrInteractionFactors;
		    }
		}
	      else
		ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
	    }
	}
    }	  
  
  // separated diagonal terms as these will be the general rule for contact interactions
  if (this->NbrDiagonalInteractionFactors > 0)
    {	  
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  Coefficient = particles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors, DiagonalInteractionFactors, DiagonalQValues);
	      ++memory;
	      ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];
	}
    }
}

// core part of the PartialFastMultiplicationMemory method involving two-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void AbstractQHEOnLatticeHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnLattice* particles, int firstComponent, int lastComponent, long& memory)
{
	int Index;
  double Coefficient;
	int Dim = particles->GetHilbertSpaceDimension();
  // deal with kinetic energy terms first!      
  int qi;
  int qf;
  for (int j = 0; j < this->NbrHoppingTerms; ++j)
    {
      qi = this->KineticQi[j];
      qf = this->KineticQf[j];
      if (fabs(this->HoppingTerms[j].Im)<LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		{
			++memory;			
		  ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];		
		}
	    }
	}
      else
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      Index = particles->AdA(i, qf, qi, Coefficient);
	      
				if (Index < Dim)
		{
		  ++memory;
		  ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];
		}
	    }
	}
    }
}

#endif
