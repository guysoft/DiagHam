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
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
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
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
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
  void EvaluateFastMultiplicationComponent(ParticleOnLattice* particles, int index, 
					   int* indexArray, unsigned short* coefficientIndexArray, long& position);

  
};


// core part of the FastMultiplication method involving 1- and 2-body terms and diagonal elements
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientIndexArray = array of the indices of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void AbstractQHEOnLatticeHamiltonian::EvaluateFastMultiplicationComponent(ParticleOnLattice* particles, int index, 
										 int* indexArray, unsigned short* coefficientIndexArray, long& position)
{
  int qi, qf;
  int Index2;
  int tmpElementPos;
  double Coefficient;
  int PosR, PosC;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  int TmpInteractionPerComponent = this->NbrRealInteractionPerComponent[position]+this->NbrComplexInteractionPerComponent[position];
  //       cout << "Interactions for component "<<i<<":  real "<< this->NbrRealInteractionPerComponent[position]<<" complex "
  // 	   << this->NbrComplexInteractionPerComponent[position] <<endl;
  this->InteractionPerComponentIndex[position] = new int [TmpInteractionPerComponent];
  this->InteractionPerComponentCoefficientIndex[position] = new unsigned short [TmpInteractionPerComponent];      
  indexArray = this->InteractionPerComponentIndex[position];
  coefficientIndexArray = this->InteractionPerComponentCoefficientIndex[position];
  PosR = 0;  // counter for position of real matrix elements
  PosC = this->NbrRealInteractionPerComponent[position];  // counter for position of complex matrix elements

  // deal with kinetic energy terms first!             
  for (int j = 0; j < NbrHoppingTerms; ++j)
    {
      qi = this->KineticQi[j];
      qf = this->KineticQf[j];
      // considering: this->HoppingTerms[j]
      Index2 = particles->AdA(index, qf, qi, Coefficient);	  
      if (Index2 < Dim)
	{
	  //cout << "Element ("<<qi<<"->"<<qf<<"): "<<Coefficient<<endl;
	  if (fabs(this->HoppingTerms[j].Im)<LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD) // real element
	    {
	      indexArray[PosR] = Index2;
	      tmpElementPos = RealInteractionCoefficients.InsertElement
		(Coefficient*this->HoppingTerms[j].Re);
	      if (tmpElementPos > USHRT_MAX )
		{
		  cout << "Error: too many different real matrix elements for fast storage"<<endl;
		  exit(1);
		}
	      coefficientIndexArray[PosR] = (unsigned short) tmpElementPos;
	      ++PosR;
	    }
	  else
	    {
	      indexArray[PosC] = Index2;
	      tmpElementPos = ComplexInteractionCoefficients.InsertElement
		(Coefficient*this->HoppingTerms[j]);
	      if (tmpElementPos > USHRT_MAX )
		{
		  cout << "Error: too many different complex matrix elements for fast storage"<<endl;
		  exit(1);
		}
	      coefficientIndexArray[PosC] = (unsigned short) tmpElementPos;
	      ++PosC;
	    }
	  //cout << "Hopping connecting :"<<Index2<<", "<<i<<": "<<Coefficient*this->HoppingTerms[j]<<endl;
	}
    }

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
		  indexArray[PosR] = Index2;
		  tmpElementPos = RealInteractionCoefficients.InsertElement
		    (Coefficient*this->InteractionFactors[j].Re);
		  if (tmpElementPos > USHRT_MAX )
		    {
		      cout << "Error: too many different real matrix elements for fast storage"<<endl;
		      exit(1);
		    }
		  coefficientIndexArray[PosR] = (unsigned short) tmpElementPos;
		  ++PosR;
		}
	      else
		{
		  indexArray[PosC] = Index2;
		  tmpElementPos = ComplexInteractionCoefficients.InsertElement
		    (Coefficient*this->InteractionFactors[j]);
		  if (tmpElementPos > USHRT_MAX )
		    {
		      cout << "Error: too many different complex matrix elements for fast storage"<<endl;
		      exit(1);
		    }
		  coefficientIndexArray[PosC] = (unsigned short) tmpElementPos;
		  ++PosC;
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
			  indexArray[PosR] = Index2;
			  tmpElementPos = RealInteractionCoefficients.InsertElement
			    (Coefficient*Coefficient2*this->InteractionFactors[ProcessedNbrInteractionFactors].Re);
			  if (tmpElementPos > USHRT_MAX )
			    {
			      cout << "Error: too many different real matrix elements for fast storage"<<endl;
			      exit(1);
			    }
			  coefficientIndexArray[PosR] = (unsigned short) tmpElementPos;
			  ++PosR;
			  if (PosR>this->NbrRealInteractionPerComponent[position])
			    cout << "Problem with count of real matrix elements"<<endl;
			}
		      else
			{
			  indexArray[PosC] = Index2;
			  tmpElementPos = ComplexInteractionCoefficients.InsertElement
			    (Coefficient*Coefficient2*this->InteractionFactors[ProcessedNbrInteractionFactors]);
			  if (tmpElementPos > USHRT_MAX )
			    {
			      cout << "Error: too many different complex matrix elements for fast storage"<<endl;
			      exit(1);
			    }
			  coefficientIndexArray[PosC] = (unsigned short) tmpElementPos;
			  ++PosC;
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
      Coefficient = particles->AdAdAADiagonal(index, NbrDiagonalInteractionFactors,
						 DiagonalInteractionFactors, DiagonalQValues);
      if (fabs(Coefficient)>LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD)
	{
	  indexArray[PosR] = index;
	  tmpElementPos = RealInteractionCoefficients.InsertElement(Coefficient);
	  if (tmpElementPos > USHRT_MAX )
	    {
	      cout << "Error: too many different real matrix elements for fast storage"<<endl;
	      exit(1);
	    }
	  coefficientIndexArray[PosR] = (unsigned short) tmpElementPos;
	  ++PosR;
	  //cout << "diag - connecting :"<<i<<", "<<i<<": "<<Coefficient<<endl;
	}	   
    }
  ++position;
}

#endif
