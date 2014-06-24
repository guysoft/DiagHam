////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//                   class of Hubbard hamiltonian associated                  //
//                to particles on a lattice with spin                         //
//                                                                            //
//                        last modification : 19/06/2014                      //
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


#ifndef PARTICLEONLATTICEWITHSPINKITAEVHEISENBERGHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINKITAEVHEISENBERGHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian : public AbstractQHEHamiltonian
{

 protected:
  
  // Hilbert space associated to the system
  ParticleOnSphereWithSpin* Particles;

  // number of particles
  int NbrParticles;

  // number of sites
  int NbrSite;
  //number of orbitals
  int LzMax;
  
  //array that contains the information about the lattice geometry
  //first index = index of the site under consideration
  //second index 0 = index of the nearest neighbor linked with x link (= NbrSite if no x link)
  //second index 1 = index of the nearest neighbor linked with y link (= NbrSite if no y link)
  //second index 2 = index of the nearest neighbor linked with z link (= NbrSite if no z link)
  int** MapNearestNeighborBonds;
  
  // Hubbard potential strength
  double UPotential;
  
  // multiplicative factor in front of the isotropic nearest neighbor hopping term
  double KineticFactorIsotropic;
  // multiplicative factor in front of the anisotropic nearest neighbor hopping term
  double KineticFactorAnisotropic;
  
  //strength of the isotropic nearest neighbor spin interaction
  double J1Factor;
  //strength of the anisotropic nearest neighbor spin interaction
  double J2Factor;

  // shift to apply to go from precalculation index to the corresponding index in the HilbertSpace
  int PrecalculationShift;

  // shift to apply to the Hamiltonian diagonal elements
  double HamiltonianShift;

  // amount of memory (in bytes) that can be used to store precalculated matrix elements
  long Memory;

  // arrays containing all interaction factors, the first index corresponds to index sum for the creation (or annhilation) operators
  // the second index is a linearized index (m1,m2) + (n1,n2) * (nbr element in current index sum) (m for creation operators, n for annhilation operators)
  // array containing all interaction factors for Ad_{ui}Ad{uj}A_{ui}A_{uj}
  double** InteractionFactorsupup;
  // array containing all interaction factors for Ad_{di}Ad{dj}A_{di}A_{dj}
  double** InteractionFactorsdowndown;
  // array containing all interaction factors for Ad_{ui}Ad{udj}A_{ui}A_{dj}
  double** InteractionFactorsupdown;
  // array containing all interaction factors for Ad_{ui}Ad{uj}A_{di}A_{dj}
  double** InteractionFactorsupupdowndown;
  // array containing all interaction factors for Ad_{ui}Ad{dj}A_{di}A_{uj}
  double** InteractionFactorsupdowndownup;

  // array that contains all one-body interaction factors for particles with spin up
  double** OneBodyInteractionFactorsupup;
  // array that contains all one-body interaction factors for particles with spin down
  double** OneBodyInteractionFactorsdowndown;
  // array that contains all one-body interaction factors for tunnelling terms for particles with different spin
  Complex** OneBodyInteractionFactorsupdown;
  
  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // step between each precalculated index
  int FastMultiplicationStep;

  // stored interactions per component
  int *NbrInteractionPerComponent;

  // indices of matrix elements per component
  int **InteractionPerComponentIndex;
  // coefficients of matrix elements per component
  Complex** InteractionPerComponentCoefficient;
  // translations of matrix elements per component
  int **InteractionPerComponentNbrTranslation;
  
  //array containing all the cosinus that are needed when computing matrix elements
//   double* CosinusTable;
  //array containing all the sinus that are needed when computing matrix elements
//   double* SinusTable;

  // flag for implementation of hermitian symmetry
  bool HermitianSymmetryFlag;

  // pointer to an optional S^2 operator in the Hamiltonian 
//   ParticleOnSphereWithSpinS2Hamiltonian* S2Hamiltonian;

  
 public:

  // default constructor
  //
  ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSites = number of sites
  // kineticFactorIntra = multiplicative factor in front of the intraspin kinetic term
  // uPotential = Hubbard potential strength
  //j1Factor = strength of the isotropic nearest neighbor spin interaction
  //j2Factor = strength of the anisotropic nearest neighbor spin interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSite, double kineticFactorIsotropic, double kineticFactorAnisotropic, double uPotential, double j1Factor, double j2Factor, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  virtual ~ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian();
  
  // ask if Hamiltonian implements hermitian symmetry operations
  //
  virtual bool IsHermitian();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

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

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
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
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
//   virtual ComplexVector& HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 						      int firstComponent, int nbrComponent);
 
  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
//   virtual ComplexVector* HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 							      int firstComponent, int nbrComponent);

  //get the index of the nearest neighbor linked with bond x
  //
  // i = index of the site under consideration
  // return value = index of the site linked to i with x bond
  int GetIndexNearestNeighborXBond(int i);
  
  //get the index of the nearest neighbor linked with bond y
  //
  // i = index of the site under consideration
  // return value = index of the site linked to i with y bond
  int GetIndexNearestNeighborYBond(int i);
  
  //get the index of the nearest neighbor linked with bond z
  //
  // i = index of the site under consideration
  // return value = index of the site linked to i with z bond
  int GetIndexNearestNeighborZBond(int i);
  
  
 protected:
 
  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
  virtual void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
						     ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
//   virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the two-body interaction for a set of vectors
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // tmpCoefficients = a temporary array whose size is nbrVectors
//   virtual void HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 							      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients);

  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the FastMultiplication method involving the one-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  virtual void EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
							    int* indexArray, Complex* coefficientArray, long& position);

  // core part of the AddMultiply method involving the one-body interaction, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
						     int step, ComplexVector& vSource, ComplexVector& vDestination);

  // core part of the AddMultiply method involving the one-body interaction for a set of vectors, including loop on vector components
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = first index vector to act on
  // lastComponent = last index vector to act on (not included)
  // step = step to go from one component to the other one
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  virtual void EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
						     int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors);

  // core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);

  // core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  virtual void EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);


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
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
//   virtual ComplexVector* HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 									 int firstComponent, int nbrComponent);

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // test the amount of memory needed for fast multiplication algorithm
  //
  // allowedMemory = amount of memory that cam be allocated for fast multiplication
  // return value = amount of memory needed
  virtual long FastMultiplicationMemory(long allowedMemory);

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
  // nbrComponent  = index of the last component that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int nbrComponent);
  
  
  //fill up the positions of the three types of nearest neighbor bonds
  //
  void PlotMapNearestNeighborBonds();

};

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  Complex Coefficient4;
  int* TmpIndices;
  double* TmpInteractionFactor;
  int Index;
  for (int j = 0; j < this->NbrSite; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
	int j2 = this->MapNearestNeighborBonds[j][k];
	if (j2 < this->NbrSite)
	{
	  Coefficient3 = particles->AuAu(index, j, j2);
	  if (Coefficient3 != 0.0)
	   {
	    Coefficient4 = vSource[index];
	    Coefficient4 *= Coefficient3;
	    
	    TmpInteractionFactor = &(this->InteractionFactorsupup[j][k]);
	    Index = particles->AduAdu(j, j2, Coefficient);
	    if (Index < Dim)
	      vDestination[index] += Coefficient * Coefficient4 * (*TmpInteractionFactor);
	    ++TmpInteractionFactor;
	  
	    TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][k]);
	    Index = particles->AddAdd(j, j2, Coefficient);
	    if (Index < Dim)
	      vDestination[index] += Coefficient * Coefficient4 * (*TmpInteractionFactor);
	    ++TmpInteractionFactor;
	   }

	  Coefficient3 = particles->AdAd(index, j, j2);
	  if (Coefficient3 != 0.0)
	   {
	    Coefficient4 = vSource[index];
	    Coefficient4 *= Coefficient3;
	      
	    TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][k]);
	    Index = particles->AddAdd(j, j2, Coefficient);
	    if (Index < Dim)
	      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
	    ++TmpInteractionFactor;
	      
	    TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][k]);
	    Index = particles->AduAdu(j, j2, Coefficient);
	    if (Index < Dim)
	      vDestination[index] += Coefficient * Coefficient4 * (*TmpInteractionFactor);
	    ++TmpInteractionFactor;
	    }
	    
	    Coefficient3 = particles->AuAd(index, j, j2);
	    if (Coefficient3 != 0.0)
	    {
	      Coefficient4 = vSource[index];
	      Coefficient4 *= Coefficient3;
	      
	      TmpInteractionFactor = &(this->InteractionFactorsupdowndownup[j][k]);
	      Index = particles->AduAdd(j2, j, Coefficient);
	      if (Index < Dim)
		vDestination[index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
	      ++TmpInteractionFactor;
	      
	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][k]);
	      Index = particles->AduAdd(j, j2, Coefficient);
	      if (Index < Dim)
		vDestination[index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
	      ++TmpInteractionFactor;
	    }
	    
// 	    Coefficient3 = particles->AdAu(index, j, j2);
// 	    if (Coefficient3 != 0.0)
// 	    {
// 	      Coefficient4 = vSource[index];
// 	      Coefficient4 *= Coefficient3;
// 	      
// 	      TmpInteractionFactor = &(this->InteractionFactorsupdowndownup[j][k]);
// 	      Index = particles->AduAdd(j, j2, Coefficient);
// 	      if (Index < Dim)
// 		vDestination[index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
// 	      ++TmpInteractionFactor;
// 	      
// 	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][k]);
// 	      Index = particles->AddAdu(j, j2, Coefficient);
// 	      if (Index < Dim)
// 		vDestination[index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
// 	      ++TmpInteractionFactor;
// 	    }
	}
      }
    if (this->UPotential != 0)
      {
	Coefficient3 = particles->AuAd(index, j, j);
	if (Coefficient3 != 0)
	{
	  Coefficient4 = vSource[index];
	  Coefficient4 *= Coefficient3;
	  TmpInteractionFactor = &(this->UPotential);
	  Index = particles->AduAdd(j, j, Coefficient);
	  if (Index < Dim)
	    vDestination[index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
	  ++TmpInteractionFactor;
	}
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

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
												      ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;
  int Index;
  
  int* TmpIndices;
  double* TmpInteractionFactor;
  for (int j = 0; j < this->NbrSite; ++j)
    {
      for (int k = 0; k < 3; ++k)
	{
	  int j2 = this->MapNearestNeighborBonds[j][k];
	  if (j2 < this->NbrSite)
	  {
	    Coefficient3 = particles->AuAu(index, j, j2);
	    if (Coefficient3 != 0.0)
	      {
		for (int p = 0; p < nbrVectors; ++p)
		  tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	      
		TmpInteractionFactor = &(this->InteractionFactorsupup[j][k]);
		Index = particles->AduAdu(j, j2, Coefficient);
		if (Index < Dim)
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		++TmpInteractionFactor;
		
		TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][k]);
		Index = particles->AddAdd(j, j2, Coefficient);
		if (Index < Dim)
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		++TmpInteractionFactor;
		}
	    
	    Coefficient3 = particles->AdAd(index, j, j2);
	    if (Coefficient3 != 0.0)
	      {
		for (int p = 0; p < nbrVectors; ++p)
		  tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		
		TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][k]);
		Index = particles->AddAdd(j, j2, Coefficient);
		if (Index < Dim)
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		++TmpInteractionFactor;
		
		TmpInteractionFactor = &(this->InteractionFactorsupupdowndown[j][k]);
		Index = particles->AduAdu(j, j2, Coefficient);
		if (Index < Dim)
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		++TmpInteractionFactor;
	      }
	      
	    Coefficient3 = particles->AuAd(index, j, j2);
	    if (Coefficient3 != 0.0)
	      {
		for (int p = 0; p < nbrVectors; ++p)
		  tmpCoefficients[p] = Coefficient3 * vSources[p][index];
		
		TmpInteractionFactor = &(this->InteractionFactorsupdown[j][k]);
		Index = particles->AduAdd(j, j2, Coefficient);
		if (Index < Dim)
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		++TmpInteractionFactor;
		
		TmpInteractionFactor = &(this->InteractionFactorsupdowndownup[j][k]);
		Index = particles->AduAdd(j2, j, Coefficient);
		if (Index < Dim)
		  for (int p = 0; p < nbrVectors; ++p)
		    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
		++TmpInteractionFactor;
	      }
	     
// 	     Coefficient3 = particles->AdAu(index, j, j2);
// 	     if (Coefficient3 != 0.0)
// 	      {
// 		for (int p = 0; p < nbrVectors; ++p)
// 		  tmpCoefficients[p] = Coefficient3 * vSources[p][index];
// 		
// 		TmpInteractionFactor = &(this->InteractionFactorsupdown[j][k]);
// 		Index = particles->AddAdu(j, j2, Coefficient);
// 		if (Index < Dim)
// 		  for (int p = 0; p < nbrVectors; ++p)
// 		    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
// 		++TmpInteractionFactor;
// 		
// 		TmpInteractionFactor = &(this->InteractionFactorsupdowndownup[j][k]);
// 		Index = particles->AduAdd(j, j2, Coefficient);
// 		if (Index < Dim)
// 		  for (int p = 0; p < nbrVectors; ++p)
// 		    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
// 		++TmpInteractionFactor;
// 	      }
	}
    }
    if (this->UPotential != 0)
    {
      Coefficient3 = particles->AuAd(index, j, j);
      if (Coefficient3 != 0)
	{
	  for (int p = 0; p < nbrVectors; ++p)
	    tmpCoefficients[p] = Coefficient3 * vSources[p][index];
	  
	  TmpInteractionFactor = &(this->UPotential);
	  Index = particles->AduAdd(j, j, Coefficient);
	  if (Index < Dim)
	    for (int p = 0; p < nbrVectors; ++p)
	      vDestinations[p][index] += (Coefficient * (*TmpInteractionFactor)) * tmpCoefficients[p];
	  ++TmpInteractionFactor;
	}
    }
  }
}

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

// inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector& vSource, ComplexVector& vDestination)
// {
//   //int Dim = particles->GetHilbertSpaceDimension();
//   double Coefficient;
//   double Coefficient3;
//   Complex Coefficient4;
//   int* TmpIndices;
//   Complex* TmpInteractionFactor;
//   Complex TmpSum = 0.0;
//   int Index;
//   for (int j = 0; j < this->NbrIntraSectorSums; ++j)
//     {
//       int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
//       TmpIndices = this->IntraSectorIndicesPerSum[j];
//       for (int i1 = 0; i1 < Lim; i1 += 2)
// 	{
// 	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
// 	  if (Coefficient3 != 0.0)
// 	    {
// 	      TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
// 	      Coefficient4 = vSource[index];
// 	      Coefficient4 *= Coefficient3;
// 	      for (int i2 = 0; i2 < Lim; i2 += 2)
// 		{
// 		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
// 		  if (Index <= index)
// 		    {
// 		      if (Index < index)
// 			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
// 		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
// 		    }
// 		  ++TmpInteractionFactor;
// 		}
// 	    }
// 	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
// 	  if (Coefficient3 != 0.0)
// 	    {
// 	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
// 	      Coefficient4 = vSource[index];
// 	      Coefficient4 *= Coefficient3;
// 	      for (int i2 = 0; i2 < Lim; i2 += 2)
// 		{
// 		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
// 		  if (Index <= index)
// 		    {
// 		      if (Index < index)
// 			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
// 		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
// 		    }
// 		  ++TmpInteractionFactor;
// 		}
// 	    }
// 	}
//     }
//   for (int j = 0; j < this->NbrInterSectorSums; ++j)
//     {
//       int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
//       TmpIndices = this->InterSectorIndicesPerSum[j];
//       for (int i1 = 0; i1 < Lim; i1 += 2)
// 	{
// 	  Coefficient3 = particles->AuAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
// 	  if (Coefficient3 != 0.0)
// 	    {
// 	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim) >> 2]);
// 	      Coefficient4 = vSource[index];
// 	      Coefficient4 *= Coefficient3;
// 	      for (int i2 = 0; i2 < Lim; i2 += 2)
// 		{
// 		  Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
// 		  if (Index <= index)
// 		    {
// 		      if (Index < index)
// 			TmpSum += vSource[Index] * (Coefficient * Coefficient3) * Conj(*TmpInteractionFactor);
// 		      vDestination[Index] += (Coefficient * (*TmpInteractionFactor)) * Coefficient4;
// 		    }
// 		  ++TmpInteractionFactor;
// 		}
// 	    }
// 	} 
//     }
//   vDestination[index] += TmpSum;
// }

// core part of the AddMultiply method involving the two-body interaction for a set of vectors
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// tmpCoefficients = a temporary array whose size is nbrVectors

// inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::HermitianEvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, ComplexVector* vSources, 
// 													       ComplexVector* vDestinations, int nbrVectors, Complex* tmpCoefficients)
// {
//   //int Dim = particles->GetHilbertSpaceDimension();
//   double Coefficient;
//   double Coefficient3;
//   int Index;
//   
//   int* TmpIndices;
//   Complex* TmpInteractionFactor;
//   Complex* TmpSum = new Complex[nbrVectors];
//   for (int j = 0; j < this->NbrIntraSectorSums; ++j)
//     {
//       int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
//       TmpIndices = this->IntraSectorIndicesPerSum[j];
//       for (int i1 = 0; i1 < Lim; i1 += 2)
// 	{
// 	  Coefficient3 = particles->AuAu(index, TmpIndices[i1], TmpIndices[i1 + 1]);
// 	  if (Coefficient3 != 0.0)
// 	    {
// 	      TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
// 	      for (int p = 0; p < nbrVectors; ++p)
// 		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
// 	      for (int i2 = 0; i2 < Lim; i2 += 2)
// 		{
// 		  Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
// 		  if (Index <= index)
// 		    {
// 		      if (Index < index)
// 			{
// 			  for (int p = 0; p < nbrVectors; ++p)
// 			    {
// 			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
// 			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor)) * vSources[p][Index];
// 			    }
// 			}
// 		      else
// 			{
// 			  for (int p = 0; p < nbrVectors; ++p)
// 			    {
// 			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
// 			    }
// 			}
// 		    }
// 		  ++TmpInteractionFactor;
// 		}
// 	    }
// 	  Coefficient3 = particles->AdAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
// 	  if (Coefficient3 != 0.0)
// 	    {
// 	      TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
// 	      for (int p = 0; p < nbrVectors; ++p)
// 		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
// 	      for (int i2 = 0; i2 < Lim; i2 += 2)
// 		{
// 		  Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
// 		  if (Index <= index)
// 		    {
// 		      if (Index < index)
// 			{
// 			  for (int p = 0; p < nbrVectors; ++p)
// 			    {
// 			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
// 			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor)) * vSources[p][Index];
// 			    }
// 			}
// 		      else
// 			{
// 			  for (int p = 0; p < nbrVectors; ++p)
// 			    {
// 			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
// 			    }
// 			}
// 		    }
// 		  ++TmpInteractionFactor;
// 		}
// 	    }
// 	}
//     }
//   for (int j = 0; j < this->NbrInterSectorSums; ++j)
//     {
//       int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
//       TmpIndices = this->InterSectorIndicesPerSum[j];
//       for (int i1 = 0; i1 < Lim; i1 += 2)
// 	{
// 	  Coefficient3 = particles->AuAd(index, TmpIndices[i1], TmpIndices[i1 + 1]);
// 	  if (Coefficient3 != 0.0)
// 	    {
// 	      TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim) >> 2]);
// 	      for (int p = 0; p < nbrVectors; ++p)
// 		tmpCoefficients[p] = Coefficient3 * vSources[p][index];
// 	      for (int i2 = 0; i2 < Lim; i2 += 2)
// 		{
// 		  Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
// 		  if (Index <= index)
// 		    {
// 		      if (Index < index)
// 			{
// 			  for (int p = 0; p < nbrVectors; ++p)
// 			    {
// 			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
// 			      TmpSum[p] += (Coefficient * Coefficient3) * Conj((*TmpInteractionFactor)) * vSources[p][Index];
// 			    }
// 			}
// 		      else
// 			{
// 			  for (int p = 0; p < nbrVectors; ++p)
// 			    {
// 			      vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * tmpCoefficients[p];
// 			    }
// 			}
// 		    }
// 		  ++TmpInteractionFactor;
// 		}
// 	    }
// 	}
//     }
//   for (int l = 0; l < nbrVectors; ++l)
//     vDestinations[l][index] += TmpSum[l];
//   delete[] TmpSum;
// }

// core part of the FastMultiplication method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
													     int* indexArray, Complex* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  double* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();

  if (this->HermitianSymmetryFlag == false)
    {
      for (int j = 0; j < this->NbrSite; ++j)
	{
	  for (int k = 0; k < 3; ++k)
	  {
	    int j2 = this->MapNearestNeighborBonds[j][k];
	    if (j2 < this->NbrSite)
	    {
	      Coefficient2 = particles->AuAu(index + this->PrecalculationShift, j, j2);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupup[j][k]);
		  Index = particles->AduAdu(j, j2, Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
		      ++position;
		    }
		  ++TmpInteractionFactor;
		  
		  TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][k]);
		  Index = particles->AddAdd(j, j2, Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
		      ++position;
		    }
		  ++TmpInteractionFactor;
		  
		}
	  
	      Coefficient2 = particles->AdAd(index + this->PrecalculationShift, j, j2);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][k]);
		  Index = particles->AddAdd(j, j2, Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
		      ++position;
		    }
		   ++TmpInteractionFactor;
		  
		  TmpInteractionFactor = &(this->InteractionFactorsupup[j][k]);
		  Index = particles->AduAdu(j, j2, Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
		      ++position;
		    }
		   ++TmpInteractionFactor;
		}
	
		  
	      Coefficient2 = particles->AuAd(index + this->PrecalculationShift, j, j2);
	      if (Coefficient2 != 0.0)
		{
		  TmpInteractionFactor = &(this->InteractionFactorsupdown[j][k]);
		  Index = particles->AduAdd(j, j2, Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
		      ++position;
		    }
		  ++TmpInteractionFactor;
		  
		  TmpInteractionFactor = &(this->InteractionFactorsupdowndownup[j][k]);
		  Index = particles->AduAdd(j2, j, Coefficient);
		  if (Index < Dim)
		    {
		      indexArray[position] = Index;
		      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
		      ++position;
		    }
		  ++TmpInteractionFactor;
		 }
		 
// 	      Coefficient2 = particles->AdAu(index + this->PrecalculationShift, j, j2);
// 	      if (Coefficient2 != 0.0)
// 		{
// 		  TmpInteractionFactor = &(this->InteractionFactorsupdown[j][k]);
// 		  Index = particles->AddAdu(j, j2, Coefficient);
// 		  if (Index < Dim)
// 		    {
// 		      indexArray[position] = Index;
// 		      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
// 		      ++position;
// 		    }
// 		  ++TmpInteractionFactor;
// 		  
// 		  TmpInteractionFactor = &(this->InteractionFactorsupdowndownup[j][k]);
// 		  Index = particles->AduAdd(j, j2, Coefficient);
// 		  if (Index < Dim)
// 		    {
// 		      indexArray[position] = Index;
// 		      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
// 		      ++position;
// 		    }
// 		  ++TmpInteractionFactor;
// 		 }
		 
		}
	    }
	 if (this->UPotential != 0)
	 {
	   Coefficient2 = particles->AuAd(index, j, j);
	   if (Coefficient2 != 0)
	    {	  
	      TmpInteractionFactor = &(this->UPotential);
	      Index = particles->AduAdd(j, j, Coefficient);
	      if (Index < Dim)
	      {
		indexArray[position] = Index;
		coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
		++position;
	      }
		
	      ++TmpInteractionFactor;
	    }
	 }
	}  
  }
//   else
//     {
//       int AbsoluteIndex = index + this->PrecalculationShift;
//       for (int j = 0; j < this->NbrIntraSectorSums; ++j)
// 	{
// 	  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
// 	  TmpIndices = this->IntraSectorIndicesPerSum[j];
// 	  for (int i1 = 0; i1 < Lim; i1 += 2)
// 	    {
// 	      Coefficient2 = particles->AuAu(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
// 	      if (Coefficient2 != 0.0)
// 		{
// 		  TmpInteractionFactor = &(this->InteractionFactorsupup[j][(i1 * Lim) >> 2]);
// 		  for (int i2 = 0; i2 < Lim; i2 += 2)
// 		    {
// 		      Index = particles->AduAdu(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
// 		      if (Index <= AbsoluteIndex)
// 			{
// 			  if (Index == AbsoluteIndex)
// 			    {
// 			      indexArray[position] = Index;
// 			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (*TmpInteractionFactor);
// 			      ++position;
// 			    }
// 			  else
// 			    {
// 			      indexArray[position] = Index;
// 			      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
// 			      ++position;
// 			    }
// 			}
// 		      ++TmpInteractionFactor;
// 		    }
// 		}
// 	      Coefficient2 = particles->AdAd(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
// 	      if (Coefficient2 != 0.0)
// 		{
// 		  TmpInteractionFactor = &(this->InteractionFactorsdowndown[j][(i1 * Lim) >> 2]);
// 		  for (int i2 = 0; i2 < Lim; i2 += 2)
// 		    {
// 		      Index = particles->AddAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
// 		      if (Index <= AbsoluteIndex)
// 			{
// 			  if (Index == AbsoluteIndex)
// 			    {
// 			      indexArray[position] = Index;
// 			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (*TmpInteractionFactor);
// 			      ++position;
// 			    }
// 			  else
// 			    {
// 			      indexArray[position] = Index;
// 			      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
// 			      ++position;
// 			    }
// 			}
// 		      ++TmpInteractionFactor;
// 		    }
// 		}
// 	    }
// 	}
//       for (int j = 0; j < this->NbrInterSectorSums; ++j)
// 	{
// 	  int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
// 	  TmpIndices = this->InterSectorIndicesPerSum[j];
// 	  for (int i1 = 0; i1 < Lim; i1 += 2)
// 	    {
// 	      Coefficient2 = particles->AuAd(AbsoluteIndex, TmpIndices[i1], TmpIndices[i1 + 1]);
// 	      if (Coefficient2 != 0.0)
// 		{
// 		  TmpInteractionFactor = &(this->InteractionFactorsupdown[j][(i1 * Lim) >> 2]);
// 		  for (int i2 = 0; i2 < Lim; i2 += 2)
// 		    {
// 		      Index = particles->AduAdd(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
// 		      if (Index <= AbsoluteIndex)
// 			{
// 			  if (Index == AbsoluteIndex)
// 			    {
// 			      indexArray[position] = Index;
// 			      coefficientArray[position] = 0.5 * Coefficient * Coefficient2 * (*TmpInteractionFactor);
// 			      ++position;
// 			    }
// 			  else
// 			    {
// 			      indexArray[position] = Index;
// 			      coefficientArray[position] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
// 			      ++position;
// 			    }
// 			}
// 		      ++TmpInteractionFactor;
// 		    }
// 		}
// 	    }
// 	}  
//     }
}


// core part of the AddMultiply method involving the one-body interaction, including loop on vector components
// 
// particles = pointer to the Hilbert space
// firstComponent = first index vector to act on
// lastComponent = last index vector to act on (not included)
// step = step to go from one component to the other one
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
												      int step, ComplexVector& vSource, ComplexVector& vDestination)
{
  int Index;
	double Coefficient;
  if (this->OneBodyInteractionFactorsupup != 0)
    if (this->OneBodyInteractionFactorsdowndown != 0)
      {
	for (int i = firstComponent; i < lastComponent; i += step)
	  { 
	    for (int j = 0; j < this->NbrSite; ++j) 
	      {
		for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		  {
		    Index = particles->AduAu(i , j, j2, Coefficient);
		    if (Index < particles->GetHilbertSpaceDimension())
		      vDestination[Index] += Coefficient * this->OneBodyInteractionFactorsupup[j][k] * vSource[i];
		  
		    Index = particles->AddAd(i , j, j2, Coefficient);
		    if (Index < particles->GetHilbertSpaceDimension())
		      vDestination[Index] += Coefficient * this->OneBodyInteractionFactorsdowndown[j][k] * vSource[i];		  
		  }
		}
	      }
	  }
      }
    else
      {
	for (int i = firstComponent; i < lastComponent; i += step)
	  { 
	    Coefficient = 0.0;
	    for (int j = 0; j < this->NbrSite; ++j) 
	      {
		for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		  {
		    Index = particles->AduAu(i, j, j2, Coefficient);
		    if (Index < particles->GetHilbertSpaceDimension())
		      vDestination[Index] += Coefficient * this->OneBodyInteractionFactorsupup[j][k] * vSource[i];
		  }
		}
	      }
	  }
      }
  else
    if (this->OneBodyInteractionFactorsdowndown != 0)
      {
	for (int i = firstComponent; i < lastComponent; i += step)
	  { 
	    Coefficient = 0.0;
	    for (int j = 0; j < this->NbrSite; ++j) 
	      {
		for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		  {
		    Index = particles->AddAd(i, j, j2, Coefficient);
		    if (Index < particles->GetHilbertSpaceDimension())
		      vDestination[Index] += Coefficient * this->OneBodyInteractionFactorsdowndown[j][k] * vSource[i];		  
		  }
		}
	      }
	  }
      }	
  for (int i = firstComponent; i < lastComponent; i += step)
    vDestination[i] += this->HamiltonianShift * vSource[i];
  if (this->OneBodyInteractionFactorsupdown != 0)
    {
      double Coefficient;
      Complex Source;
      int Dim = particles->GetHilbertSpaceDimension();
      int Index;
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  Source = vSource[i];
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
	      {
		int j2 = this->MapNearestNeighborBonds[j][k];
		if (j2 < this->NbrSite)
		{
		  Index = particles->AddAu(i + this->PrecalculationShift, j, j2, Coefficient);
		  if (Index < Dim)
		      vDestination[Index] += (Coefficient * Conj(OneBodyInteractionFactorsupdown[j][k])) * Source;
		  Index = particles->AduAd(i + this->PrecalculationShift, j, j2, Coefficient);
		  if (Index < Dim)
		      vDestination[Index] += (Coefficient * OneBodyInteractionFactorsupdown[j][k]) * Source;
		}
	      }
	    }
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

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNOneBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent,
												      int step, ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int Index;
  if (this->OneBodyInteractionFactorsupup != 0) 
    if (this->OneBodyInteractionFactorsdowndown != 0)
      {
	for (int p = 0; p < nbrVectors; ++p)
	  {
	    ComplexVector& TmpSourceVector = vSources[p];
	    ComplexVector& TmpDestinationVector = vDestinations[p];
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		for (int j = 0; j < this->NbrSite; ++j) 
		  {
		    for (int k = 0; k < 3; ++k)
		    {
		      int j2 = this->MapNearestNeighborBonds[j][k];
		      if (j2 < this->NbrSite)
		      {
			Index = particles->AduAu(i, j, j2, Coefficient);
			if (Index < Dim)
			  TmpDestinationVector[Index] += Coefficient * this->OneBodyInteractionFactorsupup[j][k] * TmpSourceVector[i];
		      
			Index = particles->AddAd(i, j, j2, Coefficient);
			if (Index < Dim)
			  TmpDestinationVector[Index] += Coefficient * this->OneBodyInteractionFactorsdowndown[j][k] * TmpSourceVector[i];
		      }
		    }
		  }
	      }
	  }
      }
    else
      {
	for (int p = 0; p < nbrVectors; ++p)
	  {
	    ComplexVector& TmpSourceVector = vSources[p];
	    ComplexVector& TmpDestinationVector = vDestinations[p];
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		for (int j = 0; j < this->NbrSite; ++j) 
		{
		  for (int k = 0; k < 3; ++k)
		  {
		    int j2 = this->MapNearestNeighborBonds[j][k];
		    if (j2 < this->NbrSite)
		    {
		      Index = particles->AduAu(i, j, j2, Coefficient);
			if (Index < Dim)
			  TmpDestinationVector[Index] += Coefficient * this->OneBodyInteractionFactorsupup[j][k] * TmpSourceVector[i];
		    }
		  }
		}
	      }
	  }
      }
  else
    if (this->OneBodyInteractionFactorsdowndown != 0)
      {
	for (int p = 0; p < nbrVectors; ++p)
	  {
	    ComplexVector& TmpSourceVector = vSources[p];
	    ComplexVector& TmpDestinationVector = vDestinations[p];
	    for (int i = firstComponent; i < lastComponent; i += step)
	      { 
		for (int j = 0; j < this->NbrSite; ++j) 
		{
		  for (int k = 0; k < 3; ++k)
		  {
		    int j2 = this->MapNearestNeighborBonds[j][k];
		    if (j2 < this->NbrSite)
		    {
		      Index = particles->AddAd(i, j, j2, Coefficient);
			if (Index < Dim)
			  TmpDestinationVector[Index] += Coefficient * this->OneBodyInteractionFactorsdowndown[j][k] * TmpSourceVector[i];
		    }
		  }
		}
	      }
	  }
      }	
//   for (int p = 0; p < nbrVectors; ++p)
//     {
//       ComplexVector& TmpSourceVector = vSources[p];
//       ComplexVector& TmpDestinationVector = vDestinations[p];
//       for (int i = firstComponent; i < lastComponent; i += step)
// 	TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
//     }
  for (int p = 0; p < nbrVectors; ++p)
    {
      ComplexVector& TmpSourceVector = vSources[p];
      ComplexVector& TmpDestinationVector = vDestinations[p];
      for (int i = firstComponent; i < lastComponent; i += step)
	TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
    }
  if (this->OneBodyInteractionFactorsupdown != 0)
    {
      for (int i = firstComponent; i < lastComponent; i += step)
	{
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
		{
		  int j2 = this->MapNearestNeighborBonds[j][k];
		  if (j2 < this->NbrSite)
		  {
		    Index = particles->AddAu(i + this->PrecalculationShift, j, j2, Coefficient);
		    if (Index < Dim)
		      {
			for (int p = 0; p < nbrVectors; ++p)
			  vDestinations[p][Index] += Coefficient * Conj(OneBodyInteractionFactorsupdown[j][k]) * vSources[p][i];
		      }
		    Index = particles->AduAd(i + this->PrecalculationShift, j, j2, Coefficient);
		    if (Index < Dim)
		    {
		      for (int p = 0; p < nbrVectors; ++p)
			vDestinations[p][Index] += Coefficient * OneBodyInteractionFactorsupdown[j][k] * vSources[p][i];
		    }
		  }
	      }
	    }
	  }
      }

}

// core part of the FastMultiplication method involving the one-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNOneBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, 
													     int* indexArray, Complex* coefficientArray, long& position)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  int Index;
  if ((this->OneBodyInteractionFactorsdowndown != 0) || (this->OneBodyInteractionFactorsupup != 0))
    {
      for (int j = 0; j < this->NbrSite; ++j)
	{
	  for (int k = 0; k < 3; ++k)
	  {
	    int j2 = this->MapNearestNeighborBonds[j][k];
	    if (j2 < this->NbrSite)
	    {
	      Index = particles->AduAu(index + this->PrecalculationShift, j, j2, Coefficient);
	      if (Index < Dim)
	      {
		indexArray[position] = Index;
		coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsupup[j][k];
		++position;
	      }
	      Index = particles->AddAd(index + this->PrecalculationShift, j, j2, Coefficient);
	      if (Index < Dim)
	      {
		indexArray[position] = Index;
		coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsdowndown[j][k];
		++position;
	      }
	    }
	  }
	}
      }
   else
   {
     if (this->OneBodyInteractionFactorsupup != 0)
     {
       for (int j = 0; j < this->NbrSite; ++j)
	{
	  for (int k = 0; k < 3; ++k)
	  {
	    int j2 = this->MapNearestNeighborBonds[j][k];
	    if (j2 < this->NbrSite)
	    {
	      Index = particles->AduAu(index + this->PrecalculationShift, j, j2, Coefficient);
	      if (Index < Dim)
	      {
		indexArray[position] = Index;
		coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsupup[j][k];
		++position;
	      }
	    }
	  }
	}
     }
     else
       if (this->OneBodyInteractionFactorsdowndown != 0)
       {
	 for (int j = 0; j < this->NbrSite; ++j)
	{
	  for (int k = 0; k < 3; ++k)
	  {
	    int j2 = this->MapNearestNeighborBonds[j][k];
	    if (j2 < this->NbrSite)
	    {
	      Index = particles->AddAd(index + this->PrecalculationShift, j, j2, Coefficient);
	      if (Index < Dim)
	      {
		indexArray[position] = Index;
		coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsdowndown[j][k];
		++position;
	      }
	    }
	  }
	}
       } 
   }
  if (this->OneBodyInteractionFactorsupdown != 0)
    {
      for (int j = 0; j < this->NbrSite; ++j)
	{
	  for (int k = 0; k < 3; ++k)
	  {
	    int j2 = this->MapNearestNeighborBonds[j][k];
	    if (j2 < this->NbrSite)
	    {
	      Index = particles->AddAu(index + this->PrecalculationShift, j, j2, Coefficient);
	      if (Index < Dim)
	      {
		indexArray[position] = Index;
		coefficientArray[position] = Coefficient * this->OneBodyInteractionFactorsupdown[j][k];
		++position;
	      }
	      Index = particles->AduAd(index + this->PrecalculationShift, j, j2, Coefficient);
	      if (Index < Dim)
	      {
		indexArray[position] = Index;
		coefficientArray[position] = Coefficient * Conj(this->OneBodyInteractionFactorsupdown[j][k]);
		++position;
	      }
	    }
	  }
	}
    }
}


// core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndices;
  // double* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      for (int k = 0; k < 3; ++k)
	      {
		int j2 = this->MapNearestNeighborBonds[j][k];
		if (j2 < this->NbrSite)
		{
		  Coefficient2 = particles->AuAu(i, j, j2);
		  if (Coefficient2 != 0.0)
		    {
		      Index = particles->AduAdu(j, j2, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
			
		      Index = particles->AddAdd(j, j2, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		    
		  Coefficient2 = particles->AdAd(i, j, j2);
		  if (Coefficient2 != 0.0)
		    {
		      Index = particles->AddAdd(j, j2, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
			
		      Index = particles->AduAdu(j, j2, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		    
		    Coefficient2 = particles->AuAd(i, j, j2);
		    if (Coefficient2 != 0.0)
		    {
		      Index = particles->AduAdd(j, j2, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
			
		      Index = particles->AduAdd(j2, j, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		  }
		}
	    }
	}
    }
}



// core part of the PartialFastMultiplicationMemory method involving two-body term and one-body terms
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateMNOneBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;

  int Dim = particles->GetHilbertSpaceDimension();
  
  if (this->HermitianSymmetryFlag == false)
    {
      if (this->OneBodyInteractionFactorsupup != 0) 
      {
	for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		  {
		    int j2 = this->MapNearestNeighborBonds[j][k];
		    if (j2 < this->NbrSite)
		    {
		      Index = particles->AduAu(i, j, j2, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		  }
		}
	    }
      }
     
     if (this->OneBodyInteractionFactorsdowndown != 0) 
      {
	for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		  {
		    int j2 = this->MapNearestNeighborBonds[j][k];
		    if (j2 < this->NbrSite)
		    {
		      Index = particles->AddAd(i, j, j2, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		    }
		  }
		}
	    }
      }
      if (this->OneBodyInteractionFactorsupdown != 0)
	{
	  for (int i = firstComponent; i < lastComponent; ++i)
	    {
	      for (int j=0; j< this->NbrSite; ++j)
		{
		  for (int k = 0; k < 3; ++k)
		  {
		    int j2 = this->MapNearestNeighborBonds[j][k];
		    if (j2 < this->NbrSite)
		    {
		      Index = particles->AddAu(i, j, j2, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		      Index = particles->AduAd(i, j, j2, Coefficient);
		      if (Index < Dim)
			{
			  ++memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
			}
		      }
		    }
		}
	    }
	}
    }
}


inline void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::PlotMapNearestNeighborBonds()
{
  int Index;
  this->MapNearestNeighborBonds = new int* [this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
    this->MapNearestNeighborBonds[i] = new int [3];
  
  for (int i = 0; i <= (this->NbrSite - 2)/4 ; ++i)
  {
    Index = 4*i;
    if (i == 0)
      this->MapNearestNeighborBonds[Index][0] = this->NbrSite;
    else
      this->MapNearestNeighborBonds[Index][0] = Index - 2;
    if (i < (this->NbrSite - 2)/4)
      this->MapNearestNeighborBonds[Index][1] = Index + 2;
    else
      this->MapNearestNeighborBonds[Index][1] = this->NbrSite;
    this->MapNearestNeighborBonds[Index][2] = Index + 1;
    
    Index = 4*i + 1;
    if (i < (this->NbrSite - 2)/4)
      this->MapNearestNeighborBonds[Index][0] = Index + 2;
    else
      this->MapNearestNeighborBonds[Index][0] = this->NbrSite;
    if (i == 0)
      this->MapNearestNeighborBonds[Index][1] = this->NbrSite;
    else
      this->MapNearestNeighborBonds[Index][1] = Index - 2;
    this->MapNearestNeighborBonds[Index][2] = Index - 1;
    
    if (i < (this->NbrSite - 2)/4)
    {
      Index = 4*i + 2;
      this->MapNearestNeighborBonds[Index][0] = Index + 2;
      this->MapNearestNeighborBonds[Index][1] = Index - 2;
      this->MapNearestNeighborBonds[Index][2] = this->NbrSite;
      
      Index = 4*i + 3;
      this->MapNearestNeighborBonds[Index][0] = Index - 2;
      this->MapNearestNeighborBonds[Index][1] = Index + 2;
      this->MapNearestNeighborBonds[Index][2] = this->NbrSite;      
    }
  }
}

inline int ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::GetIndexNearestNeighborXBond(int i)
{
  return this->MapNearestNeighborBonds[i][0];
}

inline int ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::GetIndexNearestNeighborYBond(int i)
{
  return this->MapNearestNeighborBonds[i][1];
}

inline int ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::GetIndexNearestNeighborZBond(int i)
{
  return this->MapNearestNeighborBonds[i][2];
}
#endif