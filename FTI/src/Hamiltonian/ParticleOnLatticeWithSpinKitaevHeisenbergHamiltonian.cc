////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//                   class of Hubbard model hamiltonian associated            //
//                to particles on a lattice                                   //
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


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;



// default constructor
//

ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian()
{
  this->HermitianSymmetryFlag = false;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kineticFactor = multiplicative factor in front of the kinetic term
// uPotential = Hubbard potential strength
// bandParameter = band parameter
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSite, double kineticFactorIsotropic, double kineticFactorAnisotropic, double uPotential, double j1Factor, double j2Factor, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite - 1;
  this->UPotential = 2.0 * uPotential / ((double) this->NbrParticles);
  this->HamiltonianShift = 0.0;//4.0 * uPotential;
  this->KineticFactorIsotropic = kineticFactorIsotropic;
  this->KineticFactorAnisotropic = kineticFactorAnisotropic;
  this->J1Factor = j1Factor;
  this->J2Factor = j2Factor;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->PlotMapNearestNeighborBonds();
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = false;
  
  this->EvaluateInteractionFactors();
  
//   for (int i = 0; i < this->NbrSite; ++i)
//   {
//     for (int k = 0; k < 3; ++k)
//     {
//       cout << i << " " << k << " " << this->InteractionFactorsupup[i][k] << endl;
//     }
//   }
    
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1 << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1 << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::~ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian()
{
  if (this->MapNearestNeighborBonds != 0)
  {
    for (int i = 0; i < this->NbrSite; ++i)
      delete[] this->MapNearestNeighborBonds[i];
   delete[] this->MapNearestNeighborBonds; 
  }
  if (this->InteractionFactorsupup != 0)
    {
      for (int i = 0; i < this->NbrSite; ++i)
	{
	  delete[] this->InteractionFactorsupup[i];
	  delete[] this->InteractionFactorsdowndown[i];
	  delete[] this->InteractionFactorsupdown[i];
	  delete[] this->InteractionFactorsupdowndownup[i];
	  delete[] this->InteractionFactorsupupdowndown[i];
	}
      delete[] this->InteractionFactorsupup;
      delete[] this->InteractionFactorsupdown;
      delete[] this->InteractionFactorsdowndown;
      delete[] this->InteractionFactorsupdowndownup;
      delete[] this->InteractionFactorsupupdowndown;
    }
  if (this->OneBodyInteractionFactorsupup != 0)
    {
      for (int i = 0; i < this->NbrSite; ++i)
	{
	  delete[] this->OneBodyInteractionFactorsupup[i];
	  delete[] this->OneBodyInteractionFactorsdowndown[i];
	  delete[] this->OneBodyInteractionFactorsupdown[i];
	}
	  
      delete[] this->OneBodyInteractionFactorsupup;
      delete[] this->OneBodyInteractionFactorsdowndown;
      delete[] this->OneBodyInteractionFactorsupdown;
    }
  
    
  if (this->FastMultiplicationFlag == true)
    {
      long MinIndex;
      long MaxIndex;
      this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
      int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
      this->FastMultiplicationFlag = false;
    }
}
  
// ask if Hamiltonian implements hermitian symmetry operations
//

bool ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnSphereWithSpin*) hilbertSpace;
  this->EvaluateInteractionFactors();
}
 
// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  Complex Coefficient;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      Coefficient = vSource[k];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      vDestination[k++] += this->HamiltonianShift * Coefficient;
	    }
	}
      else
	{
	  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  Complex Coefficient;
	  int j;
	  int TmpNbrInteraction;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  int Pos = firstComponent / this->FastMultiplicationStep; 
	  int PosMod = firstComponent % this->FastMultiplicationStep;
	  if (PosMod != 0)
	    {
	      ++Pos;
	      PosMod = this->FastMultiplicationStep - PosMod;
	    }
	  int l =  PosMod + firstComponent + this->PrecalculationShift;
	  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
	      Coefficient = vSource[l];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      vDestination[l] += this->HamiltonianShift * Coefficient;
	      l += this->FastMultiplicationStep;
	      ++Pos;
	    }
	  firstComponent += this->PrecalculationShift;
	  LastComponent += this->PrecalculationShift;
	  for (l = 0; l < this->FastMultiplicationStep; ++l)
	    if (PosMod != l)
	      {	
		for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
		  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
		this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
	      }
	  delete TmpParticles;
	}
    }
  return vDestination;
}
 
// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
											       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      Complex* Coefficient2 = new Complex [nbrVectors];
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
      delete[] Coefficient2;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  Complex Coefficient;
	  Complex* Coefficient2 = new Complex [nbrVectors];
	  Complex* TmpCoefficientArray; 
	  int j;
	  int Pos;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  Coefficient2[l] = vSources[l][k];
		  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
      else
	{
	  this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	}
    }
  return vDestinations;
}

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

ComplexVector* ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
												   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int* TmpIndexArray;
  Complex* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int Pos2;
  Complex Coefficient;
  int PosMod = firstComponent % this->FastMultiplicationStep;
  Complex* Coefficient2 = new Complex [nbrVectors];
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
      for (int k = 0; k < nbrVectors; ++k)
	{
	  Coefficient2[k] = vSources[k][l];
	  vDestinations[k][l] += this->HamiltonianShift * Coefficient2[k];
	}
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  Pos2 = TmpIndexArray[j];
	  Coefficient = TmpCoefficientArray[j];
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient2[k];
	}
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (l = 0; l < this->FastMultiplicationStep; ++l)
    if (PosMod != l)
      {	
	for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
	  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSources, vDestinations, nbrVectors);
      }
  delete[] Coefficient2;
  delete TmpParticles;
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

// ComplexVector& ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
// 												int firstComponent, int nbrComponent)
// {
//   int LastComponent = firstComponent + nbrComponent;
//   if (this->FastMultiplicationFlag == false)
//     {
//       ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
//       for (int i = firstComponent; i < LastComponent; ++i)
// 	this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
//       this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
//       delete TmpParticles;
//     }
//   else
//     {
//       if (this->FastMultiplicationStep == 1)
// 	{
// 	  int* TmpIndexArray;
// 	  Complex* TmpCoefficientArray; 
// 	  int j;
// 	  int TmpNbrInteraction;
// 	  int k = firstComponent;
// 	  Complex Coefficient;
// 	  firstComponent -= this->PrecalculationShift;
// 	  LastComponent -= this->PrecalculationShift;
// 	  for (int i = firstComponent; i < LastComponent; ++i)
// 	    {
// 	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
// 	      TmpIndexArray = this->InteractionPerComponentIndex[i];
// 	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
// 	      Coefficient = vSource[k];
// 	      Complex TmpSum = 0.0;
// 	      for (j = 0; j < TmpNbrInteraction; ++j)
// 		{
// 		  TmpSum += Conj(TmpCoefficientArray[j]) * vSource[TmpIndexArray[j]];
// 		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
// 		}
// 	      TmpSum += this->HamiltonianShift * Coefficient;
// 	      vDestination[k++] += TmpSum;
// 	    }
// 	}
//       else
// 	{
// 	  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
// 	  int* TmpIndexArray;
// 	  Complex* TmpCoefficientArray; 
// 	  Complex Coefficient;
// 	  int j;
// 	  int TmpNbrInteraction;
// 	  firstComponent -= this->PrecalculationShift;
// 	  LastComponent -= this->PrecalculationShift;
// 	  int Pos = firstComponent / this->FastMultiplicationStep; 
// 	  int PosMod = firstComponent % this->FastMultiplicationStep;
// 	  if (PosMod != 0)
// 	    {
// 	      ++Pos;
// 	      PosMod = this->FastMultiplicationStep - PosMod;
// 	    }
// 	  int l =  PosMod + firstComponent + this->PrecalculationShift;
// 	  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
// 	    {
// 	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
// 	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
// 	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
// 	      Coefficient = vSource[l];
// 	      Complex TmpSum = 0.0;
// 	      for (j = 0; j < TmpNbrInteraction; ++j)
// 		{
// 		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
// 		  TmpSum += Conj(TmpCoefficientArray[j]) * vSource[TmpIndexArray[j]];
// 		}
// 	      TmpSum += this->HamiltonianShift * Coefficient;	      
// 	      vDestination[l] += TmpSum;
// 	      l += this->FastMultiplicationStep;
// 	      ++Pos;
// 	    }
// 	  firstComponent += this->PrecalculationShift;
// 	  LastComponent += this->PrecalculationShift;
// 	  for (l = 0; l < this->FastMultiplicationStep; ++l)
// 	    if (PosMod != l)
// 	      {	
// 		for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
// 		  this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
// 		this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
// 	      }
// 	  delete TmpParticles;
// 	}
//     }
//   return vDestination;
// }

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

// ComplexVector* ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
// 													int firstComponent, int nbrComponent)
// {
//   int LastComponent = firstComponent + nbrComponent;
//   if (this->FastMultiplicationFlag == false)
//     {
//       Complex* Coefficient2 = new Complex [nbrVectors];
//       ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
//       for (int i = firstComponent; i < LastComponent; ++i)
// 	this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
//       this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
//       delete[] Coefficient2;
//       delete TmpParticles;
//     }
//   else
//     {
//       if (this->FastMultiplicationStep == 1)
// 	{
// 	  int* TmpIndexArray;
// 	  Complex Coefficient;
// 	  Complex* Coefficient2 = new Complex [nbrVectors];
// 	  Complex* TmpCoefficientArray; 
// 	  int j;
// 	  int Pos;
// 	  int TmpNbrInteraction;
// 	  int k = firstComponent;
// 	  firstComponent -= this->PrecalculationShift;
// 	  LastComponent -= this->PrecalculationShift;
// 	  Complex* TmpSum = new Complex [nbrVectors];
// 	  for (int i = firstComponent; i < LastComponent; ++i)
// 	    {
// 	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
// 	      TmpIndexArray = this->InteractionPerComponentIndex[i];
// 	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
// 	      for (int l = 0; l < nbrVectors; ++l)
// 		{
// 		  TmpSum[l] = 0.0;
// 		  Coefficient2[l] = vSources[l][k];
// 		}
// 	      for (j = 0; j < TmpNbrInteraction; ++j)
// 		{
// 		  Pos = TmpIndexArray[j];
// 		  Coefficient = TmpCoefficientArray[j];
// 		  for (int l = 0; l < nbrVectors; ++l)
// 		    {
// 		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
// 		      TmpSum[l] += Conj(Coefficient) * vSources[l][Pos];
// 		    }
// 		}
// 	      for (int l = 0; l < nbrVectors; ++l)
// 		{
// 		  TmpSum[l] += this->HamiltonianShift * Coefficient2[l];
// 		  vDestinations[l][k] += TmpSum[l];
// 		}
// 	      ++k;
// 	    }
// 	  delete[] Coefficient2;
// 	  delete[] TmpSum;
// 	}
//       else
// 	{
// 	  this->HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
// 	}
//     }
//   return vDestinations;
// }

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

// ComplexVector* ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
// 																     int firstComponent, int nbrComponent)
// {
//   int LastComponent = firstComponent + nbrComponent;
//   ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
//   int* TmpIndexArray;
//   Complex* TmpCoefficientArray; 
//   int j;
//   int TmpNbrInteraction;
//   firstComponent -= this->PrecalculationShift;
//   LastComponent -= this->PrecalculationShift;
//   int Pos = firstComponent / this->FastMultiplicationStep; 
//   int Pos2;
//   Complex Coefficient;
//   int PosMod = firstComponent % this->FastMultiplicationStep;
//   Complex* Coefficient2 = new Complex [nbrVectors];
//   Complex* TmpSum = new Complex [nbrVectors];
//   if (PosMod != 0)
//     {
//       ++Pos;
//       PosMod = this->FastMultiplicationStep - PosMod;
//     }
//   int l =  PosMod + firstComponent + this->PrecalculationShift;
//   for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
//     {
//       TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
//       TmpIndexArray = this->InteractionPerComponentIndex[Pos];
//       TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
//       for (int k = 0; k < nbrVectors; ++k)
// 	{
// 	  TmpSum[k] = 0.0;
// 	  Coefficient2[k] = vSources[k][l];
// 	}
//       for (j = 0; j < TmpNbrInteraction; ++j)
// 	{
// 	  Pos2 = TmpIndexArray[j];
// 	  Coefficient = TmpCoefficientArray[j];
// 	  for (int k = 0; k < nbrVectors; ++k)
// 	    {
// 	      vDestinations[k][TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient2[k];
// 	      TmpSum[k] += Conj(Coefficient) * vSources[k][Pos2];
// 	    }
// 	}
//      for (int k = 0; k < nbrVectors; ++k)
//        {
// 	 TmpSum[k] += this->HamiltonianShift * Coefficient2[k];
// 	 vDestinations[k][l] += TmpSum[k];
//        }
//       l += this->FastMultiplicationStep;
//       ++Pos;
//     }
//   firstComponent += this->PrecalculationShift;
//   LastComponent += this->PrecalculationShift;
//   for (l = 0; l < this->FastMultiplicationStep; ++l)
//     if (PosMod != l)
//       {	
// 	for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
// 	  this->HermitianEvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
// 	this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSources, vDestinations, nbrVectors);
//       }
//   delete[] TmpSum;
//   delete[] Coefficient2;
//   delete TmpParticles;
//   return vDestinations;
// }

// evaluate all nearest neighbor interaction factors
//   

void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  
  this->InteractionFactorsupup = new double* [this->NbrSite];
  this->InteractionFactorsdowndown = new double* [this->NbrSite];
  this->InteractionFactorsupdown = new double* [this->NbrSite];
  this->InteractionFactorsupupdowndown = new double* [this->NbrSite];
  this->InteractionFactorsupdowndownup = new double* [this->NbrSite];
  
  this->OneBodyInteractionFactorsupup = new double* [this->NbrSite];
  this->OneBodyInteractionFactorsdowndown = new double* [this->NbrSite];
  this->OneBodyInteractionFactorsupdown = new Complex* [this->NbrSite];
  
  for (int i = 0; i < this->NbrSite; ++i)
  {
   this->InteractionFactorsupup[i]  = new double [3];
   this->InteractionFactorsdowndown[i]  = new double [3];
   this->InteractionFactorsupdown[i]  = new double [3];
   this->InteractionFactorsupupdowndown[i] = new double [3];
   this->InteractionFactorsupdowndownup[i] = new double [3];
   
   this->OneBodyInteractionFactorsupup[i] = new double [3];
   this->OneBodyInteractionFactorsdowndown[i] = new double [3];
   this->OneBodyInteractionFactorsupdown[i] = new Complex [3];
   for (int j = 0; j < 3; ++j)
   {
    this->InteractionFactorsupup[i][j]  = 0.0; 
    this->InteractionFactorsdowndown[i][j]  = 0.0; 
    this->InteractionFactorsupdown[i][j]  = 0.0; 
    this->InteractionFactorsupupdowndown[i][j] = 0.0;
    this->InteractionFactorsupdowndownup[i][j] = 0.0;
    
    this->OneBodyInteractionFactorsupup[i][j] = 0.0;
    this->OneBodyInteractionFactorsdowndown[i][j] = 0.0;
    this->OneBodyInteractionFactorsupdown[i][j] = 0.0;
   }
  }
  
  for (int i = 0; i < this->NbrSite; ++i)
  {
   int jx = this->MapNearestNeighborBonds[i][0];
   int jy = this->MapNearestNeighborBonds[i][1];
   int jz = this->MapNearestNeighborBonds[i][2];
   
   if (jx < this->NbrSite)
   {
    this->OneBodyInteractionFactorsupup[i][0] += -this->KineticFactorIsotropic;
    this->OneBodyInteractionFactorsdowndown[i][0] += -this->KineticFactorIsotropic;
    this->OneBodyInteractionFactorsupdown[i][0] += Complex(-this->KineticFactorIsotropic - this->KineticFactorAnisotropic, 0) * 2;
    
    this->InteractionFactorsupup[i][0] += -(this->J1Factor - this->J2Factor);
    this->InteractionFactorsdowndown[i][0] += -(this->J1Factor - this->J2Factor);
    this->InteractionFactorsupdown[i][0] += (this->J1Factor - this->J2Factor) * 2;
    this->InteractionFactorsupupdowndown[i][0] += -2*this->J1Factor;
    this->InteractionFactorsupdowndownup[i][0] += 2*this->J2Factor * 2;
    
    TotalNbrInteractionFactors += 8;
   }
   
   if (jy < this->NbrSite)
   {
     this->OneBodyInteractionFactorsupup[i][1] += -this->KineticFactorIsotropic;
     this->OneBodyInteractionFactorsdowndown[i][1] += -this->KineticFactorIsotropic;
     this->OneBodyInteractionFactorsupdown[i][1] += Complex(-this->KineticFactorIsotropic,  this->KineticFactorAnisotropic) * 2;
     
     this->InteractionFactorsupup[i][0] += -(this->J1Factor - this->J2Factor);
     this->InteractionFactorsdowndown[i][0] += -(this->J1Factor - this->J2Factor);
     this->InteractionFactorsupdown[i][0] += (this->J1Factor - this->J2Factor) * 2;
     this->InteractionFactorsupupdowndown[i][0] += -2*this->J1Factor;
     this->InteractionFactorsupdowndownup[i][0] += -2*this->J2Factor * 2;
     
     TotalNbrInteractionFactors += 8;
   }
   
   if (jz < this->NbrSite)
   {
    this->OneBodyInteractionFactorsupup[i][2] += -this->KineticFactorIsotropic - this->KineticFactorAnisotropic;
    this->OneBodyInteractionFactorsdowndown[i][2] += -this->KineticFactorIsotropic + this->KineticFactorAnisotropic;
    
    this->InteractionFactorsupup[i][0] += -(this->J1Factor + this->J2Factor);
    this->InteractionFactorsdowndown[i][0] += -(this->J1Factor + this->J2Factor);
    this->InteractionFactorsupdown[i][0] += (this->J1Factor + this->J2Factor) * 2;
    this->InteractionFactorsupupdowndown[i][0] += -2*(this->J1Factor - this->J2Factor);
    
    TotalNbrInteractionFactors += 8;
   }
//    cout << i << " " << this->InteractionFactorsupup[i][0] << " " << this->InteractionFactorsupup[i][1] << " " << this->InteractionFactorsupup[i][2] << endl;
  }
 
  
  
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  
  this->NbrInteractionPerComponent = new int [EffectiveHilbertSpaceDimension];
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    this->NbrInteractionPerComponent[i] = 0;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);

  long Memory = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    Memory += this->NbrInteractionPerComponent[i];

  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(Complex*)) * EffectiveHilbertSpaceDimension;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(Complex)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(Complex)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(Complex*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(Complex*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(Complex)));
      long ResidualMemory = allowedMemory - Memory;
      if (ResidualMemory > 0)
	{
	  int TotalReducedSpaceDimension = ReducedSpaceDimension;
	  int* TmpNbrInteractionPerComponent = new int [TotalReducedSpaceDimension];
	  int i = 0;
	  int Pos = 0;
	  for (; i < ReducedSpaceDimension; ++i)
	    {
	      TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[Pos];
	      Pos += this->FastMultiplicationStep;
	    }
	  delete[] this->NbrInteractionPerComponent;
	      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
	}
      else
	{
	  int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
	  delete[] this->NbrInteractionPerComponent;
	  this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
	}
    }
  else
    {
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(Complex*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(Complex)));
      this->FastMultiplicationStep = 1;
    }

  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  return Memory;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;

  this->EvaluateMNOneBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  this->EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);

  delete TmpParticles;

  return Memory;
}

// enable fast multiplication algorithm
//

void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new Complex* [ReducedSpaceDimension];
  //ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles;
  int TotalPos = 0;
  
  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
    {
      this->InteractionPerComponentIndex[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
      this->InteractionPerComponentCoefficient[TotalPos] = new Complex [this->NbrInteractionPerComponent[TotalPos]];      
//       TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
//       TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];
//       Pos = 0l;
//       this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, Pos);
//       this->EvaluateMNOneBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, Pos);
       ++TotalPos;
    }

  QHEParticlePrecalculationOperation Operation(this, false);
  Operation.ApplyOperation(this->Architecture);

  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = index of the last component that has to be precalcualted

void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  int LastComponent = nbrComponent + firstComponent;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();

  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  long Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      long TotalPos = 0;
     
      this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i + this->PrecalculationShift, this->InteractionPerComponentIndex[Pos], 
							 this->InteractionPerComponentCoefficient[Pos], TotalPos);
      
      this->EvaluateMNOneBodyFastMultiplicationComponent(TmpParticles, i + this->PrecalculationShift, this->InteractionPerComponentIndex[Pos], 
      							 this->InteractionPerComponentCoefficient[Pos], TotalPos);
      ++Pos;
    }
  delete TmpParticles;
}


