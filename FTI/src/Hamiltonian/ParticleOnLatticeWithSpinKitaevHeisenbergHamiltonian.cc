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
  this->OneBodyGenericInteractionFactorsupup = 0;
  this->OneBodyGenericInteractionFactorsdowndown = 0;
  this->OneBodyGenericInteractionFactorsupdown = 0;
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
//   if (this->InteractionFactorsupup != 0)
//     {
//       for (int i = 0; i < this->NbrSite; ++i)
// 	{
// 	  delete[] this->InteractionFactorsupup[i];
// 	  delete[] this->InteractionFactorsdowndown[i];
// 	  delete[] this->InteractionFactorsupdown[i];
// 	  delete[] this->InteractionFactorsupupdowndown[i];
// 	}
//       delete[] this->InteractionFactorsupup;
//       delete[] this->InteractionFactorsupdown;
//       delete[] this->InteractionFactorsdowndown;
//       delete[] this->InteractionFactorsupupdowndown;
//     }
  if (this->OneBodyGenericInteractionFactorsupup != 0)
    {
      for (int i = 0; i < this->NbrSite; ++i)
	{
	  delete[] this->OneBodyGenericInteractionFactorsupup[i];
	  delete[] this->OneBodyGenericInteractionFactorsdowndown[i];
	  delete[] this->OneBodyGenericInteractionFactorsupdown[i];
	}
	  
      delete[] this->OneBodyGenericInteractionFactorsupup;
      delete[] this->OneBodyGenericInteractionFactorsdowndown;
      delete[] this->OneBodyGenericInteractionFactorsupdown;
    }
  
    
//   if (this->FastMultiplicationFlag == true)
//     {
//       long MinIndex;
//       long MaxIndex;
//       this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
//       int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
//       int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
//       if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
// 	++ReducedDim;
//       for (int i = 0; i < ReducedDim; ++i)
// 	{
// 	  delete[] this->InteractionPerComponentIndex[i];
// 	  delete[] this->InteractionPerComponentCoefficient[i];
// 	}
//       delete[] this->InteractionPerComponentIndex;
//       delete[] this->InteractionPerComponentCoefficient;
//       delete[] this->NbrInteractionPerComponent;
//       this->FastMultiplicationFlag = false;
//     }
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
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  
  this->NbrInterSectorSums = 0;
  for (int j = 0; j < this->NbrSite; ++j)
    for (int k = 0; k < 3; ++k)
      if (this->MapNearestNeighborBonds[j][k] < this->NbrSite)
	++this->NbrInterSectorSums;
   
  
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  this->NbrIntraSectorSums = this->NbrInterSectorSums;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      
  
  int index = 0;
  for (int j = 0; j < this->NbrSite; ++j)
    for (int k = 0; k < 3; ++k)  
      if (this->MapNearestNeighborBonds[j][k] < this->NbrSite)
	{
	  ++this->NbrInterSectorIndicesPerSum[index];
	  ++index;
	}
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      if (this->NbrInterSectorIndicesPerSum[i] > 0)
	{
	  this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	  this->NbrInterSectorIndicesPerSum[i] = 0;
	}
    }
  int TmpSum = 0;
  for (int j = 0; j < this->NbrSite; ++j)
    for (int k = 0; k < 3; ++k)
      if (this->MapNearestNeighborBonds[j][k] < this->NbrSite)
	  {
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = j;
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = this->MapNearestNeighborBonds[j][k];
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	    ++TmpSum;
	  }
 
  TmpSum = 0;
  for (int j = 0; j < this->NbrSite; ++j)
    for (int k = 0; k < 3; ++k)
      {
	int Index1 = j;
	int Index2 = this->MapNearestNeighborBonds[j][k];
	if (Index2 < this->NbrSite)
	{
	  if (Index1 < Index2)
	    ++this->NbrIntraSectorIndicesPerSum[TmpSum];
	  ++TmpSum;
	}
      }
  this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    {
      if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
    }
  TmpSum = 0;
  for (int j = 0;j < this->NbrSite; ++j)
    for (int k = 0; k < 3; ++k)
      {
	int Index1 = j;
	int Index2 = this->MapNearestNeighborBonds[j][k];
	if (Index2 < this->NbrSite)
	  {
	    if (Index1 < Index2)
	      {
		this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		++this->NbrIntraSectorIndicesPerSum[TmpSum];    
	      }
	    ++TmpSum;
	  }
      }
  this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
  this->InteractionFactorsupupupdown = new Complex* [this->NbrInterSectorSums];
  this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
  this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
  this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
  this->InteractionFactorsdowndownupdown = new Complex* [this->NbrInterSectorSums];
  this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
  this->InteractionFactorsupdownupdown = new Complex* [this->NbrInterSectorSums];
  this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];
  
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    {
      this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
      this->InteractionFactorsdowndowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
      this->InteractionFactorsupupdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
      this->InteractionFactorsdowndownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
      this->InteractionFactorsupdownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
      this->InteractionFactorsupdowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
      int Index = 0;
      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	{
	  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      if (this->FindBondType(Index1, Index2) == 0)
		{
		  this->InteractionFactorsupupupup[i][Index] = -(this->J1Factor - this->J2Factor);
		  this->InteractionFactorsdowndowndowndown[i][Index] = -(this->J1Factor - this->J2Factor);
		  this->InteractionFactorsdowndownupup[i][Index] = -2*this->J1Factor;
		  this->InteractionFactorsupupdowndown[i][Index] = -2*this->J1Factor;
		}
	      
	      if (this->FindBondType(Index1, Index2) == 1)
		{
		  this->InteractionFactorsupupupup[i][Index] = -(this->J1Factor - this->J2Factor);
		  this->InteractionFactorsdowndowndowndown[i][Index] = -(this->J1Factor - this->J2Factor);
		  this->InteractionFactorsdowndownupup[i][Index] = (this->J1Factor - this->J2Factor);
		  this->InteractionFactorsupupdowndown[i][Index] = (this->J1Factor - this->J2Factor);
		}
	      
	      if (this->FindBondType(Index1, Index2) == 2)
		{
		  this->InteractionFactorsupupupup[i][Index] = -(this->J1Factor + this->J2Factor);
		  this->InteractionFactorsdowndowndowndown[i][Index] = -(this->J1Factor + this->J2Factor);
		  this->InteractionFactorsdowndownupup[i][Index] = -2*(this->J1Factor - this->J2Factor);
		  this->InteractionFactorsupupdowndown[i][Index] = -2*(this->J1Factor - this->J2Factor);
		}	      	      
	      TotalNbrInteractionFactors += 4;
	      ++Index;
	    }
	}
      Index = 0;
      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	{
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      this->InteractionFactorsupdownupup[i][Index] = 0.0;
	      this->InteractionFactorsupdowndowndown[i][Index] = 0.0;
	      TotalNbrInteractionFactors += 2;	      
	      ++Index;
	    }
	}
    }
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      this->InteractionFactorsupdownupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
      this->InteractionFactorsupupupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
      this->InteractionFactorsdowndownupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
      int Index = 0;
      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	{
	  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      
	      if (this->FindBondType(Index1, Index2) == 0)
		  {
		    if (Index1 == Index3)
		      this->InteractionFactorsupdownupdown[i][Index] = (this->J1Factor - this->J2Factor);
		    else
		      this->InteractionFactorsupdownupdown[i][Index] = 2*this->J2Factor;
		  }
	      if (this->FindBondType(Index1, Index2) == 1)
		{
		  if (Index1 == Index3)
		    this->InteractionFactorsupdownupdown[i][Index] = (this->J1Factor - this->J2Factor);
		  else
		    this->InteractionFactorsupdownupdown[i][Index] = -2*this->J2Factor;
		}
	      if (this->FindBondType(Index1, Index2) == 2)
		{		    
		  this->InteractionFactorsupdownupdown[i][Index] = (this->J1Factor + this->J2Factor);    
		}
	      TotalNbrInteractionFactors += 1;
	      ++Index;
	    }
	}
      Index = 0;
      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	{
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      this->InteractionFactorsupupupdown[i][Index] = 0.0;
	      this->InteractionFactorsdowndownupdown[i][Index] = 0.0;
	      TotalNbrInteractionFactors += 2;	      
	      ++Index;
	    }
	}
    }
  
  
  this->OneBodyGenericInteractionFactorsupup = new double*[this->NbrSite];
  this->OneBodyGenericInteractionFactorsdowndown = new double*[this->NbrSite];
  this->OneBodyGenericInteractionFactorsupdown = new Complex*[this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
    {
      this->OneBodyGenericInteractionFactorsupup[i] = new double [3];
      this->OneBodyGenericInteractionFactorsdowndown[i] = new double [3];
      this->OneBodyGenericInteractionFactorsupdown[i] = new Complex [3];
      for (int j = 0; j < 3; ++j)
	{   
	  this->OneBodyGenericInteractionFactorsupup[i][j] = 0.0;
	  this->OneBodyGenericInteractionFactorsdowndown[i][j] = 0.0;
	  this->OneBodyGenericInteractionFactorsupdown[i][j] = 0.0;
	}
    }
  
  for (int i = 0; i < this->NbrSite; ++i)
    {
      int jx = this->MapNearestNeighborBonds[i][0];
      int jy = this->MapNearestNeighborBonds[i][1];
      int jz = this->MapNearestNeighborBonds[i][2];
      
      if (jx < this->NbrSite)
	{
	  this->OneBodyGenericInteractionFactorsupup[i][0] += -this->KineticFactorIsotropic;
	  this->OneBodyGenericInteractionFactorsdowndown[i][0] += -this->KineticFactorIsotropic;
	  this->OneBodyGenericInteractionFactorsupdown[i][0] += Complex(-this->KineticFactorIsotropic - this->KineticFactorAnisotropic, 0) * 2;
	  
	}
      
      if (jy < this->NbrSite)
	{
	  this->OneBodyGenericInteractionFactorsupup[i][1] += -this->KineticFactorIsotropic;
	  this->OneBodyGenericInteractionFactorsdowndown[i][1] += -this->KineticFactorIsotropic;
	  this->OneBodyGenericInteractionFactorsupdown[i][1] += Complex(-this->KineticFactorIsotropic,  this->KineticFactorAnisotropic) * 2;
	  
	}
      
      if (jz < this->NbrSite)
	{
	  this->OneBodyGenericInteractionFactorsupup[i][2] += -this->KineticFactorIsotropic - this->KineticFactorAnisotropic;
	  this->OneBodyGenericInteractionFactorsdowndown[i][2] += -this->KineticFactorIsotropic + this->KineticFactorAnisotropic;
	  
	}
      //    cout << i << " " << this->InteractionFactorsupup[i][0] << " " << this->InteractionFactorsupup[i][1] << " " << this->InteractionFactorsupup[i][2] << endl;
    }
  
  
  
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}



