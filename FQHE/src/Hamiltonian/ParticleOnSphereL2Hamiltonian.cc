////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere where       //
//   the hamiltonian is reduced to a simple total square angular momentum     //
//                                                                            //
//                        last modification : 15/02/2007                      //
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


#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// totalLz = twice the projected momentum total value
// architecture = architecture to use for precalculation
// l2Factor = multiplicative factor in front of the L^2 operator in the Hamiltonian
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereL2Hamiltonian::ParticleOnSphereL2Hamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int totalLz,
								       AbstractArchitecture* architecture, double l2Factor, long memory, bool onDiskCacheFlag,
								       char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->TotalLz = totalLz;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = true;
  this->L2Factor = l2Factor;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.25 * this->L2Factor * ((double) (this->TotalLz * this->TotalLz));
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  if (precalculationFileName == 0)
    {
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
	  if (this->DiskStorageFlag == false)
	    {
	      this->EnableFastMultiplication();
	    }
	  else
	    {
	      char* TmpFileName = this->Architecture->GetTemporaryFileName();
	      this->EnableFastMultiplicationWithDiskStorage(TmpFileName);	      
	      delete[] TmpFileName;
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
  this->L2Operator = 0;
}

// destructor
//

ParticleOnSphereL2Hamiltonian::~ParticleOnSphereL2Hamiltonian() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  if (this->FastMultiplicationFlag == true)
    {
       if (this->DiskStorageFlag == false)
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
	}
       else
	 {
	  remove (this->DiskStorageFileName);
	  delete[] this->DiskStorageFileName;
	 }
       delete[] this->NbrInteractionPerComponent;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereL2Hamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereL2Hamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereL2Hamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereL2Hamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift + 0.25 * this->L2Factor * ((double) (this->TotalLz * this->TotalLz));
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereL2Hamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Particles->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereL2Hamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ParticleOnSphereL2Hamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ParticleOnSphereL2Hamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   

void ParticleOnSphereL2Hamiltonian::EvaluateInteractionFactors()
{
  RealMatrix Coefficients (this->LzMax + 1, this->LzMax + 1);
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax + 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	Coefficients(i, j) = TmpCoefficient;
    }
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax - 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	Coefficients(j, i) *= 0.5 * TmpCoefficient;
    }
  double Factor = 2.0 * this->L2Factor;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Factor *= -1.0;
  this->NbrInteractionFactors = this->LzMax * this->LzMax;
//   this->NbrM12Indices = this->NbrInteractionFactors;
//   this->NbrM3Values = new int[this->NbrInteractionFactors];
//   this->M3Values = new int* [this->NbrInteractionFactors];
  this->M1Value = new int [this->NbrInteractionFactors];
  this->M2Value = new int [this->NbrInteractionFactors];
  this->M3Value = new int [this->NbrInteractionFactors];
  this->InteractionFactors = new double [this->NbrInteractionFactors];
  this->NbrInteractionFactors = 0;
  for (int m4 = 0; m4 < this->LzMax; ++m4)
    for (int m3 = 1; m3 <= this->LzMax; ++m3)
      {
	this->InteractionFactors[this->NbrInteractionFactors] = Factor * Coefficients(m4, m3);
	this->M1Value[this->NbrInteractionFactors] = m3 - 1;
	this->M2Value[this->NbrInteractionFactors] = m4 + 1;
	this->M3Value[this->NbrInteractionFactors] = m3;

// 	this->NbrM3Values[this->NbrInteractionFactors] = 1;
// 	this->M3Values[this->NbrInteractionFactors] = new int [1];
// 	this->M1Value[this->NbrInteractionFactors] = m3;
// 	this->M2Value[this->NbrInteractionFactors] = m4;
// 	this->M3Values[this->NbrInteractionFactors][0] = m3 - 1;
	++this->NbrInteractionFactors;
      }
  this->NbrOneBodyInteractionFactors = this->LzMax + 1;
  this->OneBodyMValues = new int[this->NbrOneBodyInteractionFactors];
  this->OneBodyNValues = new int[this->NbrOneBodyInteractionFactors];
  this->OneBodyInteractionFactors = new double[this->NbrOneBodyInteractionFactors];
  this->OneBodyMValues[0] = 0;
  this->OneBodyNValues[0] = 0;
  this->OneBodyInteractionFactors[0] = - this->L2Factor * Coefficients(0, 1);
  for (int i = 1; i < this->LzMax; ++i)
    {
      this->OneBodyMValues[i] = i;
      this->OneBodyNValues[i] = i;
      this->OneBodyInteractionFactors[i] = -this->L2Factor * (Coefficients(i, i + 1) + Coefficients(i - 1, i));
    }	  
  this->OneBodyMValues[this->LzMax] = this->LzMax;
  this->OneBodyNValues[this->LzMax] = this->LzMax;
  this->OneBodyInteractionFactors[this->LzMax] = - this->L2Factor * Coefficients(this->LzMax - 1, this->LzMax);
  cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

