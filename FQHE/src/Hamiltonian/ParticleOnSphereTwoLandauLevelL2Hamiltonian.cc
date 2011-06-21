////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of hamiltonian associated to particles on a sphere           //
//       with two Landau levels and where the hamiltonian is reduced          //
//                  to a simple total square angular momentum                 //
//                                                                            //
//                        last modification : 07/04/2010                      //
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


#include "Hamiltonian/ParticleOnSphereTwoLandauLevelL2Hamiltonian.h"
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

ParticleOnSphereTwoLandauLevelL2Hamiltonian::ParticleOnSphereTwoLandauLevelL2Hamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, int totalLz,
									     AbstractArchitecture* architecture, double l2Factor, long memory, bool onDiskCacheFlag,
									     char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax + 2;
  this->LzMaxUp = this->LzMax;
  this->LzMaxDown = this->LzMax - 1;
  this->LzFermionDownShift = 0;
  this->LzFermionUpShift = 0;   
  if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic ) 
    {
      this->LzFermionDownShift = 1;      
    }
  
  this->TotalLz = totalLz;
  this->NbrLzValue = this->LzMax + 1 ;
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
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;    
  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;
  this->NbrIntraSectorSums = 0;
  this->NbrInterSectorSums = 0;
  this->M1IntraValue = 0;
  this->M1InterValue = 0;
  
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
}

// destructor
//

ParticleOnSphereTwoLandauLevelL2Hamiltonian::~ParticleOnSphereTwoLandauLevelL2Hamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereTwoLandauLevelL2Hamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphereWithSpin*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereTwoLandauLevelL2Hamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereTwoLandauLevelL2Hamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereTwoLandauLevelL2Hamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift + 0.25 * this->L2Factor * ((double) (this->TotalLz * this->TotalLz));
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereTwoLandauLevelL2Hamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
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

Complex ParticleOnSphereTwoLandauLevelL2Hamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ParticleOnSphereTwoLandauLevelL2Hamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ParticleOnSphereTwoLandauLevelL2Hamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   

void ParticleOnSphereTwoLandauLevelL2Hamiltonian::EvaluateInteractionFactors()
{
  //set these to 0 so old method not used.
  this->NbrIntraSectorSums = 0;
  this->NbrInterSectorSums = 0;
  this->M1IntraValue = 0;
  this->M2InterValue = 0;

  //if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic ) 
    //{
      // Set the number of possible sums for each sector. 
      this->NbrUpUpSectorSums = 2 * this->LzMaxUp + 1; // Ranges from 0 to 2*LzMax for given sector.
      this->NbrDownDownSectorSums = 2 * this->LzMaxDown + 1 - 2; // the -2 is because a sum of 0 or 1 is not possible in the LLL.       
      this->NbrUpDownSectorSums = this->LzMaxUp + this->LzMaxDown; // goes from 1 to LzMaxUp + LzMaxDown
      
      //Allocate space for sums and set counts to zero.
      this->NbrUpUpSectorIndicesPerSum = new int[this->NbrUpUpSectorSums];
      this->NbrUpDownSectorIndicesPerSum = new int [this->NbrUpDownSectorSums];
      this->NbrDownDownSectorIndicesPerSum = new int [this->NbrDownDownSectorSums];
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	this->NbrUpUpSectorIndicesPerSum[i] = 0;
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
	this->NbrUpDownSectorIndicesPerSum[i] = 0;
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	this->NbrDownDownSectorIndicesPerSum[i] = 0;

      // Count number of combinations that sum to each.
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 0; m2 <= this->LzMaxUp; ++m2)
	  this->NbrUpUpSectorIndicesPerSum[m1 + m2]++;   
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  this->NbrUpDownSectorIndicesPerSum[m1 + m2 - 1]++;
      for (int m1 = 1; m1 <= this->LzMaxDown; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 2]++;

      // Allocate sapce for indices and reset counters to 0 so can be used as indices.
      this->UpUpSectorIndicesPerSum = new int* [this->NbrUpUpSectorSums];
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	if (this->NbrUpUpSectorIndicesPerSum[i] > 0)
	{
	  this->UpUpSectorIndicesPerSum[i] = new int[2 * this->NbrUpUpSectorIndicesPerSum[i]];      
	  this->NbrUpUpSectorIndicesPerSum[i] = 0;
	}
      this->UpDownSectorIndicesPerSum = new int* [this->NbrUpDownSectorSums];
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
      if (this->NbrUpDownSectorIndicesPerSum[i] > 0)
      {
	  this->UpDownSectorIndicesPerSum[i] = new int[2 * this->NbrUpDownSectorIndicesPerSum[i]];      
	  this->NbrUpDownSectorIndicesPerSum[i] = 0;
      }
      this->DownDownSectorIndicesPerSum = new int* [this->NbrDownDownSectorSums];
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	if (this->NbrDownDownSectorIndicesPerSum[i] > 0)
	{
	  this->DownDownSectorIndicesPerSum[i] = new int[2 * this->NbrDownDownSectorIndicesPerSum[i]];      
	  this->NbrDownDownSectorIndicesPerSum[i] = 0;
	}

      // set the indices.
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 0; m2 <= this->LzMaxUp; ++m2)
	  {
	    this->UpUpSectorIndicesPerSum[(m1 + m2)][this->NbrUpUpSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
	    this->UpUpSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrUpUpSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
	    ++this->NbrUpUpSectorIndicesPerSum[(m1 + m2)];
	  }
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  {
	    this->UpDownSectorIndicesPerSum[m1 + m2 - 1][this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)] << 1] = m1;
	    this->UpDownSectorIndicesPerSum[m1 + m2 - 1][1 + (this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)] << 1)] = m2;
	    ++this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)];
	  }
      for (int m1 = 1; m1 <= this->LzMaxDown; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  {
	    this->DownDownSectorIndicesPerSum[m1 + m2 - 2][this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 2] << 1] = m1;
	    this->DownDownSectorIndicesPerSum[m1 + m2 - 2][1 + (this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 2] << 1)] = m2;
	    ++this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 2];
	  }		      
    /*} 
  else if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic ) 
    {
      // Set the number of possible sums for each sector. 
      this->NbrUpUpSectorSums     = 2 * this->LzMaxUp - 1; // Ranges from 1 to 2*LzMaxUp-1 for given sector.
      this->NbrDownDownSectorSums = 2 * this->LzMaxDown - 1 - 2; // (Ranges from 3 to LzMaxDown + LzMaxDown - 1)  the -2 is because a sum of 0 or 1 is not possible in the LLL. 
      this->NbrUpDownSectorSums   = this->LzMaxUp + this->LzMaxDown ; // goes from 1 to LzMaxUp + LzMaxDown
      
      //Allocate space for sums and set counts to zero.
      this->NbrUpUpSectorIndicesPerSum = new int[this->NbrUpUpSectorSums];
      this->NbrUpDownSectorIndicesPerSum = new int [this->NbrUpDownSectorSums];
      this->NbrDownDownSectorIndicesPerSum = new int [this->NbrDownDownSectorSums];
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	this->NbrUpUpSectorIndicesPerSum[i] = 0;
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
	this->NbrUpDownSectorIndicesPerSum[i] = 0;
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	this->NbrDownDownSectorIndicesPerSum[i] = 0;

      // Count number of combinations that sum to each.
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 0; m2 <= this->LzMaxUp; ++m2)
	  if ( m1 != m2 ) this->NbrUpUpSectorIndicesPerSum[m1 + m2 - 1]++;   
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  this->NbrUpDownSectorIndicesPerSum[m1 + m2 - 1]++;
      for (int m1 = 1; m1 <= this->LzMaxDown; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  if ( m1 != m2 ) this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 3]++;

      // Allocate sapce for indices and reset counters to 0 so can be used as indices.
      this->UpUpSectorIndicesPerSum = new int* [this->NbrUpUpSectorSums];
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	if (this->NbrUpUpSectorIndicesPerSum[i] > 0)
	  {
	    this->UpUpSectorIndicesPerSum[i] = new int[2 * this->NbrUpUpSectorIndicesPerSum[i]];      
	    this->NbrUpUpSectorIndicesPerSum[i] = 0;
	  } 
	else 
	  {
	    this->UpUpSectorIndicesPerSum[i] = 0;
	  }
      this->UpDownSectorIndicesPerSum = new int* [this->NbrUpDownSectorSums];
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
        if (this->NbrUpDownSectorIndicesPerSum[i] > 0)
          {
	    this->UpDownSectorIndicesPerSum[i] = new int[2 * this->NbrUpDownSectorIndicesPerSum[i]];      
	    this->NbrUpDownSectorIndicesPerSum[i] = 0;
          }
	else 
	  {
	    this->UpDownSectorIndicesPerSum[i] = 0;
	  }
      this->DownDownSectorIndicesPerSum = new int* [this->NbrDownDownSectorSums];
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	if (this->NbrDownDownSectorIndicesPerSum[i] > 0)
	  {
	    this->DownDownSectorIndicesPerSum[i] = new int[2 * this->NbrDownDownSectorIndicesPerSum[i]];      
	    this->NbrDownDownSectorIndicesPerSum[i] = 0;
	  }
	else 
	  {
	    this->DownDownSectorIndicesPerSum[i] = 0;
	  }  
	

      // set the indices.
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 0; m2 <= this->LzMaxUp; ++m2)
	  {
	    if ( m1 != m2 ) 
	      {
		this->UpUpSectorIndicesPerSum[m1 + m2 - 1][this->NbrUpUpSectorIndicesPerSum[m1 + m2 - 1] << 1] = m1;
		this->UpUpSectorIndicesPerSum[m1 + m2 - 1][(this->NbrUpUpSectorIndicesPerSum[m1 + m2 - 1] << 1) + 1] = m2;
		++this->NbrUpUpSectorIndicesPerSum[m1 + m2 - 1];
	      }
	  }
	  
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  {
	    this->UpDownSectorIndicesPerSum[m1 + m2 - 1][this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)] << 1] = m1;
	    this->UpDownSectorIndicesPerSum[m1 + m2 - 1][(this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)] << 1) + 1] = m2 - this->LzFermionDownShift;
	    ++this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)];
	  }
	  
      for (int m1 = 1; m1 <= this->LzMaxDown; ++m1)
	{
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  {
	    if ( m1 != m2 ) 
	      {
		this->DownDownSectorIndicesPerSum[m1 + m2 - 3][this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 3] << 1] = m1 - this->LzFermionDownShift;
		this->DownDownSectorIndicesPerSum[m1 + m2 - 3][(this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 3] << 1) + 1] = m2 - this->LzFermionDownShift;
		++this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 3];
	      }
	  }
	}	
	
    }*/
    
  // Create interaction factor arrays and initialise to 0.
  
  // the three that end in UpUp
  this->InteractionFactorsUpUpUpUp = new double* [this->NbrUpUpSectorSums];
  this->InteractionFactorsUpDownUpUp = new double* [this->NbrUpUpSectorSums];
  this->InteractionFactorsDownDownUpUp = new double* [this->NbrUpUpSectorSums];
  //this is prob unnecessary but no harm.
  for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
    {
      this->InteractionFactorsUpUpUpUp[i] = 0;
      this->InteractionFactorsUpDownUpUp[i] = 0;
      this->InteractionFactorsDownDownUpUp[i] = 0;
    }
  
  // finally the three that end in Up Down
  this->InteractionFactorsUpUpUpDown = new double* [this->NbrUpDownSectorSums];
  this->InteractionFactorsUpDownUpDown = new double* [this->NbrUpDownSectorSums];
  this->InteractionFactorsDownDownUpDown = new double* [this->NbrUpDownSectorSums];
  //this is prob unnecessary but no harm.
  for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
    {
      this->InteractionFactorsUpUpUpDown[i] = 0;
      this->InteractionFactorsUpDownUpDown[i] = 0;
      this->InteractionFactorsDownDownUpDown[i] = 0;
    }
  
//now the trhee that end in Down Down
  this->InteractionFactorsDownDownDownDown = new double* [this->NbrDownDownSectorSums];
  this->InteractionFactorsUpDownDownDown = new double* [this->NbrDownDownSectorSums];
  this->InteractionFactorsUpUpDownDown = new double* [this->NbrDownDownSectorSums];
  //this is prob unnecessary but no harm.
  for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
    {
      this->InteractionFactorsUpUpDownDown[i] = 0;
      this->InteractionFactorsUpDownDownDown[i] = 0;
      this->InteractionFactorsDownDownDownDown[i] = 0;
    }
        
  this->NbrOneBodyInteractionFactorsDownDown = this->LzMaxDown;        
  this->OneBodyInteractionFactorsDownDown = new double [this->NbrOneBodyInteractionFactorsDownDown]; 
  this->OneBodyMValuesDownDown = new int[this->NbrOneBodyInteractionFactorsDownDown];
  for (int i = 0; i < this->LzMaxDown; ++i)
    {
      this->OneBodyMValuesDownDown[i] = i+1;
      this->OneBodyInteractionFactorsDownDown[i] = 0.0;
    }
    
  this->NbrOneBodyInteractionFactorsUpUp = this->LzMaxUp + 1;        
  this->OneBodyInteractionFactorsUpUp = new double [this->NbrOneBodyInteractionFactorsUpUp]; 
  this->OneBodyMValuesUpUp = new int[this->NbrOneBodyInteractionFactorsUpUp];
  for (int i = 0; i <= this->LzMaxUp; ++i)
    {
      this->OneBodyMValuesUpUp[i] = i;
      this->OneBodyInteractionFactorsUpUp[i] = 0.0;
    }
  
  long TotalNbrInteractionFactors = 0;
  
  //when not using pseudopotentials will use expression that was worked out for the delta interaction.
  
  //will start with the interactions confined completely to the LLL (Down).
  //cout << "Down Down Down Down" << endl;
  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{
	  this->InteractionFactorsDownDownDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrDownDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
	      int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
	      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  if (  m3 == (m1 + 2) && m4 == (m2 - 2))
		    {
		      double Q = (double)(this->LzMaxDown-1)/2.0;
		      double m1r = (double)m1/2.0; 
		      double m2r = (double)m2/2.0;
		      if ( m2 != m3 ) 
			{
			  this->InteractionFactorsDownDownDownDown[i][Index] = 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r + 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r - 1.0));
			}
		      else if ( m2 == m3 )
			{
			  this->InteractionFactorsDownDownDownDown[i][Index] = 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r + 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r - 1.0));
			  this->OneBodyInteractionFactorsDownDown[((m1 + this->LzMaxDown + 1) >> 1)-1] += 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r + 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r - 1.0));
			}		      
		    }
		  else if ( m3 == (m1 - 2) && m4 == (m2 + 2))
		    {
		      double Q = (double)(this->LzMaxDown-1)/2.0;
		      double m1r = (double)m1/2.0; 
		      double m2r = (double)m2/2.0;
		      if ( m2 != m3 ) 
			{
			  this->InteractionFactorsDownDownDownDown[i][Index] = 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r - 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r + 1.0));
			}
		      else if ( m2 == m3 ) 
			{
			  this->InteractionFactorsDownDownDownDown[i][Index] = 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r - 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r + 1.0));
			  this->OneBodyInteractionFactorsDownDown[((m1 + this->LzMaxDown + 1) >> 1)-1] += 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r - 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r + 1.0));
			}		      
		    }
		  else 
		    {
		      this->InteractionFactorsDownDownDownDown[i][Index] = 0.0;
		    }
		  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownDownDown[i][Index] << endl;
		  ++Index;
		  TotalNbrInteractionFactors++;
		}
	    }
	}
    }
    
  cout << "DownDownDownDown Terms" << endl;
  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{	      
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1;
	      int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1;
	      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxDown - 1;
		  int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1;
		  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsDownDownDownDown[i][Index] << endl;
		  ++Index;
		}
	    }
	}
    }
    
  for ( int i = 0 ; i < this->LzMaxDown ; i++ ) 
    {
      cout << "One body term on " << this->OneBodyMValuesDownDown[i] << ": " << this->OneBodyInteractionFactorsDownDown[i] << endl ;	  
    }
    
    /*//cout << "Up Down Down Down" << endl;
    //now we set the interaction terms where a single operator acts on the first LL and the rest on the LLL 
    for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{
	  this->InteractionFactorsUpDownDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i+1]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
	      int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
	      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i+1]; ++j2)
		{
		  int m3 = (this->UpDownSectorIndicesPerSum[i+1][j2 << 1] << 1) - this->LzMaxUp;
		  int m4 = (this->UpDownSectorIndicesPerSum[i+1][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  
		  this->InteractionFactorsUpDownDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,0.0 ,(double)m1 /2.0 ,0.0 ,(double)m2 /2.0,1.0,(double)m3/2.0,0.0,(double)m4/2.0);
		  //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpDownDownDown[i][Index] << endl;
		  ++Index;
		  TotalNbrInteractionFactors++;
		}
	    }
	}
    }
    
    //cout << "Up Up Down Down" << endl;
    for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{
	  this->InteractionFactorsUpUpDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+2]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
	      int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
	      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j2)
		{
		  int m3 = (this->UpUpSectorIndicesPerSum[i+2][j2 << 1] << 1) - this->LzMaxUp;
		  int m4 = (this->UpUpSectorIndicesPerSum[i+2][(j2 << 1) + 1] << 1) - this->LzMaxUp;
		  
		  this->InteractionFactorsUpUpDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/ 2.0 ,0.0,(double)m1/2.0 ,0.0,(double)m2 /2.0,1.0,(double)m3 /2.0 ,1.0,(double)m4/2.0);
		  //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpUpDownDown[i][Index] << endl;
		  ++Index;
		  TotalNbrInteractionFactors++;
		}
	    }
	}
    }
    
    //cout << "Down Down Up Down" << endl;
    for (int i = 1; i < this->NbrUpDownSectorSums-1; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{
	  this->InteractionFactorsDownDownUpDown[i-1] = new double[this->NbrDownDownSectorIndicesPerSum[i-1] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
	      int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
	      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i-1]; ++j2)
		{
		  int m3 = (this->DownDownSectorIndicesPerSum[i-1][j2 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  int m4 = (this->DownDownSectorIndicesPerSum[i-1][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  
		  this->InteractionFactorsDownDownUpDown[i-1][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1/2.0 ,0.0,(double)m2/2.0 ,0.0,(double)m3/2.0 ,0.0,(double)m4/2.0);
		  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownUpDown[i-1][Index] << endl;
		  ++Index;
		  TotalNbrInteractionFactors++;
		}
	    }
	}
    }*/
    
    //cout << "Up Down Up Down" << endl;
    for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{
	  this->InteractionFactorsUpDownUpDown[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
	      int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
	      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
		  int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  
		  if (  m3 == (m1 + 2) && m4 == (m2 - 2))
		    {
		      double Qu = (double)(this->LzMaxUp)/2.0;
		      double Qd = (double)(this->LzMaxDown-1)/2.0;
		      double m1r = (double)m1/2.0; 
		      double m2r = (double)m2/2.0;		      
		      this->InteractionFactorsUpDownUpDown[i][Index] =  sqrt(Qu * (Qu + 1.0) - m1r*(m1r + 1.0)) * sqrt(Qd * (Qd + 1.0) - m2r*(m2r - 1.0));			
		    }
		  else if ( m3 == (m1 - 2) && m4 == (m2 + 2))
		    {
		      double Qu = (double)(this->LzMaxUp)/2.0;
		      double Qd = (double)(this->LzMaxDown-1)/2.0;
		      double m1r = (double)m1/2.0; 
		      double m2r = (double)m2/2.0;
		      this->InteractionFactorsUpDownUpDown[i][Index] =  sqrt(Qu * (Qu + 1.0) - m1r*(m1r - 1.0)) * sqrt(Qd * (Qd + 1.0) - m2r*(m2r + 1.0));
		    }
		  else 
		    {
		      this->InteractionFactorsUpDownUpDown[i][Index] = 0.0;
		    }		  		  
		  
		  ++Index;
		  TotalNbrInteractionFactors++;
		}
	    }
	}
    }

  cout << "UpDownUpDown Terms" << endl;
  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{	      
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp;
	      int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1;
	      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
		  int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1;
		  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpDownUpDown[i][Index] << endl;
		  ++Index;
		}
	    }
	}
    }

    //cout << "Up Up Up Down" << endl;
    /*for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{
	  this->InteractionFactorsUpUpUpDown[i] = new double[this->NbrUpUpSectorIndicesPerSum[i+1] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
	      int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
	      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i+1]; ++j2)
		{
		  int m3 = (this->UpUpSectorIndicesPerSum[i+1][j2 << 1] << 1) - this->LzMaxUp;
		  int m4 = (this->UpUpSectorIndicesPerSum[i+1][(j2 << 1) + 1] << 1) - this->LzMaxUp;
		  
		  this->InteractionFactorsUpUpUpDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1 /2.0 ,0.0,(double)m2 /2.0,1.0,(double)m3 /2.0,1.0,(double)m4/2.0);
		  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpUpUpDown[i][Index] << endl;
		  ++Index;
		  TotalNbrInteractionFactors++;
		}
	    }
	}
    }
    
    //cout << "Down Down Up Up" << endl;
    for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{
	  this->InteractionFactorsDownDownUpUp[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+2]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j1)
	    {
	      int m1 = (this->UpUpSectorIndicesPerSum[i+2][j1 << 1] << 1) - this->LzMaxUp ;
	      int m2 = (this->UpUpSectorIndicesPerSum[i+2][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
	      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  
		  this->InteractionFactorsDownDownUpUp[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1 /2.0 ,1.0,(double)m2 /2.0,0.0,(double)m3 /2.0,0.0,(double)m4 /2.0);
		  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownUpUp[i][Index] << endl;
		  ++Index;
		  TotalNbrInteractionFactors++;
		}
	    }
	}
    }
    
    //cout << "Up Down Up Up" << endl;
    for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{
	  this->InteractionFactorsUpDownUpUp[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+1]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i+1]; ++j1)
	    {
	      int m1 = (this->UpUpSectorIndicesPerSum[i+1][j1 << 1] << 1) - this->LzMaxUp ;
	      int m2 = (this->UpUpSectorIndicesPerSum[i+1][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
	      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
		  int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  
		  this->InteractionFactorsUpDownUpUp[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/2.0 ,1.0,(double)m1 /2.0 ,1.0,(double)m2 /2.0,1.0,(double)m3 /2.0,0.0,(double)m4 /2.0);
		  //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpDownUpUp[i][Index] << endl;
		  ++Index;
		  TotalNbrInteractionFactors++;
		}
	    }
	}
    }*/
      
    //cout << "Up Up Up Up" << endl;    
    for (int i = 0; i < this->NbrUpUpSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrUpUpSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{
	  this->InteractionFactorsUpUpUpUp[i] = new double[this->NbrUpUpSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
	      int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
	      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
		  int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxUp ;
		  
		  if (  m3 == (m1 + 2) && m4 == (m2 - 2))
		    {
		      double Q = (double)(this->LzMaxUp)/2.0;
		      double m1r = (double)m1/2.0; 
		      double m2r = (double)m2/2.0;
		      if ( m2 != m3 ) 
			{
			  this->InteractionFactorsUpUpUpUp[i][Index] = 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r + 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r - 1.0));
			}
		      else if ( m2 == m3 )
			{
			  this->InteractionFactorsUpUpUpUp[i][Index] = 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r + 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r - 1.0));
			  this->OneBodyInteractionFactorsUpUp[((m1 + this->LzMaxUp) >> 1)] += 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r + 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r - 1.0));
			}		      
		    }
		  else if ( m3 == (m1 - 2) && m4 == (m2 + 2))
		    {
		      double Q = (double)(this->LzMaxUp)/2.0;
		      double m1r = (double)m1/2.0; 
		      double m2r = (double)m2/2.0;
		      if ( m2 != m3 ) 
			{
			  this->InteractionFactorsUpUpUpUp[i][Index] = 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r - 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r + 1.0));
			}
		      else if ( m2 == m3 ) 
			{
			  this->InteractionFactorsUpUpUpUp[i][Index] = 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r - 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r + 1.0));
			  this->OneBodyInteractionFactorsUpUp[((m1 + this->LzMaxUp) >> 1)] += 0.5 * sqrt(Q * (Q + 1.0) - m1r*(m1r - 1.0)) * sqrt(Q * (Q + 1.0) - m2r*(m2r + 1.0));
			}		      
		    }
		  else 
		    {
		      this->InteractionFactorsUpUpUpUp[i][Index] = 0.0;
		    }		  		  		  
		  ++Index;
		  TotalNbrInteractionFactors++;
		}
	    }
	}
    }	 
 
  cout << "UpUpUpUp Terms" << endl;
  for (int i = 0; i < this->NbrUpUpSectorSums; ++i) // go through the possible sums of Lz values on LLL
    {
      if (this->NbrUpUpSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	{	      
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp;
	      int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxUp;
	      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
		  int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxUp;
		  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpUpUpUp[i][Index] << endl;
		  ++Index;
		}
	    }
	}
    }
    
  for ( int i = 0 ; i <= this->LzMaxUp ; i++ ) 
    {
      cout << "One body term on " << this->OneBodyMValuesUpUp[i] << ": " << this->OneBodyInteractionFactorsUpUp[i] << endl ;	  
    }
 
  /*RealMatrix CoefficientsDownDown (this->LzMax + 1, this->LzMax + 1);
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax + 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	CoefficientsDownDown(i, j) = TmpCoefficient;
    }
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax - 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	CoefficientsDownDown(j, i) *= 0.5 * TmpCoefficient;
    }

  int LzMaxUp = this->LzMax + 2;
  RealMatrix CoefficientsUpUp (LzMaxUp + 1, LzMaxUp + 1);
  for (int i = 0; i <= LzMaxUp; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((LzMaxUp + 2) * LzMaxUp) - (((2 * i) - LzMaxUp) * ((2 * i) - LzMaxUp + 2))))));
      for (int j = 0; j <= LzMaxUp; ++j)
	CoefficientsUpUp(i, j) = TmpCoefficient;
    }
  for (int i = 0; i <= LzMaxUp; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((LzMaxUp + 2) * LzMaxUp) - (((2 * i) - LzMaxUp) * ((2 * i) - LzMaxUp - 2))))));
      for (int j = 0; j <= LzMaxUp; ++j)
	CoefficientsUpUp(j, i) *= 0.5 * TmpCoefficient;
    }

  RealMatrix CoefficientsUpDown (LzMaxUp + 1, this->LzMax + 1);
  for (int i = 0; i <= LzMaxUp; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((LzMaxUp + 2) * LzMaxUp) - (((2 * i) - LzMaxUp) * ((2 * i) - LzMaxUp + 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	CoefficientsUpDown(i, j) = TmpCoefficient;
    }
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax - 2))))));
      for (int j = 0; j <= LzMaxUp; ++j)
	CoefficientsUpDown(j, i) *= 0.5 * TmpCoefficient;
    }

  //  this->L2Factor  = 1.0;

  double Factor = 2.0 * this->L2Factor;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphereWithSpin::FermionicStatistic)
    {
      Factor *= -1.0;
      this->NbrM12IntraIndices = this->LzMax * (this->LzMax - 1) + 1;
    }
  else this->NbrM12IntraIndices = this->LzMax * this->LzMax + 1;
  this->M1IntraValue = new int [this->NbrM12IntraIndices];
  this->M2IntraValue = new int [this->NbrM12IntraIndices];
  this->M3IntraValues = new int* [this->NbrM12IntraIndices];
  this->NbrM3IntraValues  = new int [this->NbrM12IntraIndices];
  this->M12InteractionFactorsupup = new double [this->NbrM12IntraIndices];
  this->M12InteractionFactorsdowndown = new double [this->NbrM12IntraIndices];
  this->NbrM12IntraIndices = 0;

  for (int m3 = 1; m3 <= this->LzMax; ++m3)
    {
      if ((this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic) ||
	  (m3 != 2))
	{
	  this->M12InteractionFactorsupup[this->NbrM12IntraIndices] = Factor * CoefficientsUpUp(0, m3);
	  this->M12InteractionFactorsdowndown[this->NbrM12IntraIndices] = Factor * CoefficientsDownDown(0, m3);
	  this->M1IntraValue[this->NbrM12IntraIndices] = m3;
	  this->M2IntraValue[this->NbrM12IntraIndices] = 0;
	  this->M3IntraValues[this->NbrM12IntraIndices] = new int [1];
	  this->NbrM3IntraValues[this->NbrM12IntraIndices] = 1;
	  this->M3IntraValues[this->NbrM12IntraIndices][0] =  m3 - 1;
	  ++this->NbrM12IntraIndices;
	}
    }
  for (int m4 = 1; m4 < this->LzMax; ++m4)
    {
      int m3= 1;
      for (; m3 < m4; ++m3)
	{
	  if ((this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic) ||
	      (m3 != (m4 + 2)))
	    {
	      this->M12InteractionFactorsupup[this->NbrM12IntraIndices] = Factor * CoefficientsUpUp(m4, m3);
	      this->M12InteractionFactorsdowndown[this->NbrM12IntraIndices] = Factor * CoefficientsDownDown(m4, m3);
	      this->M1IntraValue[this->NbrM12IntraIndices] = m3;
	      this->M2IntraValue[this->NbrM12IntraIndices] = m4;
	      this->M3IntraValues[this->NbrM12IntraIndices] = new int [1];
	      this->NbrM3IntraValues[this->NbrM12IntraIndices] = 1;
	      this->M3IntraValues[this->NbrM12IntraIndices][0] =  m3 - 1;
	      ++this->NbrM12IntraIndices;
	    }
	}
      if (this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic)
	{
	  this->M12InteractionFactorsupup[this->NbrM12IntraIndices] = Factor * CoefficientsUpUp(m4, m3);
	  this->M12InteractionFactorsdowndown[this->NbrM12IntraIndices] = Factor * CoefficientsDownDown(m4, m3);
	  this->M1IntraValue[this->NbrM12IntraIndices] = m3;
	  this->M2IntraValue[this->NbrM12IntraIndices] = m4;
	  this->M3IntraValues[this->NbrM12IntraIndices] = new int [1];
	  this->NbrM3IntraValues[this->NbrM12IntraIndices] = 1;
	  this->M3IntraValues[this->NbrM12IntraIndices][0] = m3 - 1;
	  ++this->NbrM12IntraIndices;	  
	}
      ++m3;
      for (; m3 <= this->LzMax; ++m3)
	{
	  if ((this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic) ||
	      (m3 != (m4 + 2)))
            {
	      this->M12InteractionFactorsupup[this->NbrM12IntraIndices] = Factor * CoefficientsUpUp(m4, m3);
	      this->M12InteractionFactorsdowndown[this->NbrM12IntraIndices] = Factor * CoefficientsDownDown(m4, m3);
	      this->M1IntraValue[this->NbrM12IntraIndices] = m3;
	      this->M2IntraValue[this->NbrM12IntraIndices] = m4;
	      this->M3IntraValues[this->NbrM12IntraIndices] = new int [1];
	      this->NbrM3IntraValues[this->NbrM12IntraIndices] = 1;
	      this->M3IntraValues[this->NbrM12IntraIndices][0] =  m3 - 1;
	      ++this->NbrM12IntraIndices;
	    }
	}
    }

  this->NbrM12InterIndices = (this->LzMax + 3) * this->LzMax;
  this->M1InterValue = new int [this->NbrM12InterIndices];
  this->M2InterValue = new int [this->NbrM12InterIndices];
  this->M3InterValues = new int* [this->NbrM12InterIndices];
  this->NbrM3InterValues  = new int [this->NbrM12InterIndices];
  this->M12InteractionFactorsupdown = new double [4 * this->NbrM12InterIndices];
  this->NbrM12InterIndices = 0;
  int TmpNbrM12InterIndices = 0;
  
  //   this->L2Factor  = 0.0;
  //   Factor = 2.0 * this->L2Factor;
  //   if (this->Particles->GetParticleStatistic() == ParticleOnSphereWithSpin::FermionicStatistic)
  //     Factor *= -1.0;
  
  for (int m3 = 1; m3 <= this->LzMax; ++m3)
    {
      this->M12InteractionFactorsupdown[TmpNbrM12InterIndices++] = Factor * CoefficientsUpDown(0, m3);
      this->M1InterValue[this->NbrM12InterIndices] = m3;
      this->M2InterValue[this->NbrM12InterIndices] = 0;
      this->M3InterValues[this->NbrM12InterIndices] = new int [1];
      this->NbrM3InterValues[this->NbrM12InterIndices] = 1;
      this->M3InterValues[this->NbrM12InterIndices][0] =  m3 - 1;
      ++this->NbrM12InterIndices;
    }
  for (int m4 = 1; m4 < this->LzMax; ++m4)
    {
      int m3= 0;
      this->M12InteractionFactorsupdown[TmpNbrM12InterIndices++] = Factor * CoefficientsUpDown(m3, m4);
      this->M1InterValue[this->NbrM12InterIndices] = m3;
      this->M2InterValue[this->NbrM12InterIndices] = m4;
      this->M3InterValues[this->NbrM12InterIndices] = new int [1];
      this->NbrM3InterValues[this->NbrM12InterIndices] = 1;
      this->M3InterValues[this->NbrM12InterIndices][0] =  m3 + 1;
      ++this->NbrM12InterIndices;
      ++m3;
      for (; m3 < this->LzMax; ++m3)
	{
	  this->M12InteractionFactorsupdown[TmpNbrM12InterIndices++] = Factor * CoefficientsUpDown(m3, m4);
	  this->M12InteractionFactorsupdown[TmpNbrM12InterIndices++] = Factor * CoefficientsUpDown(m4, m3);
	  this->M1InterValue[this->NbrM12InterIndices] = m3;
	  this->M2InterValue[this->NbrM12InterIndices] = m4;
	  this->M3InterValues[this->NbrM12InterIndices] = new int [2];
	  this->NbrM3InterValues[this->NbrM12InterIndices] = 2;
	  this->M3InterValues[this->NbrM12InterIndices][0] =  m3 + 1;
	  this->M3InterValues[this->NbrM12InterIndices][1] =  m3 - 1;
	  ++this->NbrM12InterIndices;
	}
      this->M12InteractionFactorsupdown[TmpNbrM12InterIndices++] = Factor * CoefficientsUpDown(m4, m3);
      this->M1InterValue[this->NbrM12InterIndices] = m3;
      this->M2InterValue[this->NbrM12InterIndices] = m4;
      this->M3InterValues[this->NbrM12InterIndices] = new int [1];
      this->NbrM3InterValues[this->NbrM12InterIndices] = 1;
      this->M3InterValues[this->NbrM12InterIndices][0] =  m3 - 1;
      ++this->NbrM12InterIndices;
    }
  for (int m3 = 0; m3 < this->LzMax; ++m3)
    {
      this->M12InteractionFactorsupdown[TmpNbrM12InterIndices++] = Factor * CoefficientsUpDown(m3, this->LzMax);
      this->M1InterValue[this->NbrM12InterIndices] = m3;
      this->M2InterValue[this->NbrM12InterIndices] = this->LzMax;
      this->M3InterValues[this->NbrM12InterIndices] = new int [1];
      this->NbrM3InterValues[this->NbrM12InterIndices] = 1;
      this->M3InterValues[this->NbrM12InterIndices][0] =  m3 + 1;
      ++this->NbrM12InterIndices;
    }

  //   this->L2Factor  = 0.0;

  this->OneBodyInteractionFactorsupup = new double[LzMaxUp + 1];
  this->OneBodyInteractionFactorsdowndown = new double[this->LzMax + 1];
  this->OneBodyInteractionFactorsupup[0] = this->L2Factor * CoefficientsUpUp(0, 1);
  this->OneBodyInteractionFactorsdowndown[0] = this->L2Factor * CoefficientsDownDown(0, 1);
  for (int i = 1; i < LzMaxUp; ++i)
    {
      this->OneBodyInteractionFactorsupup[i] = this->L2Factor * (CoefficientsUpUp(i, i + 1) + CoefficientsUpUp(i - 1, i));
    }
  for (int i = 1; i < this->LzMax; ++i)
    {
      this->OneBodyInteractionFactorsdowndown[i] = this->L2Factor * (CoefficientsDownDown(i, i + 1) + CoefficientsDownDown(i - 1, i));
    }	  
  this->OneBodyInteractionFactorsupup[LzMaxUp] = this->L2Factor * CoefficientsUpUp(this->LzMax - 1, this->LzMax);
  this->OneBodyInteractionFactorsdowndown[this->LzMax] = this->L2Factor * CoefficientsDownDown(this->LzMax - 1, this->LzMax);  
  this->OneBodyInteractionFactorsupdown = NULL;
  cout << "nbr interaction = " << ((2 * (this->NbrM12IntraIndices + this->LzMax)) + 2 + this->NbrM12InterIndices) << endl;
  cout << "====================================" << endl;*/
  //this->NbrOneBodyInteractionFactorsUpUp = 0;
  this->NbrOneBodyInteractionFactorsUpDown = 0;
  this->NbrOneBodyInteractionFactorsDownUp = 0;
    
}

