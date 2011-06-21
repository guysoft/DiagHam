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


#ifndef PARTICLEONSPHERETWOLANDAULEVELL2HAMILTONIAN_H
#define PARTICLEONSPHERETWOLANDAULEVELL2HAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSpinFullHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereTwoLandauLevelDeltaHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnSphereTwoLandauLevelL2Hamiltonian : public AbstractQHEOnSphereWithSpinFullHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:

  // twice the projected momentum total value
  int TotalLz;

  // multiplicative factor in front of the L^2 operator in the Hamiltonian
  double L2Factor;
  
  // twice the maximum momentum a particle can reach in the upper Landau level
  int LzMaxUp;
  
  // twice the maximum momentum a particle can reach in the lower Landau level
  int LzMaxDown;
      
  // shift to compensate for shift on down level of fermion 2LL class.
  int LzFermionDownShift;
  
  // shift to compensate for shift on up level of fermion 2LL class.
  int LzFermionUpShift;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state in the lowest Landau level
  // totalLz = twice the projected momentum total value
  // architecture = architecture to use for precalculation
  // l2Factor = multiplicative factor in front of the L^2 operator in the Hamiltonian
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereTwoLandauLevelL2Hamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, int totalLz,
					      AbstractArchitecture* architecture, double l2Factor = 1.0, long memory = -1, 
					      bool onDiskCacheFlag = false, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereTwoLandauLevelL2Hamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  List<Matrix*> RightInteractionOperators();  


 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

 // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination);
  
  // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);
  
  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, int* indexArray, double* coefficientArray, long& position);
};

// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnSphereTwoLandauLevelL2Hamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;;
  double* TmpInteractionFactor;
  int Index;

  
    
  // Annihilation operators acting on first LL (UpUp)
  for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
    {
      for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
        {
          Coefficient3 = particles->AuAu(index, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
          if ( Coefficient3 != 0.0 )
            {
              Coefficient3 *= vSource[index];                            
              // first UpUpUpUp
              TmpInteractionFactor = this->InteractionFactorsUpUpUpUp[j] + (i * this->NbrUpUpSectorIndicesPerSum[j]);	
              for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
                {
                  Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		  if (Index < Dim)
		      vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
                  ++TmpInteractionFactor;
		}
	    }
	}
    }
	
  // Annihilation operators acting on both levels (UpDown)
  for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
    {
      for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
        {
          Coefficient3 = particles->AuAd(index, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
          if ( Coefficient3 != 0.0 )
            {
              Coefficient3 *= vSource[index];                                       			      
	      // now UpDownUpDown
	      TmpInteractionFactor = this->InteractionFactorsUpDownUpDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j]);	
	      for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
		{
		  Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;			
		    }
		  ++TmpInteractionFactor;
		}	  
	    }
	}
    }
    
	
	
  // Annihilation operators acting on LLL (DownDown)
  for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
    {
      for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
        {
          Coefficient3 = particles->AdAd(index, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
          if ( Coefficient3 != 0.0 )
            {
              Coefficient3 *= vSource[index];
	      // now DownDownDownDown
	      TmpInteractionFactor = this->InteractionFactorsDownDownDownDown[j] + (i * this->NbrDownDownSectorIndicesPerSum[j]);	
	      for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
		{
		  Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		  if (Index < Dim)
		    {				      
		      vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;		      
		      /*if (index == 3235 && Index == 3235 ) 
		        {
			  cout << "Coeff : " << Coefficient << ", Coeff3 : " << Coefficient3 << ", Factor: " << (*TmpInteractionFactor) << ", Dest :" << vDestination[Index] << ", Source: " << vSource[index] <<  endl;
			}*/
		    }
		  ++TmpInteractionFactor;
		} 
	    }
	}
    }

}


// core part of the PartialFastMultiplicationMemory method involving two-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnSphereTwoLandauLevelL2Hamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient3 = 0.0;;
  int Dim = particles->GetHilbertSpaceDimension();

  if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic ) 
    {  
      for (int idx = firstComponent; idx < lastComponent; ++idx)
	{
	  // Annihilation operators acting on first LL (UpUp)
	  for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AuAu(idx, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {			
		      // first UpUpUpUp
		      for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
			{
			  Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}			  		     
		    }
		}
	    }
		
	  // Annihilation operators acting on both levels (UpDown)
	  for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AuAd(idx, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {				    		      				      
			// now UpDownUpDown
			for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
			  {
			    Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			    {
				++memory;
				++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			  }	  
		    }
		}
	    }
	    
		
		
	  // Annihilation operators acting on LLL (DownDown)
	  for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AdAd(idx, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {		      
		      // now DownDownDownDown
		      for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
			{
			  Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}	  
		    }
		}
	    }
	} 
    }
    
  else if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic ) 
    {
      for (int idx = firstComponent; idx < lastComponent; ++idx)
	{
	  // Annihilation operators acting on first LL (UpUp)
	  for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AuAu(idx, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {			
		      // first UpUpUpUp
		      for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
			{
			  Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}			  		      
		    }
		}
	    }
		
	  // Annihilation operators acting on both levels (UpDown)
	  for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AuAd(idx, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {				    		      				      
			// now UpDownUpDown
			for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
			  {
			    Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			    {
				++memory;
				++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			  }	  
		    }
		}
	    }
	    		
		
	  // Annihilation operators acting on LLL (DownDown)
	  for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AdAd(idx, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {		        
		      // now DownDownDownDown
		      for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
			{
			  Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}	  
		    }
		}
	    }
	}            
    }
	   
  if ((this->OneBodyInteractionFactorsUpDown != 0) || (this->OneBodyInteractionFactorsUpUp != 0))
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  ++memory;
	  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
	}
    }
  if (this->OneBodyInteractionFactorsUpDown != 0)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactorsUpDown; ++j)
	    {
	      Index = particles->AduAd(i, this->OneBodyMValuesUpDown[j], Coefficient);
	      if (Index < Dim)
		{
		  ++memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	    }
	}
    }
  if (this->OneBodyInteractionFactorsDownUp != 0)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactorsDownUp; ++j)
	    {
	      Index = particles->AddAu(i, this->OneBodyMValuesDownUp[j], Coefficient);
	      if (Index < Dim)
		{
		  ++memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	    }
	}
    }
}


// core part of the FastMultiplication method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnSphereTwoLandauLevelL2Hamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, int* indexArray, double* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient3 = 0.0;
  double* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();
 
  
  if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic ) 
    {
      // Annihilation operators acting on first LL (UpUp)
      for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAu(index, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{              
		  // first UpUpUpUp
		  TmpInteractionFactor = this->InteractionFactorsUpUpUpUp[j] + (i * this->NbrUpUpSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;		     
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }		      		  
		}
	    }
	}
	    
      // Annihilation operators acting on both levels (UpDown)
      for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAd(index, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{		  
		    // now UpDownUpDown
		    TmpInteractionFactor = this->InteractionFactorsUpDownUpDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j]);	
		    for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
		      {
			Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			if (Index < Dim)
			  {
			    indexArray[position] = Index;
			    coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			    ++position;
			  }
			++TmpInteractionFactor;
		      }	  
		  }
	    }
	}
	
	    
	    
      // Annihilation operators acting on LLL (DownDown)
      for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AdAd(index, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		
		  // now DownDownDownDown
		  TmpInteractionFactor = this->InteractionFactorsDownDownDownDown[j] + (i * this->NbrDownDownSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }	  				    		
		}
	    }
	}
    }
  else if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic ) 
    {
      // Annihilation operators acting on first LL (UpUp)
      for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAu(index, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{              
		  // first UpUpUpUp
		  TmpInteractionFactor = this->InteractionFactorsUpUpUpUp[j] + (i * this->NbrUpUpSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;		     
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		}
	    }
	}
	    
      // Annihilation operators acting on both levels (UpDown)
      for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAd(index, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		    // now UpDownUpDown
		    TmpInteractionFactor = this->InteractionFactorsUpDownUpDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j]);	
		    for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
		      {
			Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			if (Index < Dim)
			  {
			    indexArray[position] = Index;
			    coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			    ++position;
			  }
			++TmpInteractionFactor;
		      }	  
		  }
	    }
	}
	
	    
	    
      // Annihilation operators acting on LLL (DownDown)
      for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AdAd(index, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		  // now DownDownDownDown
		  TmpInteractionFactor = this->InteractionFactorsDownDownDownDown[j] + (i * this->NbrDownDownSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }	  
		}
	    }
	}      
    }
}


#endif
