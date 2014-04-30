////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                to particles on a torus with magnetic translations          //
//                            and n-body interactions                         //
//                                                                            //
//                        last modification : 30/04/2014                      //
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
#include "Hamiltonian/AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian.h"

#include <iostream>


// destructor
//

AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian::~AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian()
{
  delete[] this->ExponentialFactors;
}

// evaluate all exponential factors
//   

void AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex[this->MaxMomentum];
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      this->ExponentialFactors[i] = Phase(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->MaxMomentum));
    }
}

// get all the indices that should appear in the annihilation/creation operators
//

void AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian::GetIndices()
{
  this->NbrSectorSums = this->NbrLzValue;
  this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
  for (int i = 0; i < this->NbrSectorSums; ++i)
    this->NbrSectorIndicesPerSum[i] = 0;      
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  if (this->NbrSectorIndicesPerSum[i]  > 0)
	    {
	      this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	      this->NbrSectorIndicesPerSum[i] = 0;
	    }
	}
      if (this->TwoBodyFlag == false)
	{
	  for (int i = 0; i < this->NbrSectorSums; ++i)
	    {
	      this->NbrSectorIndicesPerSum[i] = 0;
	      delete [] this->SectorIndicesPerSum[i];
	    }
	  this->InteractionFactors = new Complex* [this->NbrSectorSums];
	}
      else
	{
	  for (int m1 = 0; m1 < this->LzMax; ++m1)
	    for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	      {
		int TmpSum = (m1 + m2) % this->NbrLzValue;
		this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = m1;
		this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = m2;
		++this->NbrSectorIndicesPerSum[TmpSum];    
	      }
	}
    }
  else
    {
      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  if (this->NbrSectorIndicesPerSum[i]  > 0)
	    {
	      this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	      this->NbrSectorIndicesPerSum[i] = 0;
	    }
	}
      if (this->TwoBodyFlag == false)
	{
	  for (int i = 0; i < this->NbrSectorSums; ++i)
	    {
	      this->NbrSectorIndicesPerSum[i] = 0;
	      delete[] this->SectorIndicesPerSum[i];
	    }
	  this->InteractionFactors = new Complex* [this->NbrSectorSums];
	}
      else
	{
	  for (int m1 = 0; m1 <= this->LzMax; ++m1)
	    for (int m2 = m1; m2 <= this->LzMax; ++m2)
	      {
		int TmpSum = (m1 + m2) % this->NbrLzValue;
		this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = m1;
		this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = m2;
		++this->NbrSectorIndicesPerSum[TmpSum];    
	      }
	}
    }
  
  this->NbrNBodySectorSums = this->NbrLzValue;
  this->NbrNBodySectorIndicesPerSum = new int[this->NbrNBodySectorSums];
  for (int i = 0; i < this->NbrNBodySectorSums; ++i)
    this->NbrNBodySectorIndicesPerSum[i] = 0;      
  this->NBodySectorIndicesPerSum = new int* [this->NbrNBodySectorSums];
  switch (this->NBodyValue)
    {
    case 3:
      {
	if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	  {
	    for (int ky1 = 0; ky1 < this->LzMax; ++ky1)
	      for (int ky2 = ky1 + 1; ky2 < this->LzMax; ++ky2) 
		for (int ky3 = ky2 + 1; ky3 <= this->LzMax; ++ky3) 
		  {
		    int TmpSum = (ky1 + ky2 + ky3) % this->NbrLzValue;
		    ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		  }
	    for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	      {
		if (this->NbrNBodySectorIndicesPerSum[i]  > 0)
		  {
		    this->NBodySectorIndicesPerSum[i] = new int[this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i]];      
		    this->NbrNBodySectorIndicesPerSum[i] = 0;
		  }
	      }
	    for (int ky1 = 0; ky1 < this->LzMax; ++ky1)
	      for (int ky2 = ky1 + 1; ky2 < this->LzMax; ++ky2) 
		for (int ky3 = ky2 + 1; ky3 <= this->LzMax; ++ky3) 
		  {
		    int TmpSum = (ky1 + ky2 + ky3) % this->NbrLzValue;
		    this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 3] = ky1;
		    this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = ky2;
		    this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = ky3;
		    ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		  }
	  }
	else
	  {
	    for (int ky1 = 0; ky1 <= this->LzMax; ++ky1)
	      for (int ky2 = ky1; ky2 <= this->LzMax; ++ky2) 
		for (int ky3 = ky2; ky3 <= this->LzMax; ++ky3) 
		  {
		    int TmpSum = (ky1 + ky2 + ky3) % this->NbrLzValue;
		    ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		  }
	    for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	      {
		if (this->NbrNBodySectorIndicesPerSum[i]  > 0)
		  {
		    this->NBodySectorIndicesPerSum[i] = new int[this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i]];      
		    this->NbrNBodySectorIndicesPerSum[i] = 0;
		  }
	      }
	    for (int ky1 = 0; ky1 <= this->LzMax; ++ky1)
	      for (int ky2 = ky1; ky2 <= this->LzMax; ++ky2) 
		for (int ky3 = ky2; ky3 <= this->LzMax; ++ky3) 
		  {
		    int TmpSum = (ky1 + ky2 + ky3) % this->NbrLzValue;
		    this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 3] = ky1;
		    this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = ky2;
		    this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = ky3;
		    ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		  }
	  }
      }
      break;

    case 4:
      {
	if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	  {
	    for (int ky1 = 0; ky1 < this->LzMax; ++ky1)
	      for (int ky2 = ky1 + 1; ky2 < this->LzMax; ++ky2) 
		for (int ky3 = ky2 + 1; ky3 < this->LzMax; ++ky3) 
		  for (int ky4 = ky3 + 1; ky4 <= this->LzMax; ++ky4) 
		    {
		      int TmpSum = (ky1 + ky2 + ky3 + ky4) % this->NbrLzValue;
		      ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		    }
	    for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	      {
		if (this->NbrNBodySectorIndicesPerSum[i]  > 0)
		  {
		    this->NBodySectorIndicesPerSum[i] = new int[this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i]];      
		    this->NbrNBodySectorIndicesPerSum[i] = 0;
		  }
	      }
	    for (int ky1 = 0; ky1 < this->LzMax; ++ky1)
	      for (int ky2 = ky1 + 1; ky2 < this->LzMax; ++ky2) 
		for (int ky3 = ky2 + 1; ky3 < this->LzMax; ++ky3) 
		  for (int ky4 = ky3 + 1; ky4 <= this->LzMax; ++ky4) 
		    {
		      int TmpSum = (ky1 + ky2 + ky3 + ky4) % this->NbrLzValue;
		      this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 4] = ky1;
		      this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 4)] = ky2;
		      this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 4)] = ky3;
		      this->NBodySectorIndicesPerSum[TmpSum][3 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 4)] = ky4;
		      ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		    }
	  }
	else
	  {
	    for (int ky1 = 0; ky1 <= this->LzMax; ++ky1)
	      for (int ky2 = ky1; ky2 <= this->LzMax; ++ky2) 
		for (int ky3 = ky2; ky3 <= this->LzMax; ++ky3) 
		  for (int ky4 = ky3; ky4 <= this->LzMax; ++ky4) 
		    {
		      int TmpSum = (ky1 + ky2 + ky3 + ky4) % this->NbrLzValue;
		      ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		    }
	    for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	      {
		if (this->NbrNBodySectorIndicesPerSum[i]  > 0)
		  {
		    this->NBodySectorIndicesPerSum[i] = new int[this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i]];      
		    this->NbrNBodySectorIndicesPerSum[i] = 0;
		  }
	      }
	    for (int ky1 = 0; ky1 <= this->LzMax; ++ky1)
	      for (int ky2 = ky1; ky2 <= this->LzMax; ++ky2) 
		for (int ky3 = ky2; ky3 <= this->LzMax; ++ky3) 
		  for (int ky4 = ky3; ky4 <= this->LzMax; ++ky4) 
		    {
		      int TmpSum = (ky1 + ky2 + ky3 + ky4) % this->NbrLzValue;
		      this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 4] = ky1;
		      this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 4)] = ky2;
		      this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 4)] = ky3;
		      this->NBodySectorIndicesPerSum[TmpSum][3 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 4)] = ky4;
		      ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		    }
	  }
      }
      break;

    default:
      cout << "warning : " << this->NBodyValue << "-body interaction is not implemented in AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian" << endl;
    }
}
