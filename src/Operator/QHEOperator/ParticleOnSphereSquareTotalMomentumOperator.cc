////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of particle on sphere sqare of total momentum operator         //
//                                                                            //
//                        last modification : 06/03/2005                      //
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
#include "Operator/QHEOperator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


// constructor from default datas
//
// particle = hilbert space associated to the particles
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a fermion

ParticleOnSphereSquareTotalMomentumOperator::ParticleOnSphereSquareTotalMomentumOperator(ParticleOnSphere* particle, int totalLz, int lzMax)
{
  this->Particle = (ParticleOnSphere*) (particle->Clone());
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->Coefficients = RealMatrix(this->LzMax + 1, this->LzMax + 1);
  this->Shift = 0.25 * ((double) (this->TotalLz * this->TotalLz));
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax + 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	{
     	  this->Coefficients(i, j) = TmpCoefficient;
	}
    }
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax - 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	{
     	  this->Coefficients(j, i) *= 0.5 * TmpCoefficient;
	}
    }
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnSphereSquareTotalMomentumOperator::ParticleOnSphereSquareTotalMomentumOperator(const ParticleOnSphereSquareTotalMomentumOperator& oper)
{
  this->Particle = (ParticleOnSphere*) (oper.Particle->Clone());
  this->TotalLz = oper.TotalLz;
  this->LzMax = oper.LzMax;
  this->Coefficients = oper.Coefficients;
  this->Shift = oper.Shift;
}

// destructor
//

ParticleOnSphereSquareTotalMomentumOperator::~ParticleOnSphereSquareTotalMomentumOperator()
{
  delete this->Particle;
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereSquareTotalMomentumOperator::Clone ()
{
  return new ParticleOnSphereSquareTotalMomentumOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereSquareTotalMomentumOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereSquareTotalMomentumOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereSquareTotalMomentumOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereSquareTotalMomentumOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Element = 0.0;
  int Index = 0;
  double Coefficient = 0.0;
  for (int i = 0; i < Dim; ++i)
    {
      for (int j = 0; j < this->LzMax; ++j)
	{
	  for (int  k = 1; k <= this->LzMax; ++k)
	    {
	      Index = this->Particle->AdAdAA(i, k - 1, j + 1, k, j, Coefficient);
	      if (Index != this->Particle->GetHilbertSpaceDimension())
		{
		  Element += V1[Index] * 2.0 * V2[i] * Coefficient * this->Coefficients(j, k);		  
		}
	    }
	}
      Coefficient = this->Coefficients(0, 1) * this->Particle->AdA(i, 0);
      Coefficient += this->Coefficients(this->LzMax - 1, this->LzMax) * this->Particle->AdA(i, this->LzMax);
      for (int k = 1; k < this->LzMax; ++k)
	Coefficient += (this->Coefficients(k, k + 1) + this->Coefficients(k - 1, k)) * this->Particle->AdA(i, k);
      Element += V1[i] * V2[i] * (Coefficient + this->Shift); 
    }
  return Complex(Element);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereSquareTotalMomentumOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}
   
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereSquareTotalMomentumOperator::Multiply(RealVector& vSource, RealVector& vDestination, 
								  int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Index = 0;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      vDestination[i] = vSource[i] * this->Shift;
    }
  for (int i = firstComponent; i < Last; ++i)
    {
      for (int j = 0; j < this->LzMax; ++j)
	{
	  for (int k = 1; k <= this->LzMax; ++k)
	    {
	      Index = this->Particle->AdAdAA(i, k - 1, j + 1, k, j, Coefficient);
	      if (Index != this->Particle->GetHilbertSpaceDimension())
		{
		  vDestination[Index] += 2.0 * vSource[i] * Coefficient * this->Coefficients(j, k);		  
		}
	    }
	}
      Coefficient = this->Coefficients(0, 1) * this->Particle->AdA(i, 0);
      Coefficient += this->Coefficients(this->LzMax - 1, this->LzMax) * this->Particle->AdA(i, this->LzMax);
      for (int k = 1; k < this->LzMax; ++k)
	Coefficient += (this->Coefficients(k, k + 1) + this->Coefficients(k - 1, k)) * this->Particle->AdA(i, k);
      vDestination[i] += vSource[i] * Coefficient;
    }
  return vDestination;
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSphereSquareTotalMomentumOperator::Multiply(ComplexVector& vSource, ComplexVector& vDestination, 
								     int firstComponent, int nbrComponent)
{
  return vDestination;
}


