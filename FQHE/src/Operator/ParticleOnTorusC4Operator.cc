////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of particle on torus C4 operator                  //
//                                                                            //
//                        last modification : 11/02/2014                      //
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
#include "Operator/ParticleOnTorusC4Operator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
using std::cout;
using std::endl;


// constructor from default data
//
// particle = hilbert space associated to the particles
// piPiSectorFlag = indicates if we are considering the (pi,pi) momentum sector  
// c2Only = consider the C2 operator instead of the C4 operator

ParticleOnTorusC4Operator::ParticleOnTorusC4Operator(ParticleOnTorusWithMagneticTranslations* particle, bool piPiSectorFlag, bool c2Only)
{
  this->Particle= (ParticleOnTorusWithMagneticTranslations*) (particle->Clone());
  this->C2OnlyFlag = c2Only;
  this->PiPiSectorFlag = piPiSectorFlag;
}

// copy constructor
//
// oper = operator to copy
  
ParticleOnTorusC4Operator::ParticleOnTorusC4Operator(ParticleOnTorusC4Operator& oper)
{
  this->Particle = (ParticleOnTorusWithMagneticTranslations*) (oper.Particle->Clone());
  this->C2OnlyFlag = oper.C2OnlyFlag;
  this->PiPiSectorFlag = oper.PiPiSectorFlag;
}

// destructor
//

ParticleOnTorusC4Operator::~ParticleOnTorusC4Operator()
{
  delete this->Particle;
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnTorusC4Operator::Clone ()
{
  return new ParticleOnTorusC4Operator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusC4Operator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnTorusC4Operator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnTorusC4Operator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnTorusC4Operator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  Complex Element = 0.0;
  int Last = firstComponent + nbrComponent;
  for (int i = firstComponent; i < Last; ++i)
    {
//      Element += Conj(V1[this->Particle->ApplyXMagneticTranslation(i)]) * V2[i];
    }
  return Element;
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnTorusC4Operator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							      int firstComponent, int nbrComponent)
{
  int NbrTranslations;
  int TmpIndex;
  if (this->PiPiSectorFlag == false)
    {
      if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
	{
	  int Last = firstComponent + nbrComponent;;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      
	      vDestination[this->Particle->GetC2SymmetricState(i, NbrTranslations)] += vSource[i];
	    }
	}
      else
	{
	  long Last = ((long) firstComponent) + ((long) nbrComponent);
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      vDestination[this->Particle->GetC2SymmetricState(i, NbrTranslations)] += vSource[i];
	    }
	}
    }
  else
    {
      if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
	{
	  int Last = firstComponent + nbrComponent;;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      int Index = this->Particle->GetC2SymmetricState(i, NbrTranslations);
	      if (Index < this->Particle->GetHilbertSpaceDimension())
		{
		  if ((NbrTranslations & 1) == 0)
		    vDestination[Index] += vSource[i];
		  else
		    vDestination[Index] -= vSource[i];
		}
	    }
	}
      else
	{
	  long Last = ((long) firstComponent) + ((long) nbrComponent);
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      long Index = this->Particle->GetC2SymmetricState(i, NbrTranslations);
	      if (Index < this->Particle->GetLargeHilbertSpaceDimension())
		{
		  if ((NbrTranslations & 1) == 0)
		    vDestination[Index] += vSource[i];
		  else
		    vDestination[Index] -= vSource[i];
		}
	    }
	}
    }
  return vDestination;
}
  


