////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of Potts 3 chain hamiltonian                     //
//                                                                            //
//                        last modification : 04/06/2012                      //
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


#include "Hamiltonian/Potts3ChainHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// potentialTerms = array containing the coupling constants between two neighboring sites
// flipTerms = array containing the coupling constants for the on-site flip terms 
// periodicFlag = true if the chain is periodic

Potts3ChainHamiltonian::Potts3ChainHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* potentialTerms, double* flipTerms, bool periodicFlag)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->ReducedNbrSpin = this->NbrSpin -1;
  this->PeriodicFlag = periodicFlag;
  if (this->PeriodicFlag == false)
    this->PotentialTerms = new double [this->ReducedNbrSpin];
  else
    this->PotentialTerms = new double [this->NbrSpin];
  this->FlipTerms = new double [this->NbrSpin];
  for (int i = 0; i < this->ReducedNbrSpin; i++)
    {
      this->PotentialTerms[i] = potentialTerms[i];
    }
  if (this->PeriodicFlag == true)
    this->PotentialTerms[this->ReducedNbrSpin] = potentialTerms[this->ReducedNbrSpin];
   for (int i = 0; i < this->NbrSpin; i++)
    {
      this->FlipTerms[i] = flipTerms[i];
    }
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

Potts3ChainHamiltonian::~Potts3ChainHamiltonian() 
{
  delete[] this->PotentialTerms;
  delete[] this->FlipTerms;
  delete[] this->SzSzContributions;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void Potts3ChainHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* Potts3ChainHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

Potts3ChainHamiltonian& Potts3ChainHamiltonian::SetChain(AbstractSpinChain* chain)
{  
  delete[] this->SzSzContributions;
  this->Chain = chain;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  return *this;  
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int Potts3ChainHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void Potts3ChainHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& Potts3ChainHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
							int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  int pos;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      double TmpValue = vSource[i];
      vDestination[i] += this->SzSzContributions[i] * TmpValue;
      for (int j = 0; j < this->ReducedNbrSpin; ++j)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->PotentialTerms[j] * TmpValue;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += this->PotentialTerms[j] * TmpValue;
	    }
	}
    }
  if (this->PeriodicFlag == true)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  double TmpValue = this->PotentialTerms[this->ReducedNbrSpin] * vSource[i];
	  pos = this->Chain->SmiSpj(0, this->ReducedNbrSpin, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += TmpValue;
	    }
	  pos = this->Chain->SmiSpj(this->ReducedNbrSpin, 0, i, coef);
	  if (pos != dim)
	    {
	      vDestination[pos] += TmpValue;
	    }
	}
    }
  return vDestination;
}

// evaluate diagonal matrix elements
// 

void Potts3ChainHamiltonian::EvaluateDiagonalMatrixElements()
{
  int Dimension = this->Chain->GetHilbertSpaceDimension();
  double Coefficient;
  for (int i = 0; i < Dimension; i++)
    {
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < this->NbrSpin; j++)
	{
	  this->Chain->Szi(j, i, Coefficient);
	  if (Coefficient == 0.0)
	    this->SzSzContributions[i] += 2.0 * this->FlipTerms[j];
	  else
	    this->SzSzContributions[i] -= this->FlipTerms[j];
	}
    }
}

