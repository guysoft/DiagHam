////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of two dimension Heisenberg model                  //
//                             and 2d translations                            //
//                                                                            //
//                        last modification : 13/05/2018                      //
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


#include "Hamiltonian/TwoDimensionalHeisenbergAnd2DTranslationHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/RealVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// xMomentum = momentum along the x direction
// nbrSpinX = number of spin along the x direction
// yMomentum = momentum along the y direction
// nbrSpinY = number of spin along the y direction
// jFactor = Heisenberg XX coupling constant between nearest neighbors
// jzFactor = Heisenberg Z coupling constant between nearest neighbors

TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::TwoDimensionalHeisenbergAnd2DTranslationHamiltonian(AbstractSpinChain* chain, int xMomentum, int nbrSpinX, 
											   int yMomentum, int nbrSpinY, double jFactor, double jzFactor)
{
  this->Chain = chain;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->JFactor = jFactor;
  this->JzFactor = jzFactor;
  this->HalfJFactor = 0.5 * this->JFactor;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateExponentialFactors();
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::~TwoDimensionalHeisenbergAnd2DTranslationHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChain*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::ShiftHamiltonian (double shift)
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

ComplexVector& TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
											int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];
      Complex& TmpValue = vSource[i];
       for (int j = 0; j < this->NbrSpinX; j++)
 	{
 	  for (int k = 0; k < this->NbrSpinY; k++)
 	    {
	      int Index0 = this->GetLinearizedIndex(j, k);
	      int IndexX0 = this->GetSafeLinearizedIndex(j + 1, k);
	      int IndexY0 = this->GetSafeLinearizedIndex(j, k + 1);
 	      pos = this->Chain->SmiSpj(Index0, IndexX0, i, coef, NbrTranslationsX, NbrTranslationsY);
 	      if (pos != dim)
 		{
 		  vDestination[pos] += (coef * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
 		}
 	      pos = this->Chain->SmiSpj(IndexX0, Index0, i, coef, NbrTranslationsX, NbrTranslationsY);
 	      if (pos != dim)
 		{
 		  vDestination[pos] += (coef * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
 		}
 	      pos = this->Chain->SmiSpj(Index0, IndexY0, i, coef, NbrTranslationsX, NbrTranslationsY);
 	      if (pos != dim)
 		{
 		  vDestination[pos] += (coef * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
 		}
 	      pos = this->Chain->SmiSpj(IndexY0, Index0, i, coef, NbrTranslationsX, NbrTranslationsY);
 	      if (pos != dim)
 		{
 		  vDestination[pos] += (coef * this->HalfJFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValue;
 		}
	    }
	}
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
												int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  Complex* TmpValues = new Complex[nbrVectors];
  int NbrTranslationsX;
  int NbrTranslationsY;
  for (int k = 0; k < nbrVectors; ++k)
    {
      ComplexVector& TmpSource = vSources[k];
      ComplexVector& TmpDestination = vDestinations[k];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  TmpDestination[i] += this->SzSzContributions[i] * TmpSource[i];
	}
    }
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int k = 0; k < nbrVectors; ++k)
	{
	  TmpValues[k] = vSources[k][i];
	}
//       for (int j = 0; j < this->NbrSpinX; j++)
// 	{
// 	  for (int k = 0; k < this->NbrSpinY; k++)
// 	    {
// 	      pos = this->Chain->Spi(this->GetLinearizedIndex(j, k), i, coef, NbrTranslationsX, NbrTranslationsY);
// 	      if (pos != dim)
// 		{
// 		  for (int l = 0; l < nbrVectors; ++l)
// 		    vDestinations[l][pos] += (coef * this->HxFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValues[l];
// 		}
// 	      pos = this->Chain->Smi(this->GetLinearizedIndex(j, k), i, coef, NbrTranslationsX, NbrTranslationsY);
// 	      if (pos != dim)
// 		{
// 		  for (int l = 0; l < nbrVectors; ++l)
// 		    vDestinations[l][pos] += (coef * this->HxFactor) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * TmpValues[l];
// 		}
// 	    }
// 	}
    }
  delete[] TmpValues;
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  // SzSz part
  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      this->SzSzContributions[i] = 0.0;
      double Tmp = 0.0;
      for (int j = 0; j < this->NbrSpinX; j++)
	{
	  for (int k = 0; k < this->NbrSpinY; k++)
	    {
	      Tmp += this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					 this->GetSafeLinearizedIndex(j + 1, k), i);
	      Tmp += this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					 this->GetSafeLinearizedIndex(j, k + 1), i);
	    }
	}
      Tmp *= this->JzFactor; 
      this->SzSzContributions[i] += Tmp;
   }
}

// evaluate all exponential factors
//   

void TwoDimensionalHeisenbergAnd2DTranslationHamiltonian::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex*[this->NbrSpinX];
  for (int i = 0; i < this->NbrSpinX; ++i)
    { 
      this->ExponentialFactors[i] = new Complex[this->NbrSpinY];
      for (int j = 0; j < this->NbrSpinY; ++j)
	{ 
	  this->ExponentialFactors[i][j] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->NbrSpinX))
							       + (this->YMomentum * ((double) j) / ((double) this->NbrSpinY))));
	}
    }
}
