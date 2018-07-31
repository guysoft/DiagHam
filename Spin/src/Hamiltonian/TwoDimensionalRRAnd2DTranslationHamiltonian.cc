////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of two dimension spin model that could host             //
//                a Read-Rezayi Z3 phase with 2d translations                 //
//                                                                            //
//                        last modification : 27/07/2018                      //
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


#include "Hamiltonian/TwoDimensionalRRAnd2DTranslationHamiltonian.h"
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


// default constructor
//

TwoDimensionalRRAnd2DTranslationHamiltonian::TwoDimensionalRRAnd2DTranslationHamiltonian()
{
  this->J1Factor = 0.0;
  this->J2Factor = 0.0;
  this->J3Factor = 0.0;
  this->JcFactor = 0.0;
}

// constructor from default data
//
// chain = pointer to Hilbert space of the associated system
// xMomentum = momentum along the x direction
// nbrSpinX = number of spin along the x direction
// yMomentum = momentum along the y direction
// nbrSpinY = number of spin along the y direction
// j1Factor = amplitude of the Heisenberg coupling between nearest neighbors
// j2Factor = amplitude of the (S_i S_j)^2 nearest neighbor coupling
// j3Factor = amplitude of the (S_i S_j)^3 nearest neighbor coupling
// jcFactor = amplitude of the chiral term

TwoDimensionalRRAnd2DTranslationHamiltonian::TwoDimensionalRRAnd2DTranslationHamiltonian(AbstractSpinChain* chain, int xMomentum, int nbrSpinX, 
											 int yMomentum, int nbrSpinY, double j1Factor, double j2Factor, double j3Factor, double jcFactor)
{
  this->Chain = chain;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->J1Factor = j1Factor;
  this->JFactor = j1Factor;
  this->JzFactor = j1Factor;
  this->J2Factor = j2Factor;
  this->J3Factor = j3Factor;
  this->JcFactor = jcFactor;
  this->HalfJFactor = 0.5 * this->JFactor;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateExponentialFactors();
  this->EvaluateDiagonalMatrixElements();
}

// destructor
//

TwoDimensionalRRAnd2DTranslationHamiltonian::~TwoDimensionalRRAnd2DTranslationHamiltonian() 
{
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& TwoDimensionalRRAnd2DTranslationHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
	      int TmpIndex1 = this->GetLinearizedIndex(j, k);
	      int TmpIndex2 = this->GetSafeLinearizedIndex(j + 1, k);
	      int TmpIndex3 = this->GetSafeLinearizedIndex(j, k + 1);
	      int TmpIndex4 = this->GetSafeLinearizedIndex(j + 1, k + 1);
	      
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex2, i, dim, vDestination, TmpValue);
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex2, TmpIndex1, i, dim, vDestination, TmpValue);
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex1, TmpIndex3, i, dim, vDestination, TmpValue);
	      this->EvaluateOffDiagonalPowerHeisenbergContribution(TmpIndex3, TmpIndex1, i, dim, vDestination, TmpValue);
	      
	      this->EvaluateOffDiagonalChiralContribution(TmpIndex1, TmpIndex2, TmpIndex4, i, dim, vDestination, TmpValue);
	      this->EvaluateOffDiagonalChiralContribution(TmpIndex2, TmpIndex4, TmpIndex3, i, dim, vDestination, TmpValue);
	      this->EvaluateOffDiagonalChiralContribution(TmpIndex4, TmpIndex3, TmpIndex1, i, dim, vDestination, TmpValue);
	      this->EvaluateOffDiagonalChiralContribution(TmpIndex3, TmpIndex1, TmpIndex2, i, dim, vDestination, TmpValue);
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

ComplexVector* TwoDimensionalRRAnd2DTranslationHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
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

void TwoDimensionalRRAnd2DTranslationHamiltonian::EvaluateDiagonalMatrixElements()
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
	      double Tmp2;
	      Tmp2 = this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					 this->GetSafeLinearizedIndex(j, k + 1), i);
	      Tmp += Tmp2 * (Tmp2 * (Tmp2 * this->J3Factor + this->J2Factor) + this->J1Factor);
	      Tmp2 = this->Chain->SziSzj(this->GetLinearizedIndex(j, k), 
					 this->GetSafeLinearizedIndex(j + 1, k), i);
	      Tmp += Tmp2 * (Tmp2 * (Tmp2 * this->J3Factor + this->J2Factor) + this->J1Factor);
	    }
	}
      this->SzSzContributions[i] += Tmp;
   }
}

