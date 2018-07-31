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


#ifndef TWODIMENSIONALRRAND2DTRANSLATIONISINGHAMILTONIAN_H
#define TWODIMENSIONALRRAND2DTRANSLATIONISINGHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/TwoDimensionalHeisenbergAnd2DTranslationHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class TwoDimensionalRRAnd2DTranslationHamiltonian : public TwoDimensionalHeisenbergAnd2DTranslationHamiltonian
{

 protected:
  
  // amplitude of the Heisenberg coupling between nearest neighbors
  double J1Factor;
  // amplitude of the (S_i S_j)^2 nearest neighbor coupling
  double J2Factor;
  // amplitude of the (S_i S_j)^3 nearest neighbor coupling
  double J3Factor;
  // amplitude of the chiral term
  double JcFactor;

 public:

  // default constructor
  //
  TwoDimensionalRRAnd2DTranslationHamiltonian();

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
  TwoDimensionalRRAnd2DTranslationHamiltonian(AbstractSpinChain* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, 
					      double j1Factor, double j2Factor, double j3Factor, double jcFactor);

  // destructor
  //
  ~TwoDimensionalRRAnd2DTranslationHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						     int firstComponent, int nbrComponent);

 protected:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

  // evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
  //
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // index = index of the many-body state to act on
  // dimension = total Hilbert space dimension
  // vDestination = vector at which result has to be added
  // coefficient = global multiplicative coefficient
  virtual void EvaluateOffDiagonalPowerHeisenbergContribution(int i, int j, int index, int dimension, ComplexVector& vDestination, Complex& coefficient);

  // evaluate the off-diagonal chiral contribution for a single term ( S_i (S_j ^ S_k) )
  //
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // k = linearized position of the second spin
  // index = index of the many-body state to act on
  // dimension = total Hilbert space dimension
  // vDestination = vector at which result has to be added
  // coefficient = global multiplicative coefficient
  virtual void EvaluateOffDiagonalChiralContribution(int i, int j, int k, int index, int dimension, ComplexVector& vDestination, Complex& coefficient);

};

// evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
//
// i = linearized position of the first spin
// j = linearized position of the second spin
// index = index of the many-body state to act on
// dimension = total Hilbert space dimension
// vDestination = vector at which result has to be added
// coefficient = global multiplicative coefficient

inline void TwoDimensionalRRAnd2DTranslationHamiltonian::EvaluateOffDiagonalPowerHeisenbergContribution(int i, int j, int index, int dimension, ComplexVector& vDestination, Complex& coefficient)
{
  double TmpCoefficient;
  double TmpCoefficient2;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int pos = this->Chain->SmiSpj(i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J1Factor * TmpValue2;
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      TmpValue2 *= TmpValue3;
      vDestination[pos] += (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
    }

  pos = this->Chain->SmiSpjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      vDestination[pos] +=  this->J3Factor * TmpValue2 * TmpValue3;
    }
  pos = this->Chain->SmiSpjSmkSpl(j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      vDestination[pos] +=  this->J3Factor * TmpValue2 * TmpValue3;
    }
  pos = this->Chain->SziSzjSmkSpl(i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      vDestination[pos] +=  this->J3Factor * TmpValue2 * TmpValue3;
    }

  pos = this->Chain->SmiSpjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 =  (0.125 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J3Factor * TmpValue2;	      
    }
  pos = this->Chain->SmiSpjSmkSplSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 =  (0.125 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J3Factor * TmpValue2;	      
    }
  pos = this->Chain->SmiSpjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 =  (0.125 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J3Factor * TmpValue2;	      
    }
  pos = this->Chain->SmiSpjSmkSplSmmSpn(j, i, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 =  (0.125 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] += this->J3Factor * TmpValue2;	      
    }

  pos = this->Chain->SziSzjSmkSplSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] +=  this->J3Factor * TmpValue2;	      
    }
  pos = this->Chain->SziSzjSmkSplSmmSpn(i, j, j, i, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] +=  this->J3Factor * TmpValue2;	      
    }

  pos = this->Chain->SmiSpjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] +=  this->J3Factor * TmpValue2;	      
    }
  pos = this->Chain->SmiSpjSzkSzlSmmSpn(j, i, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] +=  this->J3Factor * TmpValue2;	      
    }

  pos = this->Chain->SziSzjSzkSzlSmmSpn(i, j, i, j, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.5 * TmpCoefficient) * this->ExponentialFactors[NbrTranslationsX][NbrTranslationsY] * coefficient;
      vDestination[pos] +=  this->J3Factor * TmpValue2;	      
    }
}

// evaluate the off-diagonal chiral contribution for a single term ( S_i (S_j ^ S_k) )
//
// i = linearized position of the first spin
// j = linearized position of the second spin
// k = linearized position of the second spin
// index = index of the many-body state to act on
// dimension = total Hilbert space dimension
// vDestination = vector at which result has to be added
// coefficient = global multiplicative coefficient

inline void TwoDimensionalRRAnd2DTranslationHamiltonian::EvaluateOffDiagonalChiralContribution(int i, int j, int k, int index, int dimension, ComplexVector& vDestination, Complex& coefficient)
{
  double TmpCoefficient;
  int NbrTranslationsX;
  int NbrTranslationsY;
  int pos = this->Chain->SpiSmjSzk(i, j, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.5;
      TmpCoefficient *= this->JcFactor;
      Complex& Tmp = vDestination[pos];
      vDestination[pos].Re -= TmpCoefficient * coefficient.Im;
      vDestination[pos].Im += TmpCoefficient * coefficient.Re;
    }
  pos = this->Chain->SpiSmjSzk(i, k, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.5;
      TmpCoefficient *= this->JcFactor;
      Complex& Tmp = vDestination[pos];
      vDestination[pos].Re += TmpCoefficient * coefficient.Im;
      vDestination[pos].Im -= TmpCoefficient * coefficient.Re;
    }
  pos = this->Chain->SpiSmjSzk(j, k, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.5;
      TmpCoefficient *= this->JcFactor;
      Complex& Tmp = vDestination[pos];
      vDestination[pos].Re -= TmpCoefficient * coefficient.Im;
      vDestination[pos].Im += TmpCoefficient * coefficient.Re;
    }
  pos = this->Chain->SpiSmjSzk(k, j, i, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.5;
      TmpCoefficient *= this->JcFactor;
      Complex& Tmp = vDestination[pos];
      vDestination[pos].Re += TmpCoefficient * coefficient.Im;
      vDestination[pos].Im -= TmpCoefficient * coefficient.Re;
    }
  pos = this->Chain->SpiSmjSzk(k, i, j, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.5;
      TmpCoefficient *= this->JcFactor;
      Complex& Tmp = vDestination[pos];
      vDestination[pos].Re -= TmpCoefficient * coefficient.Im;
      vDestination[pos].Im += TmpCoefficient * coefficient.Re;
    }
  pos = this->Chain->SpiSmjSzk(j, i, k, index, TmpCoefficient, NbrTranslationsX, NbrTranslationsY);
  if (pos != dimension)
    {
      TmpCoefficient *= 0.5;
      TmpCoefficient *= this->JcFactor;
      Complex& Tmp = vDestination[pos];
      vDestination[pos].Re += TmpCoefficient * coefficient.Im;
      vDestination[pos].Im -= TmpCoefficient * coefficient.Re;
    }
}

#endif
