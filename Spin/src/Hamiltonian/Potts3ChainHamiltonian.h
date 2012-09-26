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


#ifndef POTTS3CHAINHAMILTONIAN_H
#define POTTS3CHAINHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class Potts3ChainHamiltonian : public AbstractHamiltonian
{

 protected:
  
  AbstractSpinChain* Chain;

  // array containing the coupling constants between two neighboring sites
  double* PotentialTerms;

  // array containing the coupling constants for the on-site flip terms
  double* FlipTerms;

  // number of sites
  int NbrSpin;
  // number of sites minus one
  int ReducedNbrSpin;

  double* SzSzContributions;
 
  // true if the chain is periodic
  bool PeriodicFlag;

 public:

  // constructor from default datas
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spin
  // potentialTerms = array containing the coupling constants between two neighboring sites
  // flipTerms = array containing the coupling constants for the on-site flip terms 
  // periodicFlag = true if the chain is periodic
  Potts3ChainHamiltonian(AbstractSpinChain* chain, int nbrSpin, double* potentialTerms, double* flipTerms, bool periodicFlag);

  // destructor
  //
  ~Potts3ChainHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set chain
  // 
  // chain = pointer on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  Potts3ChainHamiltonian& SetChain(AbstractSpinChain* chain);

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

/*   // evaluate matrix element */
/*   // */
/*   // V1 = vector to left multiply with current matrix */
/*   // V2 = vector to right multiply with current matrix */
/*   // return value = corresponding matrix element */
/*   Complex MatrixElement (RealVector& V1, RealVector& V2); */
  
/*   // evaluate matrix element */
/*   // */
/*   // V1 = vector to left multiply with current matrix */
/*   // V2 = vector to right multiply with current matrix */
/*   // return value = corresponding matrix element */
/*   Complex MatrixElement (ComplexVector& V1, ComplexVector& V2); */

/*   // multiply a vector by the current hamiltonian and store result in another vector */
/*   // low level function (no architecture optimization) */
/*   // */
/*   // vSource = vector to be multiplied */
/*   // vDestination = vector where result has to be stored */
/*   // return value = reference on vectorwhere result has been stored */
/*   RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination); */

/*   // multiply a vector by the current hamiltonian for a given range of idinces  */
/*   // and store result in another vector, low level function (no architecture optimization) */
/*   // */
/*   // vSource = vector to be multiplied */
/*   // vDestination = vector where result has to be stored */
/*   // firstComponent = index of the first component to evaluate */
/*   // nbrComponent = number of components to evaluate */
/*   // return value = reference on vector where result has been stored */
/*   RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination,  */
/* 			       int firstComponent, int nbrComponent); */

/*   // multiply a vector by the current hamiltonian for a given range of indices  */
/*   // and add result to another vector, low level function (no architecture optimization) */
/*   // */
/*   // vSource = vector to be multiplied */
/*   // vDestination = vector at which result has to be added */
/*   // return value = reference on vectorwhere result has been stored */
/*   RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination); */

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
				  int firstComponent, int nbrComponent);

/*   // multiply a vector by the current hamiltonian and store result in another vector */
/*   // low level function (no architecture optimization) */
/*   // */
/*   // vSource = vector to be multiplied */
/*   // vDestination = vector where result has to be stored */
/*   // return value = reference on vectorwhere result has been stored */
/*   ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination); */

/*   // multiply a vector by the current hamiltonian for a given range of indices  */
/*   // and store result in another vector, low level function (no architecture optimization) */
/*   // */
/*   // vSource = vector to be multiplied */
/*   // vDestination = vector where result has to be stored */
/*   // firstComponent = index of the first component to evaluate */
/*   // nbrComponent = number of components to evaluate */
/*   // return value = reference on vector where result has been stored */
/*   ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,  */
/* 				  int firstComponent, int nbrComponent); */

/*   // multiply a vector by the current hamiltonian for a given range of indices  */
/*   // and add result to another vector, low level function (no architecture optimization) */
/*   // */
/*   // vSource = vector to be multiplied */
/*   // vDestination = vector at which result has to be added */
/*   // return value = reference on vectorwhere result has been stored */
/*   ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination); */

/*   // multiply a vector by the current hamiltonian for a given range of indices  */
/*   // and add result to another vector, low level function (no architecture optimization) */
/*   // */
/*   // vSource = vector to be multiplied */
/*   // vDestination = vector at which result has to be added */
/*   // firstComponent = index of the first component to evaluate */
/*   // nbrComponent = number of components to evaluate */
/*   // return value = reference on vector where result has been stored */
/*   ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination,  */
/* 				     int firstComponent, int nbrComponent); */
 
 private:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
