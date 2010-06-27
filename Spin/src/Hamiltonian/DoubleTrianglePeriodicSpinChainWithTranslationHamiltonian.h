////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of double triangle spin chain  hamiltonian           //
//                with periodic boundary conditions and translations          //
//                                                                            //
//                        last modification : 25/06/2010                      //
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


#ifndef DOUBLETRIANGLEPERIODICSPINCHAINWITHTRANSLATIONHAMILTONIAN_H
#define DOUBLETRIANGLEPERIODICSPINCHAINWITHTRANSLATIONHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian : public AbstractHamiltonian
{

 protected:
  
  // pointer to the Hilbert space
  AbstractSpinChain* Chain;

  // nearest neighbour coupling constant and half its value
  double J1;
  double HalfJ1;
  // second nearest neighbour coupling constant and half its value
  double J2;
  double HalfJ2;

  // number of spins
  int NbrSpin;
  // number of triangles
  int NbrTriangles;

  // temporary stoage for the Hamiltonian diagonal elements
  double* SzSzContributions;

  //array containing all the cosinus that are needed when computing matrix elements
  double* CosinusTable;
  //array containing all the sinus that are needed when computing matrix elements
  double* SinusTable;

 public:

  // constructor from default datas
  //
  // chain = reference on Hilbert space of the associated system
  // nbrSpin = number of spins
  // j1 = nearest neighbour coupling constant
  // j2 = second nearest neighbour coupling constant
  DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian(AbstractSpinChain* chain, int nbrSpin, double j1, double j2);

  // destructor
  //
  ~DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian();

  // set chain
  // 
  // chain = pointer on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian& SetChain(AbstractSpinChain* chain);

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

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
  ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors,
					     int firstComponent, int nbrComponent);


 private:
 
  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

  // evaluate all cosinus/sinus that are needed when computing matrix elements
  //
  void EvaluateCosinusTable();

};

#endif
