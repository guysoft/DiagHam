////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract fractional quantum Hall hamiltonian          //
//              associated to particles with SU(2) spin on a sphere           //
//                                                                            //
//                        last modification : 07/06/2007                      //
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


#ifndef ABSTRACTQHEONSPHEREWITHSPINHAMILTONIAN_H
#define ABSTRACTQHEONSPHEREWITHSPINHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"

#include <iostream>


using std::ostream;


class AbstractQHEOnSphereWithSpinHamiltonian : public AbstractQHEOnSphereHamiltonian
{

 protected:
  
  // number of index sum for the creation (or annhilation) operators for the intra spin sector
  int NbrIntraSectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the intra spin sector
  int* NbrIntraSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators for the intra spin sector
  int** IntraSectorIndicesPerSum;

  // number of index sum for the creation (or annhilation) operators for the inter spin sector
  int NbrInterSectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the inter spin sector
  int* NbrInterSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators for the inter spin sector
  int** InterSectorIndicesPerSum;

  // alternative method used to store indices attached to each interaction factor
  //  arrays of m1 and m2 intra-sector indices (i.e. indices for the a_m1_s a_m2_s factors)
  int* M1IntraValue;
  int* M2IntraValue;  
  //  number of element in each m1 or m2 intra-sector index array
  int NbrM12IntraIndices;
  //  the m3 intra-sector index array (i.e. index for the ad_m3_s a_(m1+m3-m3)_s factors) for each set of (m1,m2)
  int** M3IntraValues;
  // number of possible m3 intra-sector index for each set of (m1,m2)
  int* NbrM3IntraValues;

  // alternative method used to store indices attached to each interaction factor
  //  arrays of m1 and m2 inter-sector indices (i.e. indices for the a_m1_s a_m2_-s factors)
  int* M1InterValue;
  int* M2InterValue;  
  //  number of element in each m1 or m2 inter-sector index array
  int NbrM12InterIndices;
  //  the m3 inter-sector index array (i.e. index for the ad_m3_s a_(m1+m3-m3)_-s factors) for each set of (m1,m2)
  int** M3InterValues;
  // number of possible m3 inter-sector index for each set of (m1,m2)
  int* NbrM3InterValues;

  // arrays containing all interaction factors, the first index correspond correspond to index sum for the creation (or annhilation) operators
  // the second index is a linearized index (m1,m2) + (n1,n2) * (nbr element in current index sum) (m for creation operators, n for annhilation operators)
  // array containing all interaction factors for spin up and spin up
  double** InteractionFactorsupup;
  // array containing all interaction factors for spin down and spin down
  double** InteractionFactorsdowndown;
  // array containing all interaction factors for spin up and spin down
  double** InteractionFactorsupdown;

  // arrays containing all interaction factors needed for the alternative method
  // array containing all interaction factors for spin up and spin up
  double* M12InteractionFactorsupup;
  // array containing all interaction factors for spin down and spin down
  double* M12InteractionFactorsdowndown;
  // array containing all interaction factors for spin up and spin down
  double* M12InteractionFactorsupdown;

  // array that contains all one-body interaction factors for particles with spin up
  double* OneBodyInteractionFactorsupup;
  // array that contains all one-body interaction factors for particles with spin down
  double* OneBodyInteractionFactorsdowndown;
  

 public:

  // destructor
  //
  virtual ~AbstractQHEOnSphereWithSpinHamiltonian() = 0;

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					  int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual  RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						   int firstComponent, int nbrComponent);

 protected:
 
  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  RealVector* LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							     int firstComponent, int nbrComponent);

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication();

  // enable fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm
  //
  // fileName = prefix of the name of the file where temporary matrix elements will be stored
  virtual void EnableFastMultiplicationWithDiskStorage(char* fileName);

};

#endif
