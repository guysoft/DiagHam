////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of abstract spin chain                       //
//                                                                            //
//                        last modification : 18/04/2001                      //
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


#include "HilbertSpace/AbstractSpinChain.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"

#include <iostream>

using std::cout;
using std::endl;


// virtual destructor
//

AbstractSpinChain::~AbstractSpinChain () 
{
}


// return matrix representation of Sx
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& AbstractSpinChain::Sxi (int i, Matrix& M)
{
  return M;
}

// return matrix representation of i * Sy
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& AbstractSpinChain::Syi (int i, Matrix& M)
{
  return M;
}

// return matrix representation of Sz
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& AbstractSpinChain::Szi (int i, Matrix& M)
{
  return M;
}

// return index of resulting state from application of S+_i S-_j Sz_k operator on a given state
//
// i = position of S+ operator
// j = position of S- operator
// k = position of Sz operator
// state = index of the state to be applied on S+_i S-_j Sz_k operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::SpiSmjSzk (int i, int j, int k, int state, double& coefficient)
{
  cout << "warning : SpiSmjSzk is not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix AbstractSpinChain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  return this->EvaluatePartialDensityMatrix(nbrSites, szSector, 0, groundState, architecture);
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix AbstractSpinChain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  return this->EvaluatePartialDensityMatrix(nbrSites, szSector, 0, groundState, architecture);
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  return this->EvaluatePartialEntanglementMatrix(nbrSites, szSector, 0, groundState, architecture);
}	

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  return this->EvaluatePartialEntanglementMatrix(nbrSites, szSector, 0, groundState, architecture);
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// shift = position of the A part leftmost site within the full system
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix AbstractSpinChain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, int shift, RealVector& groundState, AbstractArchitecture* architecture)
{
  return RealSymmetricMatrix();
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// shift = position of the A part leftmost site within the full system
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix AbstractSpinChain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, int shift, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  return HermitianMatrix();
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// shift = position of the A part leftmost site within the full system
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, int shift, RealVector& groundState, AbstractArchitecture* architecture)
{
  return RealMatrix();
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// shift = position of the A part leftmost site within the full system
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, int shift, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  return ComplexMatrix();
}
	
