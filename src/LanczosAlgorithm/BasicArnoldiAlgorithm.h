////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of basic  Arnoldi algorithm                     //
//                         for non symmetric matrices                         //
//                                                                            //
//                        last modification : 17/11/2012                      //
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


#ifndef BASICARNOLDIALGORITHM_H
#define BASICARNOLDIALGORITHM_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/ComplexVector.h"
#include "GeneralTools/GarbageFlag.h"


class BasicArnoldiAlgorithm : public AbstractLanczosAlgorithm
{

 protected:

  // array where the vectors of the Krylov subspace
  ComplexVector* ArnoldiVectors;

  // maximum  number of iterations
  int MaximumNumberIteration;
  // current iteration index
  int Index;

  // garbage flag to avoid duplicating memory area
  GarbageFlag Flag;

  // Hessenberg matrix of the Hamiltonian in the Krylov subspace
  ComplexMatrix ReducedMatrix;
  // temporary matrix used to duplicated ReducedMatrix before diagonalize it
  ComplexMatrix TemporaryReducedMatrix;

  ComplexDiagonalMatrix ComplexDiagonalizedMatrix;

  // number of wanted eigenvalues
  int NbrEigenvalue;
  // value of the last wanted eigenvalue at previous Arnoldi iteration
  double PreviousLastWantedEigenvalue;
  // value of the wanted eigenvalue at previous Arnoldi iteration
  double* PreviousWantedEigenvalues;
  // flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  bool StrongConvergenceFlag;

  // array used to store temporary scalar products
  Complex* TemporaryCoefficients;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues
  // maxIter = an approximation of maximal number of iteration
  // strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  BasicArnoldiAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int maxIter = 100,
			bool strongConvergence = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  BasicArnoldiAlgorithm(const BasicArnoldiAlgorithm& algorithm);

  // destructor
  //
  ~BasicArnoldiAlgorithm();

  // initialize Arnoldi algorithm with a random vector
  //
  void InitializeLanczosAlgorithm();
  
  // initialize Arnoldi algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  void InitializeLanczosAlgorithm(const Vector& vector);

  // get last produced vector
  //
  // return value = reference on last produced vector
  Vector& GetGroundState();

  // get the n first eigenstates
  //
  // nbrEigenstates = number of needed eigenstates
  // return value = array containing the eigenstates
  Vector* GetEigenstates(int nbrEigenstates);

  // run current Arnoldi algorithm (continue from previous results if Arnoldi algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  void RunLanczosAlgorithm (int nbrIter);
  
  // test if convergence has been reached
  //
  // return value = true if convergence has been reached
  bool TestConvergence ();

 protected:

  // diagonalize tridiagonalized matrix and find ground state energy
  //
  virtual void Diagonalize();

};

#endif
