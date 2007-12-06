////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar M�ller                  //
//                                                                            //
//                                                                            //
//       class implementing generalized Halperin states on the sphere //
//                                                                            //
//                        last modification : 06/12/2007                      //
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


#ifndef EXTENDEDHALPERINWAVEFUNCTION_H
#define EXTENDEDHALPERINWAVEFUNCTION_H


#ifdef HAVE_LAPACK
// you may opt out of using the LAPACK determinant routines by quoting the definition below
#define USE_LAPACK_CFCB
#endif


// use the following flag to move jastrow factors inside the Cauchy determinant if possible
#define JASTROW_INSIDE

#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereOrbitals.h"

class ExtendedHalperinWavefunction: public Abstract1DComplexFunction
{

 protected:

  int NbrParticles;
  int NbrParticlesPerLayer;
  int K_outside;
  int K_inside;
  int M_outside;
  int M_inside;
  int Q;
  int R;
  
  // flag to indicate whether Jastrow factors should be moved inside Cauchy determinant.
  bool JastrowInside;

  // single-particle Jastrow factors
  Complex *J11;
  Complex *J12;
  Complex *J21;
  Complex *J22;


  bool HaveDeterminant;

  // factor used to multiply each element of the Cauchy matrix
  double DeterminantNorm;
  // factor used to multiply each element of the Cauchy determinant
  double PermanentNorm;
  // factor used to multiply each term within Jastrow-factors
  double JastrowNorm;

  // temporary arrays used to store (u_i v_j - u_j v_i) factors
  Complex** JastrowFactorElements;
  // temporary array used to store u spinor coordinates
  Complex* SpinorUCoordinates;
  // temporary array used to store v spinor coordinates
  Complex* SpinorVCoordinates;

  // matrix for calculation of cauchy determinant
#ifdef USE_LAPACK_CFCB
  ComplexLapackDeterminant *Matrix;
#else
  ComplexMatrix *Matrix;
#endif  

  // matrix for calculation of cauchy permanent
  ComplexMatrix *Matrix2;
  
 public:

  // default constructor
  // 
  ExtendedHalperinWavefunction();

  // constructor
  //
  // nbrParticles = number of particles per layer
  // k = Jastrow factors between same species
  // m = Jastrow factors between different species
  // q = power inside Cauchy-determinant
  // r = power inside Cauchy-permanent
  // moveJastrowInside = move jastrow-factors inside cauchy determinant
  ExtendedHalperinWavefunction(int nbrParticles, int k, int m, int q, int r=0, bool moveJastrowInside = false );

  // copy constructor
  //
  // function = reference on the wave function to copy
  ExtendedHalperinWavefunction(const ExtendedHalperinWavefunction& function);

  // destructor
  //
  ~ExtendedHalperinWavefunction();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  virtual Complex operator ()(RealVector& x);

  // normalize the wave-function to one for the given particle positions
  // x = point where the function has to be evaluated
  void AdaptNorm(RealVector& x);

  // normalize the wave-function over an average number of MC positions
  void AdaptAverageMCNorm(int thermalize=100, int average=250);

 private:

  // do all precalculation operations required for a new set of positions

  void EvaluateTables(RealVector& x);

};


#endif  //  EXTENDEDHALPERINWAVEFUNCTION_H
