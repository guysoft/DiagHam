////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                 class of quatum Hall hamiltonian associated                //
//   to particles with contact interactions on a lattice in magnetic field    //
//                                                                            //
//                      last modification : 13/02/2008                        //
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


#ifndef GROSSPITAEVSKIIONLATTICESTATE_H
#define GROSSPITAEVSKIIONLATTICESTATE_H


#include "config.h"
#include "Tools/FQHESpectrum/LatticePhases.h"
#include "MathTools/Complex.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include <iostream>


using std::ostream;

class AbstractRandomNumberGenerator;

class GrossPitaevskiiOnLatticeState
{
 protected:
  // number of sites
  int NbrSites;

  // number of hopping terms
  int NbrHoppingTerms;

  // values of hopping terms
  Complex *HoppingTerms;
  int* KineticQi;
  int* KineticQf;

  // number of 2-particle interactions
  int NbrInteractionFactors;

  // values of 2-particle interactions
  Complex *InteractionFactors;
  int* Q1Value;
  int* Q2Value;
  int* Q3Value;
  int* Q4Value;

  // filename for external one-particle terms
  char *OneParticleTerms;
  
  // filename for external two-particle interactions
  char *TwoParticleTerms;

  // lattice structure (if known)
  LatticePhases *LatticeGeometry;

  // variational parameters
  RealVector VariationalParameters;

  // internal coefficients
  ComplexVector VariationalCoefficients;

  // number of evaluations in optimization
  int NbrEvaluations;

  // random number generator
  AbstractRandomNumberGenerator *RandomNumbers;

  // value of the chemical potential
  double ChemicalPotential;

  
 public:

  // constructor for contact interactions on a square lattice
  //
  // nbrStates = number of quantum states
  // oneParticleTerms = file describing single particle terms
  // twoParticleTerms = file describing two-particle terms
  GrossPitaevskiiOnLatticeState(int nbrStates, const char* oneParticleTerms, const char* twoParticleTerms, LatticePhases *latticeGeometry = NULL, RealVector *variationalParameters=NULL);

  // destructor
  //
  ~GrossPitaevskiiOnLatticeState();

  // get the parameters of the Many-Body state that was last calculated
  // return = state
  RealVector & GetVariationalParameters() {return this->VariationalParameters;}

  // set trial parameters
  void SetVariationalParameters(RealVector &variationalParameters);

  // set chemical potential
  void SetChemicalPotential(double chemicalPotential) {this->ChemicalPotential=chemicalPotential;}

  // set parameters to a random initial distribution (random phase)
  // amplitude = amplitude determining the density
  void SetToRandomPhase(double amplitude=1.0);

  // get expectation value of the energy
  double GetEnergy();

  // get the total number of particles corresponding to the last configuration
  double GetNbrParticles();

  // optimize wavefunction starting from present settings of VariationalParameters
  // tolerance = final tolerance on the variational parameters
  // maxIter = maximal number of function evaluations
  //
  double Optimize(double tolerance, int maxIter);

 private:

  // target function for optimizer routine:
  double EvaluateEnergy(int nbrParameters, double *x);
  
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();
  
};

#endif // GROSSPITAEVSKIIONLATTICESTATE_H
