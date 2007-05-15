////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//         Copyright (C) 2001-2002 Nicolas Regnault and Gunnar Moeller        //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 16/09/2004                      //
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


#ifndef JAINCFFILLEDLEVELONSPHEREWAVEFUNCTION_H
#define JAINCFFILLEDLEVELONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Matrix/ComplexMatrix.h"
#include "MathTools/Complex.h"

#include "Matrix/Double3DArray.h"
#include "Matrix/Complex3DArray.h"
#include "Matrix/Complex4DArray.h"

class DerivativeProductFactor;
class DerivativeProduct;
class SumDerivativeProduct;

class JainCFOnSphereOrbitals
{

 private:
  
  friend class DerivativeProductFactor;

 protected:

  // number of particles
  int NbrParticles;

  // number of orbitals
  int NbrOrbitals;

  // number of Landau levels filled with composite fermions
  int NbrLandauLevels;

  // twice the value of the momentum in the lowest pseudo-Landau level
  int TwiceS;
  // falg to indicate if flux  is reversed
  bool ReverseFluxFlag;

  // power to which the Jastrow factor has to be raised
  int ActualJastrowPower;

  // for fermionic case, half the power to which the Jastrow factor has to be raised, for bosonic case, half the power + 1 to which the Jastrow factor has to be raised
  int JastrowPower;

  // number of spinor powers to be calculated (differs for positive and negative flux)
  int MaxSpinorPower;

  // number of derivatives to be calculated (differs for positive and negative flux)
  int MaxDerivativeNum;

  // number of derivatives to be calculated for the LLL
  int LLLDerivativeNum;

  // maximal power to which the derivative has to be raised
  int **MaxDerivativePower;

  // array with factorial and sign prefactors to the DerivativeFactors
  double* FactorialSignFactors;

  // array containing prefactors of each projected monopole harmonic
  // indices: [LL-index n][Momentum m (0=-S-n)]
  double** NormalizationPrefactors;

  // array containing constant factors that appears in the sum of projected monopole harmonic (except LLL)
  // indices: [LL-index n][s from Jains notation(summation index)][momentum m (0=-S-n)]
  // factor nonzero for n-s<= m <= MaxMomentum - s
  Double3DArray SumPrefactors;

  // garbage flag to avoid duplication of precalculation array
  GarbageFlag Flag;

  // temporary array used to store (u_i v_j - u_j v_i)^-1 factors
  Complex** JastrowFactorElements;
  // temporary array  used to store f(a,b) = S_k' (u_i v_k - u_k v_i)^-(a+b) * (u_i v_k)^a (u_k v_i)^b factors that appear in the CF monopole spherical harmonic
  // indices are: [#u][#v][Power][particle#]
  Complex4DArray DerivativeFactors;
  
  // a duplicate array of DerivativeFactors used for precalculations with different organisation:
  // indices are: [particle#][#u][#v]
  Complex3DArray DerivativeFactors2;

  // temporary array used to store u spinor coordinates
  Complex* SpinorUCoordinates;
  // temporary array used to store v spinor coordinates
  Complex* SpinorVCoordinates;
  // temporary array used to store powers of the u spinor coordinates
  Complex** SpinorUCoordinatePower;
  // temporary array used to store powers of the v spinor coordinates
  Complex** SpinorVCoordinatePower;

  // Structure for calculating and storing the form of the derivatives
  // indices: [LL-index n][derivatives / u]
  SumDerivativeProduct** DerivativeStructure;

  // Matrix to store return values:  
  ComplexMatrix *Orbitals;

  // vectors used in Evaluate Orbitals for calculation of columns in Orbitals
  Complex *OrbitalVector;
  Complex *OrbitalVector2;
  
 public:

  // default constructor
  //
  JainCFOnSphereOrbitals();

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // jastrowPower = power to which the Jastrow factor has to be raised
  JainCFOnSphereOrbitals(int nbrParticles, int nbrLandauLevels, int nbrEffectiveFlux, int jastrowPower);
  
  // copy constructor
  //
  // function = reference on the wave function to copy
  JainCFOnSphereOrbitals(const JainCFOnSphereOrbitals& function);

  // destructor
  //
  ~JainCFOnSphereOrbitals();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = Matrix with CF Orbitals, 
  ComplexMatrix& operator ()(RealVector& x);


  // access for data members:
  int GetNbrParticles(){return NbrParticles;}

  int GetNbrOrbitals(){return NbrOrbitals;}

 protected:

  // evaluate precalculation tables used during wave function evaluation (called at each evaluation)
  //
  // x = point where the function has to be evaluated
  // derivativeFlag = indicate if precalculation tables invloved in derivative evaluation have to be calculated
  // return value = value of the Jastrow factor
  Complex EvaluateTables(RealVector& x, bool derivativeFlag = true);

  // evaluate normalization factors of projected monopole harmonics
  //
  void EvaluateNormalizationPrefactors();

  // evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
  //
  void EvaluateSumPrefactors();

  // evaluate Structure of Derivatives
  //
  void EvaluateDerivativeStructure();

  SumDerivativeProduct VRecursion(int VDerivatives);
  
  SumDerivativeProduct URecursion(int UDerivatives, int VDerivatives);

  // evaluates one line of CF orbitals and stores in internal matrix Orbitals
  // Index = line where to store results
  // momentum = monopole spherical harmonic Lz momentum
  // landauLevel = index of the pseudo Landau level
  // maximumMomentum = maxixum momentum that can be reached in the current pseudo Landau level
  void EvaluateOrbitals (int Index, int momentum, int landauLevel, int maximumMomentum);

};

#endif
