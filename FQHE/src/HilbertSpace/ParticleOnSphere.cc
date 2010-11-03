////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of particle on sphere                        //
//                                                                            //
//                        last modification : 15/07/2002                      //
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


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"

#include <iostream>


using std::cout;
using std::endl;


// virtual destructor
//

ParticleOnSphere::~ParticleOnSphere ()
{
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void ParticleOnSphere::SetTargetSpace(ParticleOnSphere* targetSpace)
{
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int ParticleOnSphere::GetTargetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id

int ParticleOnSphere::GetHilbertSpaceAdditionalSymmetry()
{
  return ParticleOnSphere::NoSymmetry;
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2), safe version i.e. works with any numbers of particles
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdAdAASafe (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->AdAdAA(index, m1, m2, n1, n2, coefficient);
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  double Coefficient = this->AA(index, n1, n2);
  int Index = this->AdAd(m1, m2, coefficient);
  coefficient *= Coefficient;
  return Index;
}

// apply a^+_m1 a^+_m2 a^+_m3 a_n1 a_n2 a_n3 operator to a given state (with m1+m2+m3=n1+n2+n3)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// m3 = third index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// n3 = third index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdAdAdAAA (int index, int m1, int m2, int m3, int n1, int n2, int n3, double& coefficient)
{
  int TmpM[3];
  int TmpN[3];
  TmpM[0] = m1;
  TmpM[1] = m2;
  TmpM[2] = m3;
  TmpN[0] = n1;
  TmpN[1] = n2;
  TmpN[2] = n3;
  return ProdAdProdA(index, TmpM, TmpN, 3, coefficient);
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1 a^+_m2 a^+_m3 a^+_m4 a_n1 a_n2 a_n3 a_n4 operator to a given state (with m1+m2+m3+m4=n1+n2+n3+n4)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// m3 = third index for creation operator
// m4 = fourth index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// n3 = third index for annihilation operator
// n4 = fourth index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdAdAdAdAAAA (int index, int m1, int m2, int m3, int m4, int n1, int n2, int n3, int n4, double& coefficient)
{
  int TmpM[4];
  int TmpN[4];
  TmpM[0] = m1;
  TmpM[1] = m2;
  TmpM[2] = m3;
  TmpM[3] = m4;
  TmpN[0] = n1;
  TmpN[1] = n2;
  TmpN[2] = n3;
  TmpN[3] = n4;
  return ProdAdProdA(index, TmpM, TmpN, 4, coefficient);
}

// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// m = Lz value of particle to be added
// coefficient = reference on the double where the multiplicative factor has to be stored
unsigned long ParticleOnSphere::Ad (unsigned long state, int m, double& coefficient)
{
  cout << "Attention: calling placeholder function ParticleOnSphereWithSpin::Ad - please override in inherited class!" <<endl;
  return 0x0l;
}


// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphere::AA (int index, int n1, int n2)
{
  return 0.0;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double ParticleOnSphere::ProdA (int index, int* n, int nbrIndices)
{
  return 0.0;
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdAd (int m1, int m2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}


// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphere::AdA (int index, int m)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdA (int index, int m, int n, double& coefficient)
{
  if (m != n)
    return this->HilbertSpaceDimension;
  else
    {
      coefficient = this->AdA(index, m);
      return index;
    }
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphere::AdA (long index, int m)
{
  return this->LargeHilbertSpaceDimension;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

long ParticleOnSphere::AdA (long index, int m, int n, double& coefficient)
{
  if (m != n)
    return this->HilbertSpaceDimension;
  else
    {
      coefficient = this->AdA(index, m);
      return index;
    }
}

// check whether HilbertSpace implements ordering of operators
//
bool ParticleOnSphere::HaveOrder ()
{
  return false;
}

// check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
int ParticleOnSphere::CheckOrder (int* m, int* n, int nbrIndices)
{
  return 0;
}


// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex ParticleOnSphere::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// return value = wave function evaluated at the given location

Complex ParticleOnSphere::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, 
						     this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location
Complex ParticleOnSphere::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
						int firstComponent, int nbrComponent)
{
  return Complex(0.0, 0.0);
}


// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex ParticleOnSphere::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, 
								 int nextCoordinates, int firstComponent, 
								 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void ParticleOnSphere::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
                                    
// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem

RealSymmetricMatrix ParticleOnSphere::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState)
{
  RealSymmetricMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem

HermitianMatrix ParticleOnSphere::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, ComplexVector& groundState)
{
  HermitianMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}
  
// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particle. The geometrical cut is a stripe.
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax+shitedCut to -Lzmax+shitedCut+subsytemSize-1)
// shiftedCut = first orbital belonging to the subsystem (with angular momentum -Lzmax+shitedCut)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphere::EvaluateShiftedPartialDensityMatrix (int subsytemSize, int nbrShiftedOrbitals, int nbrBosonSector, int lzSector, RealVector& groundState)
{
  if (nbrShiftedOrbitals == 0)
    return this->EvaluatePartialDensityMatrix(subsytemSize, nbrBosonSector, lzSector, groundState);
  RealSymmetricMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particle. The geometrical cut is a stripe.
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax+shitedCut to -Lzmax+shitedCut+subsytemSize-1)
// shiftedCut = first orbital belonging to the subsystem (with angular momentum -Lzmax+shitedCut)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix ParticleOnSphere::EvaluateShiftedPartialDensityMatrix (int subsytemSize, int nbrShiftedOrbitals, int nbrBosonSector, int lzSector, ComplexVector& groundState)
{
  if (nbrShiftedOrbitals == 0)
    return this->EvaluatePartialDensityMatrix(subsytemSize, nbrBosonSector, lzSector, groundState);
  HermitianMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphere::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState)
{
  RealSymmetricMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using real space partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// thetaTop =  inclination angle defining one edge of the cut in degrees
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphere::EvaluatePartialDensityMatrixRealSpacePartition (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, RealVector& groundState)
{
  RealSymmetricMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState)
{
  RealMatrix PartialEntanglementMatrix;
  return PartialEntanglementMatrix;  
}

// evaluate coeffecicents requested to compute the real space partition
//
// lzMax = twice the maximum angular momentum
// thetaTop = inclination angle defining one edge of the cut in radians
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in radians
// incompleteBetaThetaTop = reference on the pointer to the array (allocation done by the method) where the top part coefficients will be stored
// incompleteBetaThetaBotton = reference on the pointer to the array (allocation done by the method) where the bottom part coefficients will be stored

void ParticleOnSphere::EvaluatePartialDensityMatrixRealSpacePartitionCoefficient(int lzMax, double thetaTop, double thetaBottom, double*& incompleteBetaThetaTop, double*& incompleteBetaThetaBottom)
{
  double * LogFactorialsNphi = new double[lzMax +2];
  LogFactorialsNphi[0] = 0.0;
  LogFactorialsNphi[1] = 0.0;
  for (int i = 2 ; i < lzMax+2; ++i)
    LogFactorialsNphi[i] = LogFactorialsNphi[i - 1] + log((double) i); 
  
  
  double LogSinThetaTop = 0;
  
  if (thetaTop > 1e-10)
    LogSinThetaTop = 2.0*log(sin(thetaTop/2.));
  
  double LogCosThetaTop = 2.0*log(cos(thetaTop/2.));
  
  double LogSinThetaBot = 2.0*log(sin(thetaBottom/2.));
  double LogCosThetaBot = 2.0*log(cos(thetaBottom/2.));
  
  // Compute the incomplete beta function for x=sin^2(thetaTop/2) and x=sin^2(thetaBottom/2)
  incompleteBetaThetaTop = new double[lzMax + 1];
  incompleteBetaThetaBottom = new double[lzMax + 1];
  if ( thetaTop < 1e-10)
    incompleteBetaThetaTop [lzMax] = 0.0;
  else
    incompleteBetaThetaTop [lzMax] = exp((lzMax+1)*LogSinThetaTop);
  incompleteBetaThetaBottom [lzMax] =  exp((lzMax+1)*LogSinThetaBot);
  
  for (int i=lzMax-1 ; i>= 0; --i)
    {
      if (thetaTop < 1e-10)
	incompleteBetaThetaTop[i] = 0.0;
      else
	incompleteBetaThetaTop[i] = exp(LogFactorialsNphi[lzMax+1] - LogFactorialsNphi[i+1] - LogFactorialsNphi[lzMax-i] + (i+1.0)*LogSinThetaTop + (lzMax - i)*LogCosThetaTop) + incompleteBetaThetaTop[i+1] ;
      
      incompleteBetaThetaBottom[i] = exp(LogFactorialsNphi[lzMax+1] - LogFactorialsNphi[i+1] - LogFactorialsNphi[lzMax-i] + (i+1.0)*LogSinThetaBot + (lzMax - i)*LogCosThetaBot) + incompleteBetaThetaBottom[i+1] ;
    }
  delete [] LogFactorialsNphi;
}


// compute part of the Schmidt decomposition, allowing cut in the reduced denisty matrix eigenvalue space
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// eigenvalueCut = discard all contribution from the reduced density matrix whose eigenvalues is lower than eigenvalueCut
// rebuiltSchmidtGroundState = reference on the state to whose current sector contribution to the Schmidt decomposition will be added 
// diagonalizedDensityMatrix = reference on the diagonalized reduced density matrix associated to the current sector (with down ordered diagonal elements)
// transformationMatrix =  reference on the transformation matric that diagonalizes the reduced density matrix
// return value = reference on rebuiltSchmidtGroundState

RealVector& ParticleOnSphere::EvaluatePartialSchmidtDecomposition(int subsytemSize, int nbrFermionSector, int lzSector, double eigenvalueCut,
								  RealVector& groundState, RealVector& rebuiltSchmidtGroundState,
								  RealDiagonalMatrix& diagonalizedDensityMatrix, RealMatrix& transformationMatrix)
{
  return rebuiltSchmidtGroundState;
}

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int ParticleOnSphere::FindStateIndex(char* stateDescription)
{
  return -1;
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component has to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& ParticleOnSphere::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& ParticleOnSphere::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  return state;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& ParticleOnSphere::PrintStateMonomial (ostream& Str, int state)
{
  return Str;
}

// fuse two states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// leftVector = reference on the vector whose Hilbert space will be fuse to the left
// rightVector = reference on the vector whose Hilbert space will be fuse to the right
// padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
// leftSpace = point to the Hilbert space that will be fuse to the left
// rightSpace = point to the Hilbert space that will be fuse to the right
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

RealVector& ParticleOnSphere::FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
					  ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace, bool symmetrizedFlag, double coefficient)
{
  return outputVector;
}

// use product rule to produce part of the components of a system from a smaller one
//
// outputVector = reference on the vector which will contain the product rule state  (without zeroing components which do not occur in the fusion)
// inputVector = reference on the vector associated to the smaller system
// inputSpace = pointer to the Hilbert space of the smaller system
// commonPattern = array describing the shared leftmost pattern between the n-body states in both the smaller and larger system sizes
// commonPatterSize = number of elements in the commonPattern array
// addedPattern = array describing the pattern that has to be inserted to go from the smaller system to the larger one
// addedPatterSize = number of elements in the addedPattern array
// coefficient = multiplicqtive fqctor to go fron the component of the smaller system to the larger one
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// return value = reference on the product rule state

RealVector& ParticleOnSphere::ProductRules (RealVector& outputVector, RealVector& inputVector, ParticleOnSphere* inputSpace, 
					    int* commonPattern, int commonPatterSize, int* addedPattern, int addedPatterSize,
					    double coefficient, bool symmetrizedFlag)
{
  return outputVector;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double ParticleOnSphere::JackSqrNormalization (RealVector& outputVector, long minIndex, long nbrComponents)
{
  return 0.0;
}

// remove part of each Fock state, discarding component if the Fock state does not a given pattern
//
// inputVector = state to truncate
// reducedSpace = Hilbert space where the truncated state will lie
// pattern = array describing the pattern 
// patternSize = pattern size
// patternShift = indicate where the pattern has to be applied
// return value = trucated state

RealVector ParticleOnSphere::TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnSphere* reducedSpace, int* pattern, int patternSize, int patternShift)
{
  RealVector TmpVector;
  return TmpVector;
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals
bool ParticleOnSphere::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  return false;
}


// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the Lz component
int ParticleOnSphere::GetLzValue(int j)
{
  return 0;
}

// get Sz component of the spin
//
// j = index of the vector in Hilbert space
// return value = Sz component

int ParticleOnSphere::GetSzValue(int j)
{
  return 0;
} 

// transform a vector belonging to this vector space in the lz->-lz
//
// finalSpace = the space obtained after the lz->-lz operation
// initialVector = vector on which the operation will be apply
// return value = vector resulting of the operation

RealVector ParticleOnSphere::GetLzSymmetricVector(ParticleOnSphere* finalSpace, RealVector& initialVector)
{
  RealVector Tmp (this->LargeHilbertSpaceDimension, true);
  return Tmp;
}
