////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 19/10/2004                        //
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
#include "Complex.h"
#include "Vector/ComplexVector.h"
#include "Hamiltonian/QuantumDotHamiltonian/PeriodicElectronHole3DHamiltonian.h"
#include "Tools/QuantumDot/Potential/ThreeDConstantCellPotential.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PeriodicThreeDOneParticle.h"

#include <iostream>
#include <math.h>

using std::ostream;
using std::cout;
using std::endl;

#define PERIODIC_HAMILTONIAN_FACTOR 150.4
#define COULOMBIAN_FACTOR 180.79
/*
// each state is coded with 5 bits
#define NbrBitZ2 5
#define NbrBitY2 5
#define NbrBitX2 5
#define NbrBitZ1 5
#define NbrBitY1 5
#define NbrBitX1 5

#define Hex2 0x7fff 
#define Hex1 0x7fff
*/

// each state is coded with 3 bits
#define NbrBitZ2 3
#define NbrBitY2 3
#define NbrBitX2 3
#define NbrBitZ1 3
#define NbrBitY1 3
#define NbrBitX1 3

#define Hex2 0x1ff
#define Hex1 0x1ff

#define NbrBit2 (NbrBitX2 + NbrBitY2 + NbrBitZ2)
#define NbrBit1 (NbrBitX1 + NbrBitY1 + NbrBitZ1)

#define ShiftZ2 0
#define ShiftY2 (ShiftZ2 + NbrBitZ2)
#define ShiftX2 (ShiftY2 + NbrBitY2)
#define ShiftZ1 (ShiftX2 + NbrBitX2)
#define ShiftY1 (ShiftZ1 + NbrBitZ1)
#define ShiftX1 (ShiftY1 + NbrBitY1)
#define ShiftZ1bis 0
#define ShiftY1bis (ShiftZ1bis + NbrBitZ1)
#define ShiftX1bis (ShiftY1bis + NbrBitZ2)


// constructor
//
// space = pointer to the Hilbert space of two particles
// Mex, Mey, Mez = effective masses in three directions of electron (in vacuum electron mass unit)
// Mhx, Mhy, Mhz = effective masses in three directions of hole (in vacuum electron mass unit)
// potentialElectron = pointer to the potential for electron
// potentialHole = pointer to the potential for hole
// xSize, ySize, zSize = sizes of the sample in three direction (in Angstrom unit)
// dielectric = dielectric constant in the sample

PeriodicElectronHole3DHamiltonian::PeriodicElectronHole3DHamiltonian (PeriodicThreeDTwoParticles* space, double Mex, double Mey, double Mez, double Mhx, double Mhy, double Mhz, ThreeDConstantCellPotential* potentialElectron, ThreeDConstantCellPotential* potentialHole, double xSize, double ySize, double zSize, double dielectric)
{
  this->Space = space;
  PeriodicThreeDOneParticle* firstParticle = (PeriodicThreeDOneParticle*) this->Space->GetFirstParticleSpace ();
  PeriodicThreeDOneParticle* secondParticle = (PeriodicThreeDOneParticle*) this->Space->GetSecondParticleSpace ();
  this->NbrState1X = firstParticle->GetNbrStateX ();
  this->NbrState1Y = firstParticle->GetNbrStateY ();
  this->NbrState1Z = firstParticle->GetNbrStateZ ();
  this->NbrState2X = secondParticle->GetNbrStateX ();
  this->NbrState2Y = secondParticle->GetNbrStateY ();
  this->NbrState2Z = secondParticle->GetNbrStateZ ();

  if (((this->NbrState1X * 2 - 1) > pow(2, NbrBitX1)) || ((this->NbrState1Y * 2 - 1) > pow(2, NbrBitY1)) || ((this->NbrState1Z * 2 - 1) > pow(2, NbrBitZ1)) || ((this->NbrState2X * 2 - 1) > pow(2, NbrBitX2)) || ((this->NbrState2Y * 2 - 1) > pow(2, NbrBitY2)) || ((this->NbrState2Z * 2 - 1) > pow(2, NbrBitZ2)))
    {
      cout << "At least a number of states in a direction is too big.!" << endl;
      exit (1);
    }
  this->MakeConversionTable ();
  cout << "Evaluating the kinetic terms ..." << endl;
  this->EvaluateKineticTerm (Mex, Mey, Mez, Mhx, Mhy, Mhz, firstParticle, secondParticle, xSize, ySize, zSize);
  cout << "Evaluating the electron confinement terms ..." << endl;
  this->EvaluateConfinementTerm (potentialElectron, firstParticle, 1, this->RealElectronConfinement, this->ImaginaryElectronConfinement);
  cout << "Evaluating the hole confinement terms ..." << endl;
  this->EvaluateConfinementTerm (potentialHole, secondParticle, 2, this->RealHoleConfinement, this->ImaginaryHoleConfinement);
  cout << "Evaluating the Coulombian term ..." << endl;
  this->EvaluateCoulombianTerm (xSize, ySize, zSize, dielectric);
  delete firstParticle; 
  delete secondParticle;
  cout << "End of the evaluation." << endl;
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

PeriodicElectronHole3DHamiltonian::PeriodicElectronHole3DHamiltonian(const PeriodicElectronHole3DHamiltonian& hamiltonian)
{
  this->Space = hamiltonian.Space;
  this->NbrState1X = NbrState1X;
  this->NbrState1Y = NbrState1Y;
  this->NbrState1Z = NbrState1Z;
  this->NbrState2X = NbrState2X;
  this->NbrState2Y = NbrState2Y;
  this->NbrState2Z = NbrState2Z;
  this->IToX = hamiltonian.IToX;
  this->IToX1 = hamiltonian.IToX1;
  this->IToX2 = hamiltonian.IToX2;
  this->KineticTerm = hamiltonian.KineticTerm;
  this->RealElectronConfinement = hamiltonian.RealElectronConfinement;
  this->ImaginaryElectronConfinement = hamiltonian.ImaginaryElectronConfinement;
  this->RealHoleConfinement = hamiltonian.RealHoleConfinement;
  this->ImaginaryHoleConfinement = hamiltonian.ImaginaryHoleConfinement;
  this->CoulombianTerm = hamiltonian.CoulombianTerm;
}

// destructor
//

PeriodicElectronHole3DHamiltonian::~ PeriodicElectronHole3DHamiltonian()
{
  cout << "PeriodicElectronHole3DHamiltonian destructor is being called." << endl;
  delete   this->Space;
  cout << "Destructor of Hilbert space is being called" << endl;
  delete[] this->IToX;
  delete[] this->IToX1;
  delete[] this->IToX2;
  delete[] this->KineticTerm;
  delete[] this->RealElectronConfinement;
  delete[] this->ImaginaryElectronConfinement;
  delete[] this->RealHoleConfinement;
  delete[] this->ImaginaryHoleConfinement;
  delete[] this->CoulombianTerm;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* PeriodicElectronHole3DHamiltonian::Clone ()
{
  return new PeriodicElectronHole3DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PeriodicElectronHole3DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void PeriodicElectronHole3DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->KineticTerm[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicElectronHole3DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
{
  double x = 0.0;
  int dim = this->Space->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicElectronHole3DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& PeriodicElectronHole3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Space->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of idinces 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& PeriodicElectronHole3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
						       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination.Re(i) = 0.0;
      vDestination.Im(i) = 0.0;
    }
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& PeriodicElectronHole3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Space->GetHilbertSpaceDimension());
}
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& PeriodicElectronHole3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{ 
  int dimension = this->Space->GetHilbertSpaceDimension ();
  int lastComponent = firstComponent + nbrComponent;
  int XY1 = 0, XY2 = 0; int Sum1 = 0;
  int X1 = 0, Y1 = 0, X2 = 0, Y2 = 0;
  double* tmpCoulombian; double* tmpRealElectron; double* tmpImaginaryElectron;
  double* tmpRealHole; double* tmpImaginaryHole;

  for (int Index1 = firstComponent; Index1 < lastComponent; ++Index1)
    {
      XY1 = this->IToX[Index1];
      Y1 = XY1 & Hex2; X1 = (XY1 >> NbrBit2) & Hex1;
      Sum1 = X1 + Y1;
      // kinetic term
      vDestination.Re(Index1) += this->KineticTerm[Index1] * vSource.Re(Index1);
      vDestination.Im(Index1) += this->KineticTerm[Index1] * vSource.Im(Index1);
      
      tmpCoulombian = this->CoulombianTerm[X1];
      tmpRealElectron = this->RealElectronConfinement[X1]; tmpImaginaryElectron = this->ImaginaryElectronConfinement[X1];
      tmpRealHole = this->RealHoleConfinement[Y1]; tmpImaginaryHole = this->ImaginaryHoleConfinement[Y1];
      for (int Index2 = 0; Index2 < dimension; ++Index2)
	{
	  XY2 = this->IToX[Index2];
	  Y2 = XY2 & Hex2; X2 = (XY2 >> NbrBit2) & Hex1;
	  
	  // Coulombian term	  
	  if ((X2 + Y2) == Sum1)
	    {
	      // cout << "Coulomb: " << Index1 << " " << Index2 << " " << XY1 << " " << XY2 << " " << X1 << " " << Y1 << " " << X2 << " " << Y2 << endl;
	      vDestination.Re(Index1) += tmpCoulombian[X2] * vSource.Re(Index2);
	      vDestination.Im(Index1) += tmpCoulombian[X2] * vSource.Im(Index2);
	    }
	  	  
	  // confinement term for the electrons	  
	  if (Y1 == Y2)
	    {
	      //cout << "Electron: " << Index1 << " " << Index2 << " " << XY1 << " " << XY2 << " " << X1 << " " << Y1 << " " << X2 << " " << Y2 << endl;
	      vDestination.Re(Index1) += (tmpRealElectron[X2] * vSource.Re(Index2) - tmpImaginaryElectron[X2] * vSource.Im(Index2));
	      vDestination.Im(Index1) += (tmpRealElectron[X2] * vSource.Im(Index2) + tmpImaginaryElectron[X2] * vSource.Re(Index2));
	    }
	  	  
	  // confinement term for the holes
	  if (X1 == X2)
	    {
	      //cout << "Hole: " << Index1 << " " << Index2 << " " << XY1 << " " << XY2 << " " << X1 << " " << Y1 << " " << X2 << " " << Y2 << endl;
	      vDestination.Re(Index1) += (tmpRealHole[Y2] * vSource.Re(Index2) - tmpImaginaryHole[Y2] * vSource.Im(Index2));
	      vDestination.Im(Index1) += (tmpRealHole[Y2] * vSource.Im(Index2) + tmpImaginaryHole[Y2] * vSource.Re(Index2));
	    }	  
	}
    }
  return vDestination;
}
 
// determine the maximal value of the kenetic elements
//
// return = the wanted value

double PeriodicElectronHole3DHamiltonian::MaxKineticElement()
{
  double tmp = this->KineticTerm[0];
  for (int i = 1; i < this->Space->GetHilbertSpaceDimension(); ++i)
    if (tmp < this->KineticTerm[i])
      tmp = this->KineticTerm[i];
  return tmp;
}

// make the conversion table to hexadecimal indices
//

void PeriodicElectronHole3DHamiltonian::MakeConversionTable ()
{
  int dimension = this->Space->GetHilbertSpaceDimension ();
  
  this->IToX = new int [dimension];
  int index = 0;
  //cout << "I to X:" << endl;
  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    for (int n1 = 0; n1 < this->NbrState1Y; ++n1)      
      for (int p1 = 0; p1 < this->NbrState1Z; ++p1)	    	    
	for (int m2 = 0; m2 < this->NbrState2X; ++m2)
	  for (int n2 = 0; n2 < this->NbrState2Y; ++n2)      
	    for (int p2 = 0; p2 < this->NbrState2Z; ++p2)	    
	      {
		this->IToX[index] = (m1 << ShiftX1) | (n1 << ShiftY1) | (p1 << ShiftZ1) | (m2 << ShiftX2) | (n2 << ShiftY2) | (p2 << ShiftZ2);
		//cout << "Index: " << index << " " << m1 << " " << n1 << " " << p1 << " " << m2 << " " << n2 << " " << p2 << " " << this->IToX[index] << endl;
		++index;
	      }

  this->IToX1 = new int** [this->NbrState1X];
  //cout << "I to X 1 :" << endl;
  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    {
      this->IToX1[m1] = new int* [this->NbrState1Y];
      for (int n1 = 0; n1 < this->NbrState1Y; ++n1)   
	{
	  this->IToX1[m1][n1] = new int [this->NbrState1Z];
	  for (int p1 = 0; p1 < this->NbrState1Z; ++p1)	    
	    {
	      this->IToX1[m1][n1][p1] = (m1 << ShiftX1bis) | (n1 << ShiftY1bis) | (p1 << ShiftZ1bis);
	      //cout <<  m1 << " " << n1 << " " << p1 << " " << this->IToX1[m1][n1][p1] << endl;	  
	    }
	}
    }

  this->IToX2 = new int** [this->NbrState2X];
  for (int m2 = 0; m2 < this->NbrState2X; ++m2)
    {
      this->IToX2[m2] = new int* [this->NbrState2Y];
      for (int n2 = 0; n2 < this->NbrState2Y; ++n2)   
	{
	  this->IToX2[m2][n2] = new int [this->NbrState2Z];
	  for (int p2 = 0; p2 < this->NbrState2Z; ++p2)
	    {
	      this->IToX2[m2][n2][p2] = (m2 << ShiftX2) | (n2 << ShiftY2) | (p2 << ShiftZ2);
	      //cout <<  m2 << " " << n2 << " " << p2 << " " << this->IToX2[m2][n2][p2] << endl;	  
	    }
	}
    }

}

// evaluate the kinetic term
//
// mex, mey, mez = effective masses in three directions of electron (in vacuum electron mass unit)
// mhx, mhy, mhz = effective masses in three directions of hole (in vacuum electron mass unit)
// firstParticle = pointer to the first particle's Hilbert space
// secondParticle = pointer to the second particle's Hilbert space
// xSize, ySize, zSize = sizes of the sample in three direction (in Angstrom unit)

void PeriodicElectronHole3DHamiltonian::EvaluateKineticTerm (double mex, double mey, double mez, double mhx, double mhy, double mhz, PeriodicThreeDOneParticle* firstParticle, PeriodicThreeDOneParticle* secondParticle, double xSize, double ySize, double zSize)
{
  int dimension = this->Space->GetHilbertSpaceDimension ();
  this->KineticTerm = new double [dimension];
  int LowerImpulsion1X = firstParticle->GetLowerImpulsionX();
  int LowerImpulsion1Y = firstParticle->GetLowerImpulsionY();
  int LowerImpulsion1Z = firstParticle->GetLowerImpulsionZ();
  int LowerImpulsion2X = secondParticle->GetLowerImpulsionX();
  int LowerImpulsion2Y = secondParticle->GetLowerImpulsionY();
  int LowerImpulsion2Z = secondParticle->GetLowerImpulsionZ();

  //cout << "LowerImpulsion: " << LowerImpulsion1X << " " << LowerImpulsion1Y << " " << LowerImpulsion1Z << " " << LowerImpulsion2X << " " << LowerImpulsion2Y << " " << LowerImpulsion2Z << endl;
  
  double InvXFactor1 = PERIODIC_HAMILTONIAN_FACTOR / (mex * xSize * xSize);
  double InvYFactor1 = PERIODIC_HAMILTONIAN_FACTOR / (mey * ySize * ySize);
  double InvZFactor1 = PERIODIC_HAMILTONIAN_FACTOR / (mez * zSize * zSize);
  
  double InvXFactor2 = PERIODIC_HAMILTONIAN_FACTOR / (mhx * xSize * xSize);
  double InvYFactor2 = PERIODIC_HAMILTONIAN_FACTOR / (mhy * ySize * ySize);
  double InvZFactor2 = PERIODIC_HAMILTONIAN_FACTOR / (mhz * zSize * zSize);
  
  double* factorX1 = new double [this->NbrState1X]; double* factorY1 = new double [this->NbrState1Y]; double* factorZ1 = new double [this->NbrState1Z];
  double* factorX2 = new double [this->NbrState2X]; double* factorY2 = new double [this->NbrState2Y]; double* factorZ2 = new double [this->NbrState2Z];
  
  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    factorX1[m1] = InvXFactor1 * ((double) (m1 + LowerImpulsion1X)) * ((double) (m1 + LowerImpulsion1X));
  for (int n1 = 0; n1 < this->NbrState1Y; ++n1)
    factorY1[n1] = InvYFactor1 * ((double) (n1 + LowerImpulsion1Y)) * ((double) (n1 + LowerImpulsion1Y));
  for (int p1 = 0; p1 < this->NbrState1Z; ++p1)
    factorZ1[p1] = InvZFactor1 * ((double) (p1 + LowerImpulsion1Z)) * ((double) (p1 + LowerImpulsion1Z));
  for (int m2 = 0; m2 < this->NbrState2X; ++m2)
    factorX2[m2] = InvXFactor2 * ((double) (m2 + LowerImpulsion2X)) * ((double) (m2 + LowerImpulsion2X));
  for (int n2 = 0; n2 < this->NbrState2Y; ++n2)
    factorY2[n2] = InvYFactor2 * ((double) (n2 + LowerImpulsion2Y)) * ((double) (n2 + LowerImpulsion2Y));
  for (int p2 = 0; p2 < this->NbrState2Z; ++p2)
    factorZ2[p2] = InvZFactor2 * ((double) (p2 + LowerImpulsion2Z)) * ((double) (p2 + LowerImpulsion2Z));

  //cout << "Kinetic terms: " << endl;
  int index = 0; double FactorX1, FactorY1, FactorZ1, FactorX2, FactorY2;
  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    {
      FactorX1 = factorX1[m1];
      for (int n1 = 0; n1 < this->NbrState1Y; ++n1)	
	{
	  FactorY1 = FactorX1 + factorY1[n1];
	  for (int p1 = 0; p1 < this->NbrState1Z; ++p1)
	    {
	      FactorZ1 = FactorY1 + factorZ1[p1];	      
	      for (int m2 = 0; m2 < this->NbrState2X; ++m2)
		{
		  FactorX2 = FactorZ1 + factorX2[m2];
		  for (int n2 = 0; n2 < this->NbrState2Y; ++n2)
		    {
		      FactorY2 = FactorX2 + factorY2[n2];	  
		      for (int p2 = 0; p2 < this->NbrState2Z; ++p2)
			{
			  this->KineticTerm[index] = FactorY2 + factorZ2[p2];	 
			  //cout << index << " " <<  m1 << " " << n1 << " " << p1 << " " << m2 << " " << n2 << " " << p2 << " " << this->KineticTerm[index] << endl;
			  ++index;
			}
		    }
		}
	    }
	}
    }

  delete[] factorX1; delete[] factorY1; delete[] factorZ1;
  delete[] factorX2; delete[] factorY2; delete[] factorZ2;
}

// evaluate the confinement terms for electrons and holes
//
// potential = pointer to the potential for the considered carrier
// particle = pointer to the Hilbertspace for the considered carrier
// type = type of the carrier, 1 for electron, 2 for hole
// realConfinement = reference to 2D array of real elements of the wanted terms
// imaginaryConfinement = reference to 2D array of imaginary elements of the wanted terms

void PeriodicElectronHole3DHamiltonian::EvaluateConfinementTerm (ThreeDConstantCellPotential* potential, PeriodicThreeDOneParticle* particle, int type, double** &realConfinement, double** &imaginaryConfinement)
{
  int NbrStateX = particle->GetNbrStateX (), NbrStateY = particle->GetNbrStateY (), NbrStateZ = particle->GetNbrStateZ ();
  int NbrCellX = potential->GetNbrCellX (), NbrCellY = potential->GetNbrCellY (), NbrCellZ = potential->GetNbrCellZ ();

  double** RealWaveFunctionOverlapX; double** RealWaveFunctionOverlapY; double** RealWaveFunctionOverlapZ; 
  double** ImaginaryWaveFunctionOverlapX; double** ImaginaryWaveFunctionOverlapY; double** ImaginaryWaveFunctionOverlapZ; 

  this->EvaluateWaveFunctionOverlap (NbrCellX, NbrStateX, RealWaveFunctionOverlapX, ImaginaryWaveFunctionOverlapX);
  this->EvaluateWaveFunctionOverlap (NbrCellY, NbrStateY, RealWaveFunctionOverlapY, ImaginaryWaveFunctionOverlapY);
  this->EvaluateWaveFunctionOverlap (NbrCellZ, NbrStateZ, RealWaveFunctionOverlapZ, ImaginaryWaveFunctionOverlapZ);

  int LengthX = (NbrStateX - 1) * 2 + 1; int LengthY = (NbrStateY - 1) * 2 + 1; int LengthZ = (NbrStateZ - 1) * 2 + 1;

  double*** TmpReal = new double** [LengthX];
  double*** TmpImaginary = new double** [LengthX];

  double TmpRe, TmpIm;
  double TmpRe2, TmpIm2;
  double* TmpRealWaveFunctionOverlapX;
  double* TmpImaginaryWaveFunctionOverlapX;
  double* TmpRealWaveFunctionOverlapY;
  double* TmpImaginaryWaveFunctionOverlapY;
  double* TmpRealPrecalculatedHamiltonian;
  double* TmpImaginaryPrecalculatedHamiltonian;

  for (int m = 0; m < LengthX; ++m)
    {
      TmpReal[m] = new double* [LengthY];
      TmpImaginary[m] = new double* [LengthY];
      TmpRealWaveFunctionOverlapX = RealWaveFunctionOverlapX[m];
      TmpImaginaryWaveFunctionOverlapX = ImaginaryWaveFunctionOverlapX[m];	      	  
      for (int n = 0; n < LengthY; ++n)
	{	  
	  TmpReal[m][n] = new double [NbrCellZ];
	  TmpImaginary[m][n] = new double [NbrCellZ];
	  TmpRealWaveFunctionOverlapY = RealWaveFunctionOverlapY[n];
	  TmpImaginaryWaveFunctionOverlapY = ImaginaryWaveFunctionOverlapY[n];	  
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];		  
	  for (int CellZ = 0; CellZ < NbrCellZ; ++CellZ)
	    {
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellY = 0; CellY < NbrCellY; ++CellY)
		{
		  TmpRe2 = TmpRealWaveFunctionOverlapY[CellY];
		  TmpIm2 = TmpImaginaryWaveFunctionOverlapY[CellY];
		  for (int CellX = 0; CellX < NbrCellX; ++CellX)
		    {		      
		      TmpRe += potential->GetPotential (CellX, CellY, CellZ) * (TmpRealWaveFunctionOverlapX[CellX] * TmpRe2 - TmpImaginaryWaveFunctionOverlapX[CellX] * TmpIm2);
		      TmpIm += potential->GetPotential (CellX, CellY, CellZ) * (TmpRealWaveFunctionOverlapX[CellX] * TmpIm2 + TmpImaginaryWaveFunctionOverlapX[CellX] * TmpRe2);		      
		    }
		}
	      TmpRealPrecalculatedHamiltonian[CellZ] = TmpRe;  
	      TmpImaginaryPrecalculatedHamiltonian[CellZ] = TmpIm;  
	    }
	}
    }

  double*** RealPrecalculatedHamiltonian = new double** [LengthX];
  double*** ImaginaryPrecalculatedHamiltonian = new double** [LengthX];

  double* TmpRealWaveFunctionOverlapZ;
  double* TmpImaginaryWaveFunctionOverlapZ;
  for (int m = 0; m < LengthX; ++m)
    {
      RealPrecalculatedHamiltonian[m] = new double* [LengthY];      
      ImaginaryPrecalculatedHamiltonian[m] = new double* [LengthY]; 
      for (int n = 0; n < LengthY; ++n)
	{
	  RealPrecalculatedHamiltonian[m][n] = new double [LengthZ];      
	  ImaginaryPrecalculatedHamiltonian[m][n] = new double [LengthZ]; 
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];
	  for (int p = 0; p < LengthZ; ++p)
	    {
	      TmpRealWaveFunctionOverlapZ = RealWaveFunctionOverlapZ[p];
	      TmpImaginaryWaveFunctionOverlapZ = ImaginaryWaveFunctionOverlapZ[p];
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellZ = 0; CellZ < NbrCellZ; ++CellZ)
		{
		  TmpRe += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ] - TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ]);
		  TmpIm += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ] + TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ]);
		}
	      RealPrecalculatedHamiltonian[m][n][p] = TmpRe;
	      ImaginaryPrecalculatedHamiltonian[m][n][p] = TmpIm;
	    }
	}
    }
  
  int nbrBit = 0; int*** iToX;
  if (type == 1)
    {
      nbrBit = NbrBit1;
      iToX = this->IToX1;
    }
  else    
    if (type == 2)
      {
	nbrBit = NbrBit2;
	iToX = this->IToX2;
      }
    else
      {
	cout << "This type of particle is not taken into account. Exit now!" << endl;
	exit (1);
      }

  int number = (1 << nbrBit);
  realConfinement = new double* [number];
  imaginaryConfinement = new double* [number];
  for (int i1 = 0; i1 < number; ++i1)
    {
      realConfinement[i1] = new double [number];
      imaginaryConfinement[i1] = new double [number];
      for (int i2 = 0; i2 < number; ++i2)
	{
	  realConfinement[i1][i2] = 0.0;
	  imaginaryConfinement[i1][i2] = 0.0;
	}
    }

  int OriginX = NbrStateX - 1, OriginY = NbrStateY - 1, OriginZ = NbrStateZ - 1;
  double* real; double* imaginary; int tmpIndex;
  for (int m1 = 0; m1 < NbrStateX; ++m1)
    for (int n1 = 0; n1 < NbrStateY; ++n1)
      for (int p1 = 0; p1 < NbrStateZ; ++p1)
	{
	  int X1 = iToX[m1][n1][p1];	  
	  for (int m3 = 0; m3 < NbrStateX; ++m3)
	    for (int n3 = 0; n3 < NbrStateY; ++n3)
	      {
		real = RealPrecalculatedHamiltonian[-m1 + m3 + OriginX][-n1 + n3 + OriginY];
		imaginary = ImaginaryPrecalculatedHamiltonian[-m1 + m3 + OriginX][-n1 + n3 + OriginY];
		tmpIndex = OriginZ - p1;
		for (int p3 = 0; p3 < NbrStateZ; ++p3)	
		  {
		    int X2 = iToX[m3][n3][p3];
		    realConfinement[X1][X2] = real[tmpIndex];
		    imaginaryConfinement[X1][X2] = imaginary[tmpIndex];		      	       
		    ++tmpIndex;
		    //cout << m1 << " " << n1 << " " << p1 << " " << m3 << " " << n3 << " " << p3 << " " << realConfinement[X1][X2] << " " << imaginaryConfinement[X1][X2] << endl;
		  }
	      }
	}

  delete[] RealWaveFunctionOverlapX; delete[] RealWaveFunctionOverlapY; delete[] RealWaveFunctionOverlapZ;
  delete[] ImaginaryWaveFunctionOverlapX; delete[] ImaginaryWaveFunctionOverlapY; delete[] ImaginaryWaveFunctionOverlapZ;
  delete[] TmpReal; delete[] TmpImaginary;
  delete[] RealPrecalculatedHamiltonian; delete[] ImaginaryPrecalculatedHamiltonian; 
}

// evaluate the Coulombian term
//
// xSize, ySize, zSize = sizes of the sample in three direction (in Angstrom unit)
// dielectric = dielectric constant in the sample

void PeriodicElectronHole3DHamiltonian::EvaluateCoulombianTerm (double xSize, double ySize, double zSize, double dielectric)
{
  double factor = COULOMBIAN_FACTOR / dielectric;
  double squareX = xSize * xSize;
  double squareY = ySize * ySize;
  double squareZ = zSize * zSize;
  double volume = xSize * ySize * zSize;

  int number = (1 << NbrBit1);
  this->CoulombianTerm = new double* [number];
  for (int i1 = 0; i1 < number; ++i1)
    this->CoulombianTerm[i1] = new double [number];

  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    for (int n1 = 0; n1 < this->NbrState1Y; ++n1)      
      for (int p1 = 0; p1 < this->NbrState1Z; ++p1)
	{
	  int X1 = this->IToX1[m1][n1][p1];	  
	  for (int m3 = 0; m3 < this->NbrState1X; ++m3)
	    for (int n3 = 0; n3 < this->NbrState1Y; ++n3)      
	      for (int p3 = 0; p3 < this->NbrState1Z; ++p3)
		{
		  int X2 = this->IToX1[m3][n3][p3];
		  double tmp = (m1 - m3) * (m1 - m3) / squareX + (n1 - n3) * (n1 - n3) / squareY + (p1 - p3) * (p1 - p3) / squareZ;
		  if ((m1 == m3) && (n1 == n3) && (p1 ==p3))       	       
		    this->CoulombianTerm[X1][X2] = 0;
		  else
		    this->CoulombianTerm[X1][X2] = factor / (tmp * volume);
		}
	}
}

// evaluate the wave function overlap
//
// nbrStep = number of steps in the given direction
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool PeriodicElectronHole3DHamiltonian::EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray)
{
  double Diff = 0.0;
  double Tmp = 0.0;
  double Tmp1 = 1.0 / double (nbrStep);
  int Length = (nbrState - 1) * 2 + 1;
  realArray = new double* [Length];
  imaginaryArray = new double* [Length];  
  int Origin = nbrState - 1;
  for (int delta = 0; delta < Origin; ++delta)
    {
      realArray[delta] = new double [nbrStep];
      imaginaryArray[delta] = new double [nbrStep];
      Diff = 2.0 * M_PI * double (delta - Origin);
      Tmp = Diff / nbrStep;	
      Diff = 1.0 / Diff;	
      for (int i = 0; i < nbrStep; ++i)
	{
	  realArray[delta][i] = Diff * (sin(Tmp * (i + 1)) - sin(Tmp * i));
	  imaginaryArray[delta][i] = Diff * (cos(Tmp * (i + 1)) - cos(Tmp * i));
	}    
    } 
  realArray[Origin] = new double [nbrStep];
  imaginaryArray[Origin] = new double [nbrStep];
  for (int i = 0; i < nbrStep; ++i)
    {
      realArray[Origin][i] = Tmp1;
      imaginaryArray[Origin][i] = 0.0;
    }
  for (int delta = Origin + 1; delta < Length; ++delta)
    {
      realArray[delta] = new double [nbrStep];
      imaginaryArray[delta] = new double [nbrStep];      
      for (int i = 0; i < nbrStep; ++i)
	{
	  realArray[delta][i] = realArray[Length - 1 - delta][i];
	  imaginaryArray[delta][i] = -imaginaryArray[Length - 1 - delta][i];
	}    	
    }
  return true;
}
