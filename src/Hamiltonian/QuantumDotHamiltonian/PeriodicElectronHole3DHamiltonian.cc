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
#define COULOMBIAN_FACTOR 100 // factor of multiplication for coulombian interaction, to modify

#define NbrBitZ2 5
#define NbrBitY2 5
#define NbrBitX2 5
#define NbrBitZ1 5
#define NbrBitY1 5
#define NbrBitX1 5

#define NbrBit2 (NbrBitX2 + NbrBitY2 + NbrBitZ2)
#define NbrBit1 (NbrBitX1 + NbrBitY1 + NbrBitZ1)

#define Hex2 0x7fff 
#define Hex1 0x3fff8000



#define ShiftZ2 0
#define ShiftY2 (ShiftZ2 + NbrBitZ2)
#define ShiftX2 (ShiftY2 + NbrBitY2)
#define ShiftZ1 (ShiftX2 + NbrBitX2)
#define ShiftY1 (ShiftZ1 + NbrBitZ1)
#define ShiftX1 (ShiftY1 + NbrBitY1)

/*
#define Z2Hex 0x1f
#define Y2Hex 0xfffffe0
#define X2Hex 0x
#define Z1Hex 0x
#define Y1Hex 0x
#define X1Hex 0x
*/

// constructor

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
  this->EvaluateKineticTerm (Mex, Mey, Mez, Mhx, Mhy, Mhz, firstParticle, secondParticle, xSize, ySize, zSize);
  this->EvaluateConfinementTerm (potentialElectron, potentialHole);
  this->EvaluateCoulombianTerm (xSize, ySize, zSize, dielectric);
  delete firstParticle; delete secondParticle;
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

PeriodicElectronHole3DHamiltonian::PeriodicElectronHole3DHamiltonian(const PeriodicElectronHole3DHamiltonian& hamiltonian)
{
  
}

// destructor
//

PeriodicElectronHole3DHamiltonian::~ PeriodicElectronHole3DHamiltonian()
{

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
	      vDestination.Re(Index2) += tmpCoulombian[X2] * vSource.Re(Index1);
	      vDestination.Im(Index2) += tmpCoulombian[X2] * vSource.Im(Index1);
	    }
	  // confinement term for the electrons
	  if (Y1 == Y2)
	    {
	      vDestination.Re(Index2) += (tmpRealElectron[X2] * vSource.Re(Index1) - tmpImaginaryElectron[X2] * vSource.Im(Index1));
	      vDestination.Im(Index2) += (tmpRealElectron[X2] * vSource.Im(Index1) + tmpImaginaryElectron[X2] * vSource.Re(Index1));
	    }
	  // confinement term for the holes
	  if (X1 == X2)
	    {
	      vDestination.Re(Index2) += (tmpRealHole[Y2] * vSource.Re(Index1) - tmpImaginaryHole[Y2] * vSource.Im(Index1));
	      vDestination.Im(Index2) += (tmpRealHole[Y2] * vSource.Im(Index1) + tmpImaginaryHole[Y2] * vSource.Re(Index1));	      
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
  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    for (int n1 = 0; n1 < this->NbrState1Y; ++n1)      
      for (int p1 = 0; p1 < this->NbrState1Z; ++p1)	    	    
	for (int m2 = 0; m2 < this->NbrState2X; ++m2)
	  for (int n2 = 0; n2 < this->NbrState2Y; ++n2)      
	    for (int p2 = 0; p2 < this->NbrState2Z; ++p2)	    
	      {
		this->IToX[index] = (m1 << ShiftX1) | (n1 << ShiftY1) | (p1 << ShiftZ1) | (m2 << ShiftX2) | (n2 << ShiftY2) | (p2 << ShiftZ2);
		++index;
	      }
  this->IToX1 = new int** [this->NbrState1X];
  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    {
      this->IToX1[m1] = new int* [this->NbrState1Y];
      for (int n1 = 0; n1 < this->NbrState1Y; ++n1)   
	{
	  this->IToX1[m1][n1] = new int [this->NbrState1Z];
	  for (int p1 = 0; p1 < this->NbrState1Z; ++p1)	    
	    this->IToX1[m1][n1][p1] = (m1 << ShiftX1) | (n1 << ShiftY1) | (p1 << ShiftZ1);
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
	    this->IToX2[m2][n2][p2] = (m2 << ShiftX2) | (n2 << ShiftY2) | (p2 << ShiftZ2);
	}
    }

}

// evaluate the kinetic term
//

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
  
  double InvXFactor1 = PERIODIC_HAMILTONIAN_FACTOR / (mex * xSize * xSize);
  double InvYFactor1 = PERIODIC_HAMILTONIAN_FACTOR / (mey * ySize * ySize);
  double InvZFactor1 = PERIODIC_HAMILTONIAN_FACTOR / (mez * zSize * zSize);

  double InvXFactor2 = PERIODIC_HAMILTONIAN_FACTOR / (mhx * xSize * xSize);
  double InvYFactor2 = PERIODIC_HAMILTONIAN_FACTOR / (mhy * ySize * ySize);
  double InvZFactor2 = PERIODIC_HAMILTONIAN_FACTOR / (mhz * zSize * zSize);
  
  int index = 0; double FactorX1, FactorY1, FactorZ1, FactorX2, FactorY2;
  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    {
      FactorX1 = InvXFactor1 * ((double) (m1 + LowerImpulsion1X)) * ((double) (m1 + LowerImpulsion1X));
      for (int n1 = 0; n1 < this->NbrState1Y; ++n1)	
	{
	  FactorY1 = FactorX1 + InvYFactor1 * ((double) (n1 + LowerImpulsion1Y)) * ((double) (n1 + LowerImpulsion1Y));
	  for (int p1 = 0; p1 < this->NbrState1Z; ++p1)
	    {
	      FactorZ1 = FactorY1 + InvZFactor1 * ((double) (p1 + LowerImpulsion1Z)) * ((double) (p1 + LowerImpulsion1Z));	      
	      for (int m2 = 0; m2 < this->NbrState2X; ++m2)
		{
		  FactorX2 = FactorZ1 + InvXFactor2 * ((double) (m2 + LowerImpulsion2X)) * ((double) (m2 + LowerImpulsion2X));		
		  for (int n2 = 0; n2 < this->NbrState2Y; ++n2)
		    {
		      FactorY2 = FactorX2 + InvYFactor2 * ((double) (n2 + LowerImpulsion2Y)) * ((double) (n2 + LowerImpulsion2Y));		  
		      for (int p2 = 0; p2 < this->NbrState2Z; ++p2)
			{
			  this->KineticTerm[index] = FactorY2 + InvZFactor2 * ((double) (p2 + LowerImpulsion2Z)) * ((double) (p2 + LowerImpulsion2Z));
			  ++index;
			}
		    }
		}
	    }
	}
    }
}

// evaluate the confinement terms for electrons and holes
//

void PeriodicElectronHole3DHamiltonian::EvaluateConfinementTerm (ThreeDConstantCellPotential* potentialElectron, ThreeDConstantCellPotential* potentialHole)
{
  double*** RealOverlapEX; double*** RealOverlapEY; double*** RealOverlapEZ;
  double*** ImaginaryOverlapEX; double*** ImaginaryOverlapEY; double*** ImaginaryOverlapEZ;
  double*** RealOverlapHX; double*** RealOverlapHY; double*** RealOverlapHZ;
  double*** ImaginaryOverlapHX; double*** ImaginaryOverlapHY; double*** ImaginaryOverlapHZ;
  this->EvaluateWaveFunctionOverlap (potentialElectron->GetNbrCellX (), this->NbrState1X, RealOverlapEX, ImaginaryOverlapEX);
  this->EvaluateWaveFunctionOverlap (potentialElectron->GetNbrCellY (), this->NbrState1Y, RealOverlapEY, ImaginaryOverlapEY);
  this->EvaluateWaveFunctionOverlap (potentialElectron->GetNbrCellZ (), this->NbrState1Z, RealOverlapEZ, ImaginaryOverlapEZ);
  this->EvaluateWaveFunctionOverlap (potentialHole->GetNbrCellX (), this->NbrState2X, RealOverlapHX, ImaginaryOverlapHX);
  this->EvaluateWaveFunctionOverlap (potentialHole->GetNbrCellY (), this->NbrState2Y, RealOverlapHY, ImaginaryOverlapHY);
  this->EvaluateWaveFunctionOverlap (potentialHole->GetNbrCellZ (), this->NbrState2Z, RealOverlapHZ, ImaginaryOverlapHZ);

  int number1 = pow (2, NbrBit1);
  this->RealElectronConfinement = new double* [number1];
  this->ImaginaryElectronConfinement = new double* [number1];
  for (int i1 = 0; i1 < number1; ++i1)
    {
      this->RealElectronConfinement[i1] = new double [number1];
      this->ImaginaryElectronConfinement[i1] = new double [number1];
      for (int i2 = 0; i2 < number1; ++i2)
	{
	  this->RealElectronConfinement[i1][i2] = 0.0;
	  this->ImaginaryElectronConfinement[i1][i2] = 0.0;
	}
    }
  double* tmpRealX; double* tmpRealY; double* tmpRealZ;
  double* tmpImaginaryX; double* tmpImaginaryY; double* tmpImaginaryZ;
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
		  tmpRealX = RealOverlapEX[m1][m3]; 
		  tmpRealY = RealOverlapEY[m1][m3];
		  tmpRealZ = RealOverlapEZ[m1][m3]; 
		  tmpImaginaryX = ImaginaryOverlapEX[m1][m3];
		  tmpImaginaryY = ImaginaryOverlapEY[m1][m3];
		  tmpImaginaryZ = ImaginaryOverlapEZ[m1][m3];
		  double tmpRe1 = 0.0, tmpIm1 = 0.0; 
		  for (int k = 0; k < potentialElectron->GetNbrCellZ (); ++k)
		    {
		      tmpRe1 = 0; tmpIm1 = 0;
		      for (int j = 0; j < potentialElectron->GetNbrCellY (); ++j)
			{
			  for (int i = 0; i < potentialElectron->GetNbrCellX (); ++i)
			    {
			      tmpRe1 += potentialElectron->GetPotential (i, j, k) * (tmpRealX[i] * tmpRealY[j] - tmpImaginaryX[i] * tmpImaginaryY[j]);
			      tmpIm1 += potentialElectron->GetPotential (i, j, k) * (tmpRealX[i] * tmpImaginaryY[j] + tmpRealY[j] * tmpImaginaryX[i]);
			    }
			}
		      this->RealElectronConfinement[X1][X2] += (tmpRe1 * tmpRealZ[k] - tmpIm1 * tmpImaginaryZ[k]);		      
		      this->ImaginaryElectronConfinement[X1][X2] += (tmpRe1 * tmpImaginaryZ[k] + tmpIm1 * tmpRealZ[k]);
		      
		    }
		}
	}

  // hole confinement
  int number2 = pow (2, NbrBit2);
  this->RealHoleConfinement = new double* [number2];
  this->ImaginaryHoleConfinement = new double* [number2];
  for (int i1 = 0; i1 < number2; ++i1)
    {
      this->RealHoleConfinement[i1] = new double [number2];
      this->ImaginaryHoleConfinement[i1] = new double [number2];
      for (int i2 = 0; i2 < number2; ++i2)
	{
	  this->RealHoleConfinement[i1][i2] = 0.0;
	  this->ImaginaryHoleConfinement[i1][i2] = 0.0;
	}
    }
  for (int m2 = 0; m2 < this->NbrState2X; ++m2)
    for (int n2 = 0; n2 < this->NbrState2Y; ++n2)
      for (int p2 = 0; p2 < this->NbrState2Z; ++p2)
	{
	  int Y1 = this->IToX2[m2][n2][p2];	  
	  for (int m4 = 0; m4 < this->NbrState2X; ++m4)
	    for (int n4 = 0; n4 < this->NbrState2Y; ++n4)
	      for (int p4 = 0; p4 < this->NbrState2Z; ++p4)	
		{
		  int Y2 = this->IToX2[m4][n4][p4];
		  tmpRealX = RealOverlapHX[m2][m4]; 
		  tmpRealY = RealOverlapHY[m2][m4];
		  tmpRealZ = RealOverlapHZ[m2][m4]; 
		  tmpImaginaryX = ImaginaryOverlapHX[m2][m4];
		  tmpImaginaryY = ImaginaryOverlapHY[m2][m4];
		  tmpImaginaryZ = ImaginaryOverlapHZ[m2][m4];
		  double tmpRe1 = 0.0, tmpIm1 = 0.0;
		  for (int k = 0; k < potentialHole->GetNbrCellZ (); ++k)
		    {
		      tmpRe1 = 0; tmpIm1 = 0;
		      for (int j = 0; j < potentialHole->GetNbrCellY (); ++j)
			{
			  for (int i = 0; i < potentialHole->GetNbrCellX (); ++i)
			    {
			      tmpRe1 += potentialHole->GetPotential (i, j, k) * (tmpRealX[i] * tmpRealY[j] - tmpImaginaryX[i] * tmpImaginaryY[j]);
			      tmpIm1 += potentialHole->GetPotential (i, j, k) * (tmpRealX[i] * tmpImaginaryY[j] + tmpRealY[j] * tmpImaginaryX[i]);
			    }
			}
		      this->RealHoleConfinement[Y1][Y2] += (tmpRe1 * tmpRealZ[k] - tmpIm1 * tmpImaginaryZ[k]);		      
		      this->ImaginaryHoleConfinement[Y1][Y2] += (tmpRe1 * tmpImaginaryZ[k] + tmpIm1 * tmpRealZ[k]);

		    }
		}
	}

  delete[] RealOverlapEX; delete[] RealOverlapEY; delete[] RealOverlapEZ;
  delete[] RealOverlapHX; delete[] RealOverlapHY; delete[] RealOverlapHZ;
  delete[] ImaginaryOverlapEX; delete[] ImaginaryOverlapEY; delete[] ImaginaryOverlapEZ;
  delete[] ImaginaryOverlapHX; delete[] ImaginaryOverlapHY; delete[] ImaginaryOverlapHZ;
}

// evaluate the Coulombian term
//

void PeriodicElectronHole3DHamiltonian::EvaluateCoulombianTerm (double xSize, double ySize, double zSize, double dielectric)
{
  double factor = COULOMBIAN_FACTOR / dielectric;
  double squareX = xSize * xSize;
  double squareY = ySize * ySize;
  double squareZ = zSize * zSize;

  int number = pow (2, NbrBit1);
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
		    this->CoulombianTerm[X1][X2] = factor / tmp;
		}
	}
}

// evaluate the wave function overlap
//
// nbrStep = number of steps in the given direction
// nbrState = number of states chosen for this direction
// realArray = 3D array containing the real elements of the overlap
// imaginaryArray = 3D array containing the imaginary elements of the overlap

bool PeriodicElectronHole3DHamiltonian::EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double*** &realArray, double*** &imaginaryArray)
{
  double Diff = 0.0;
  double Tmp = 0.0;
  double Tmp1 = 1.0 / ((double) nbrStep);
  realArray = new double** [nbrState];
  imaginaryArray = new double** [nbrState];  
  for (int m1 = 0; m1 < nbrState; ++m1)
    {
      realArray[m1] = new double* [nbrState];
      imaginaryArray[m1] = new double* [nbrState];
      for (int m2 = 0; m2 < nbrState; ++m2)
	{
	  realArray[m1][m2] = new double [nbrStep];
	  imaginaryArray[m1][m2] = new double [nbrStep];
	  int delta = -m1 + m2;
	  if (delta != 0)
	    {
	      Diff = 2.0 * M_PI * ((double) delta);
	      Tmp = Diff / nbrStep;	
	      Diff = 1.0 / Diff;	
	      for (int i = 0; i < nbrStep; ++i)
		{
		  realArray[m1][m2][i] = Diff * (sin(Tmp * (i + 1)) - sin(Tmp * i));
		  imaginaryArray[m1][m2][i] = Diff * (cos(Tmp * (i + 1)) - cos(Tmp * i));
		}
	    }
	  else
	    for (int i = 0; i < nbrStep; ++i)
	      {
		realArray[m1][m2][i] = Tmp1;
		imaginaryArray[m1][m2][i] = 0.0;
	      }	
	}
    }
  return true;
}
