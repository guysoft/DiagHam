////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//      class of potential in three directions with constant cylinders        //
//                                                                            //
//                        last modification : 04/22/2004                      //
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
#include "Hamiltonian/QuantumDotHamiltonian/CylindricalHamiltonianInMagneticField.h"
#include "Complex.h"
#include "Vector/ComplexVector.h"
#include "Tools/QuantumDot/Potential/ThreeDConstantCylinderPotential.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <math.h>

using std::ostream;
using std::cout;
using std::endl;

#define PERIODIC_HAMILTONIAN_FACTOR 150.4
#define ENERGY_FACTOR 1.16e-4
#define LENGTH_FACTOR 256.8


// constructor from data
//
// space = Hilbert space
// mur = effective mass in plane
// muz = effective mass in Z direction
// bz = Z magnetic field component
// PotentialInput = pointer to a 3D potential with constant value in a cell

CylindricalHamiltonianInMagneticField::CylindricalHamiltonianInMagneticField(VerticalPeriodicParticleInMagneticField* space, double mur, double muz, double bz, ThreeDConstantCylinderPotential* PotentialInput)
{
  this->Space = space;
  int nbrCylinder = PotentialInput->GetNbrCylinderZ();
  for (int k = 0; k < nbrCylinder; ++k)
    this->ZSize += PotentialInput->GetHeight(k);            
  this->Mur = mur;
  this->Muz = muz;
  this->NbrStateR = this->Space->GetNbrStateR();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();

  cout << "Evaluating confinement potential ..." << endl;
  this->EvaluateInteractionFactors(bz, PotentialInput);
  cout << "End of confinement potential evaluation." << endl;
}


// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

CylindricalHamiltonianInMagneticField::CylindricalHamiltonianInMagneticField(const CylindricalHamiltonianInMagneticField& hamiltonian)
{
  this->Space = hamiltonian.Space;
  this->ZSize = hamiltonian.ZSize;
  this->Mur = hamiltonian.Mur;
  this->Muz = hamiltonian.Muz;
  this->NbrStateR = this->Space->GetNbrStateR();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();
  this->PartialDiagonalElement = hamiltonian.PartialDiagonalElement;
  this->RealHamiltonian =  hamiltonian.RealHamiltonian;
  this->ImaginaryHamiltonian =  hamiltonian.ImaginaryHamiltonian;
}

// destructor
//

CylindricalHamiltonianInMagneticField::~ CylindricalHamiltonianInMagneticField()
{
  delete[] this->PartialDiagonalElement;
  delete[] this->RealHamiltonian;
  delete[] this->ImaginaryHamiltonian;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* CylindricalHamiltonianInMagneticField::Clone ()
{
  return new CylindricalHamiltonianInMagneticField(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void CylindricalHamiltonianInMagneticField::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void CylindricalHamiltonianInMagneticField::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->PartialDiagonalElement[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex CylindricalHamiltonianInMagneticField::MatrixElement (RealVector& V1, RealVector& V2)
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

Complex CylindricalHamiltonianInMagneticField::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& CylindricalHamiltonianInMagneticField::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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

ComplexVector& CylindricalHamiltonianInMagneticField::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
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

ComplexVector& CylindricalHamiltonianInMagneticField::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  int* TotalIndex = new int [this->NbrStateR];
  int Index1 = 0; int Index2 = 0;
  for (int n = 0; n < this->NbrStateR; ++n)
    {
      TotalIndex[n] = n * this->NbrStateZ;
      for (int p = 0; p < this->NbrStateZ; ++p)
	{
	  vDestination.Re(Index1) += vSource.Re(Index1) * this->PartialDiagonalElement[Index1];
	  vDestination.Im(Index1) += vSource.Im(Index1) * this->PartialDiagonalElement[Index1];
	  ++Index1;
	}
    }

  int n1, n2, p1, IndexZ, LimitZ;
  double TmpRe = 0.0, TmpIm = 0.0;
  double* TmpRealHamiltonian; double* TmpImaginaryHamiltonian;
  for (n1 = 0; n1 < this->NbrStateR; ++n1)
    {
      for (n2 = 0; n2 < n1; ++n2)
	{
	  Index1 = TotalIndex[n1];
	  TmpRealHamiltonian = this->RealHamiltonian[n1][n2];
	  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[n1][n2];
	  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
	    {
	      Index2 = TotalIndex[n2];
	      for (IndexZ = p1; IndexZ > 0; --IndexZ)
		{
		  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));	 
		  ++Index2;
		}
	      LimitZ = this->NbrStateZ - p1;
	      for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
		{
		  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));	 
		  ++Index2;		  
		}
	      vDestination.Re(Index1) += TmpRe;
	      vDestination.Im(Index1) += TmpIm;
	      ++Index1;
	    }
	}
      for (n2 = n1; n2 < this->NbrStateR; ++n2)
	{
	  Index1 = TotalIndex[n1];
	  TmpRealHamiltonian = this->RealHamiltonian[n2][n1];
	  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[n2][n1];
	  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
	    {
	      Index2 = TotalIndex[n2];
	      for (IndexZ = p1; IndexZ > 0; --IndexZ)
		{
		  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));	 
		  ++Index2;
		}
	      LimitZ = this->NbrStateZ - p1;
	      for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
		{
		  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));	 
		  ++Index2;		  
		}
	      vDestination.Re(Index1) += TmpRe;
	      vDestination.Im(Index1) += TmpIm;
	      ++Index1;
	    }
	}	
    }
  delete[] TotalIndex;
  return vDestination;    
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& CylindricalHamiltonianInMagneticField::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{ 
  /*
  if ((firstComponent == 0) && (nbrComponent == this->Space->GetHilbertSpaceDimension()))
    return this->LowLevelAddMultiply(vSource, vDestination);
  else
    {
      int lastComponent = firstComponent + nbrComponent;
      int OriginX = this->NbrStateX - 1; int OriginY = this->NbrStateY - 1; int OriginZ = this->NbrStateZ - 1;
      int m1, m2, n1, n2, p1;
      int IndexX, IndexY, IndexZ;
      double* TmpRealPrecalculatedHamiltonian;
      double* TmpImaginaryPrecalculatedHamiltonian;
      double TmpRe = 0.0; double TmpIm = 0.0;
      
      int Index1 = firstComponent; int Index2 = 0;
      int ReducedIndex1 = Index1 / this->NbrStateZ;
      p1 = Index1 - ReducedIndex1 * this->NbrStateZ;
      m1 = ReducedIndex1 / this->NbrStateY;
      n1 = ReducedIndex1 - m1 * this->NbrStateY;
      
      // need to add the paramagnetic term here!!!

      for (; Index1 < lastComponent; ++Index1)
	{
	  TmpRe = 0.0; TmpIm = 0.0;
	  TmpRe += vSource.Re(Index1) * this->KineticElements[Index1];
	  TmpIm += vSource.Im(Index1) * this->KineticElements[Index1];
	  Index2 = 0;
	  for (IndexX = m1 + OriginX; IndexX >= m1; --IndexX)
	    for (IndexY = n1 + OriginY; IndexY >= n1; --IndexY)
	      {
		TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		for (IndexZ = p1 + OriginZ; IndexZ >= p1; --IndexZ)
		  {
		    TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
		    TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]);
		    ++Index2;
		    
		  }
	      }	  
	  vDestination.Re(Index1) += TmpRe;
	  vDestination.Im(Index1) += TmpIm;
	  ++p1;
	  if (p1 == this->NbrStateZ)
	    {
	      p1 = 0;
	      ++n1;
	      if (n1 == this->NbrStateY)
		{
		  n1 = 0;
		  ++m1;
		}
	    }
	}
      
      return vDestination;
    }
  */return vDestination;
}

// evaluate all interaction factors
//
// Bz = magnetic field component in Z direction
// potential = pointer to the potential

void CylindricalHamiltonianInMagneticField::EvaluateInteractionFactors(double Bz, ThreeDConstantCylinderPotential* &potential)
{
  double PlaneQuantum = ENERGY_FACTOR * Bz / this->Mur;
  double OrbitRadius = LENGTH_FACTOR / sqrt(Bz);

  double** RealWaveFunctionOverlapZ; double** ImaginaryWaveFunctionOverlapZ;
  if (!this->EvaluatePlaneWaveFunctionOverlap(potential, this->NbrStateZ, RealWaveFunctionOverlapZ, ImaginaryWaveFunctionOverlapZ))
    cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;
  
  int nbrCylinder = potential->GetNbrCylinderZ();
  double** Integral = new double* [nbrCylinder];
  double tmp = 0.0; double radius = 0.0;
  double tk = 0.0;
  for (int k = 0; k < nbrCylinder; ++k)
    {
      Integral[k] = new double [this->NbrStateR * 2 - 1];
      radius = potential->GetRadius(k);
      if (radius > 0.0)
	{
	  tk = (radius * radius) / (2.0 * OrbitRadius * OrbitRadius);
	  tmp = exp(-tk);
	  Integral[k][0] = 1.0 - tmp;
	  for (int n = 1; n < (this->NbrStateR * 2 - 1); ++n)
	    {
	      tmp *= tk;
	      Integral[k][n] = double(n) * Integral[k][n - 1] - tmp;
	    }
	}
    }

  double** LaguerreCoefficient = new double* [this->NbrStateR];
  FactorialCoefficient TmpLaguerre = 1;
  for (int n = 0; n < this->NbrStateR; ++n)
    {
      LaguerreCoefficient[n] = new double [n + 1];      
      LaguerreCoefficient[n][0] = 1.0;
      TmpLaguerre = 1;
      for (int m = 1; m < (n + 1); ++m)
	{
	  TmpLaguerre *= (n + 1 -m);
	  TmpLaguerre /= m;
	  TmpLaguerre /= m;
	  if ((m % 2) == 0)
	    LaguerreCoefficient[n][m] = TmpLaguerre.GetNumericalValue();
	  else
	    LaguerreCoefficient[n][m] = -TmpLaguerre.GetNumericalValue();	    
	}
    }
  
  this->RealHamiltonian = new double** [this->NbrStateR];
  this->ImaginaryHamiltonian = new double** [this->NbrStateR];
  double TmpRe = 0.0; double TmpIm = 0.0; double TmpSum = 0.0;

  // layers having non-constant potential
  for (int n1 = 0; n1 < this->NbrStateR; ++n1)
    {
      this->RealHamiltonian[n1] = new double* [n1 + 1];
      this->ImaginaryHamiltonian[n1] = new double* [n1 + 1];
      for (int n2 = 0; n2 < (n1 + 1); ++n2)
	{
	  this->RealHamiltonian[n1][n2] = new double [this->NbrStateZ];
	  this->ImaginaryHamiltonian[n1][n2] = new double [this->NbrStateZ];
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int k = 0; k < nbrCylinder; ++k)
		{
		  radius = potential->GetRadius(k);
		  if (radius > 0.0)
		    {
		      TmpSum = 0.0;
		      for (int m1 = 0; m1 <= n1; ++m1)
			for (int m2 = 0; m2 <= n2; ++m2) 			  
			  TmpSum += (LaguerreCoefficient[n1][m1] * LaguerreCoefficient[n2][m2] * Integral[k][m1 + m2]);
		      TmpRe += (potential->GetPotential(k) * RealWaveFunctionOverlapZ[p][k] * TmpSum);
		      TmpIm += (potential->GetPotential(k) * ImaginaryWaveFunctionOverlapZ[p][k] * TmpSum);
		    }
		}
	      this->RealHamiltonian[n1][n2][p] = TmpRe;
	      this->ImaginaryHamiltonian[n1][n2][p] = TmpIm;
	    }
	}
    }

  // layers having constant potential
  for (int n = 0; n < this->NbrStateR; ++n)
    for (int p = 0; p < this->NbrStateZ; ++p)
      {	
	for (int k = 0; k < nbrCylinder; ++k)
	  {
	    radius = potential->GetRadius(k);
	    if (radius < 0.0)
	      {
		this->RealHamiltonian[n][n][p] += (potential->GetPotential(k) * RealWaveFunctionOverlapZ[p][k]);
		this->ImaginaryHamiltonian[n][n][p] += (potential->GetPotential(k) * ImaginaryWaveFunctionOverlapZ[p][k]);
	      }
	  }
      }
  // partial diagonal terms
  this->PartialDiagonalElement = new double [this->Space->GetHilbertSpaceDimension()];
  double TmpE = 0.0; int Index = 0;
  double InvZFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muz * this->ZSize * this->ZSize);
  for (int n = 0; n < this->NbrStateR; ++n)
    {
      TmpE = PlaneQuantum * (0.5 + double(n));
      for (int p = 0; p < this->NbrStateZ; ++p)
	{
	  this->PartialDiagonalElement[Index] = TmpE *  double((p + this->LowerImpulsionZ) * (p + this->LowerImpulsionZ)) * InvZFactor; 
  	  ++Index;
	}
    }
}

// evaluate the plane wave function overlap
//
// potential = pointer to the potential
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool CylindricalHamiltonianInMagneticField::EvaluatePlaneWaveFunctionOverlap(ThreeDConstantCylinderPotential* &potential, int nbrState, double** &realArray, double** &imaginaryArray)
{
  int nbrCylinder = potential->GetNbrCylinderZ();
  double* ZPosition = new double [nbrCylinder + 1];
  ZPosition[0] = 0.0;
  for (int k = 0; k < nbrCylinder; ++k)
    ZPosition[k + 1] = ZPosition[k] + potential->GetHeight(k);      
      
  realArray = new double* [nbrState];
  imaginaryArray = new double* [nbrState];   

  realArray[0] = new double [nbrCylinder];
  imaginaryArray[0] = new double [nbrCylinder];
  for (int k = 0; k < nbrCylinder; ++k)
    {
      realArray[0][k] = (ZPosition[k + 1] - ZPosition[k]) / this->ZSize;
      imaginaryArray[0][k] = 0.0;     
    }

  double Diff = 0.0, Tmp = 0.0;
  for (int delta = 1; delta < nbrState; ++delta)
    {
      realArray[delta] = new double [nbrCylinder];
      imaginaryArray[delta] = new double [nbrCylinder];
      Diff = 2.0 * M_PI * double(delta);
      Tmp = Diff / this->ZSize;
      Diff = 1.0 / Diff;
      for (int k = 0; k < nbrCylinder; ++k)
	{
	  realArray[delta][k] = Diff * (sin(Tmp * ZPosition[k + 1]) - sin(Tmp * ZPosition[k]));
	  imaginaryArray[delta][k] = Diff * (cos(Tmp * ZPosition[k]) - cos(Tmp * ZPosition[k + 1]));     
	}
    }
  return true;
}

