////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2003 Duc-Phuong Nguyen                    //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 24/11/2003                        //
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
#include "Hamiltonian/QuantumDotHamiltonian/NewPeriodicQuantumDots3DHamiltonian.h"
#include "Complex.h"
#include "Vector/ComplexVector.h"
#include "Tools/QuantumDot/Potential/ThreeDPotential.h"

#include <iostream>
#include <math.h>

using std::ostream;

//using std::ostream;
using std::cout;
using std::endl;

#define PERIODIC_HAMILTONIAN_FACTOR 150.4


// constructor from default datas
//
// space = Hilbert space associated to the system
// xSize = system dimension in the x direction (in Angstrom unit)
// ySize = system dimension in the y direction (in Angstrom unit)
// zSize = system dimension in the z direction (in Angstrom unit)
// preConstantRegionSize = region size in the z direction where potential is constant in every direction (region before gradiant zone)
// postConstantRegionSize = region size in the z direction where potential is constant in every direction (region after gradiant zone)
// postConstantRegionPotential = value of the potential in the region after the gradiant zone
// mux = effective mass in the x direction (in electron mass unit)
// muy = effective mass in the y direction (in electron mass unit)
// muz = effective mass in the z direction (in electron mass unit)
// nbrCellX = number of cells in the x direction
// nbrCellY = number of cells in the y direction
// nbrCellZ = number of cells in the z direction
// overlapingFactors = tridimensionnal array where overlaping factors are stored
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

NewPeriodicQuantumDots3DHamiltonian::NewPeriodicQuantumDots3DHamiltonian(Periodic3DOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDPotential* PotentialInput)
{
  this->Space = space;
  this->XSize = xSize;
  this->YSize = ySize;
  this->ZSize = zSize;
  this->Mux = mux;
  this->Muy = muy;
  this->Muz = muz;
  this->NbrCellX = nbrCellX;
  this->NbrCellY = nbrCellY;
  this->NbrCellZ = nbrCellZ;
  this->NbrStateX = this->Space->GetNbrStateX();
  this->LowerImpulsionX = this->Space->GetLowerImpulsionX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->LowerImpulsionY = this->Space->GetLowerImpulsionY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();
  this->InteractionFactors = PotentialInput->Potential;
  this->EvaluateInteractionFactors();
}


// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

NewPeriodicQuantumDots3DHamiltonian::NewPeriodicQuantumDots3DHamiltonian(const NewPeriodicQuantumDots3DHamiltonian& hamiltonian)
{
  this->Space = hamiltonian.Space;
  this->XSize = hamiltonian.XSize;
  this->YSize = hamiltonian.YSize;
  this->ZSize = hamiltonian.ZSize;
  this->Mux = hamiltonian.Mux;
  this->Muy = hamiltonian.Muy;
  this->Muz = hamiltonian.Muz;
  this->NbrCellX = hamiltonian.NbrCellX;
  this->NbrCellY = hamiltonian.NbrCellY;
  this->NbrCellZ = hamiltonian.NbrCellZ;
  this->NbrStateX = this->Space->GetNbrStateX();
  this->LowerImpulsionX = this->Space->GetLowerImpulsionX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->LowerImpulsionY = this->Space->GetLowerImpulsionY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();
  this->KineticElements = hamiltonian.KineticElements;
  this->NbrPrecalculatedDimension = hamiltonian.NbrPrecalculatedDimension;
  this->InteractionFactors = hamiltonian.InteractionFactors;
  this->RealWaveFunctionOverlapX = hamiltonian.RealWaveFunctionOverlapX;
  this->ImaginaryWaveFunctionOverlapX = hamiltonian.ImaginaryWaveFunctionOverlapX;
  this->RealWaveFunctionOverlapY = hamiltonian.RealWaveFunctionOverlapY;
  this->ImaginaryWaveFunctionOverlapY = hamiltonian.ImaginaryWaveFunctionOverlapY;
  this->RealWaveFunctionOverlapZ = hamiltonian.RealWaveFunctionOverlapZ;
  this->ImaginaryWaveFunctionOverlapZ = hamiltonian.ImaginaryWaveFunctionOverlapZ;
  this->RealPrecalculatedHamiltonian =  hamiltonian.RealPrecalculatedHamiltonian;
  this->ImaginaryPrecalculatedHamiltonian =  hamiltonian.ImaginaryPrecalculatedHamiltonian;
}

// destructor
//

NewPeriodicQuantumDots3DHamiltonian::~ NewPeriodicQuantumDots3DHamiltonian()
{
  delete[] this->KineticElements;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* NewPeriodicQuantumDots3DHamiltonian::Clone ()
{
  return new NewPeriodicQuantumDots3DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void NewPeriodicQuantumDots3DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void NewPeriodicQuantumDots3DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->KineticElements[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex NewPeriodicQuantumDots3DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
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

Complex NewPeriodicQuantumDots3DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian for a given range of idinces 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& NewPeriodicQuantumDots3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
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

// only in MonoProcessor Mode!!!!!
ComplexVector& NewPeriodicQuantumDots3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{ 
  int OriginX = this->NbrStateX - 1; int OriginY = this->NbrStateY - 1; int OriginZ = this->NbrStateZ - 1;

  int m1, m2, n1, n2, p1;
  int IndexX, IndexY, IndexZ;
  int** TotalIndex = new int* [this->NbrStateX]; int TmpIndex = 0;
  for (m1 = 0; m1 < this->NbrStateX; ++m1) 
    {
      TotalIndex[m1] = new int [this->NbrStateY];
      for (n1 = 0; n1 < this->NbrStateY; ++n1)	
	{
	  TotalIndex[m1][n1] = (m1 * this->NbrStateY + n1) * this->NbrStateZ;
	  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
	    {	      
	      vDestination.Re(TmpIndex) += vSource.Re(TmpIndex) * this->KineticElements[TmpIndex];
	      vDestination.Im(TmpIndex) += vSource.Im(TmpIndex) * this->KineticElements[TmpIndex];
	      ++TmpIndex;
	    }
	}
    }

  int Index1, Index2;
  double* TmpRealPrecalculatedHamiltonian;
  double* TmpImaginaryPrecalculatedHamiltonian;
  double TmpRe = 0.0; double TmpIm = 0.0;
  int* TmpTotalIndex1; int* TmpTotalIndex2;
  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {
      IndexX = m1 + OriginX;
      TmpTotalIndex1 = TotalIndex[m1];
      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];	  
	  for (n1 = 0; n1 < this->NbrStateY; ++n1)
	    {
	      IndexY = n1 + OriginY;
	      for (n2 = 0; n2 < this->NbrStateY; ++n2)
		{
		  TmpRealPrecalculatedHamiltonian = this->RealPrecalculatedHamiltonian[IndexX][IndexY];
		  TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPrecalculatedHamiltonian[IndexX][IndexY];
		  Index1 = TmpTotalIndex1[n1];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      IndexZ = p1 + OriginZ;
		      TmpRe = 0.0; TmpIm = 0.0;
		      Index2 = TmpTotalIndex2[n2];
		      for (; IndexZ >= p1; --IndexZ, ++Index2)
			{
			  TmpRe += (vSource.Re(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ] - vSource.Im(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ]);
			  TmpIm += (vSource.Re(Index2) * TmpImaginaryPrecalculatedHamiltonian[IndexZ] + vSource.Im(Index2) * TmpRealPrecalculatedHamiltonian[IndexZ]);  	  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		    }
   		  --IndexY;
		}
	    }
	  --IndexX;
	}
    }
  delete[] TotalIndex;
  return vDestination;
}
  
void NewPeriodicQuantumDots3DHamiltonian::EvaluateInteractionFactors()
{
  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellX, this->NbrStateX, this->RealWaveFunctionOverlapX, this->ImaginaryWaveFunctionOverlapX))
    cout << "Error in evaluation of function overlap in X direction. Stop!" << endl;  
  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellY, this->NbrStateY, this->RealWaveFunctionOverlapY, this->ImaginaryWaveFunctionOverlapY))
    cout << "Error in evaluation of function overlap in Y direction. Stop!" << endl;
  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellZ, this->NbrStateZ, this->RealWaveFunctionOverlapZ, this->ImaginaryWaveFunctionOverlapZ))
    cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;

  double InvXFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Mux * this->XSize * this->XSize);
  double InvYFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muy * this->YSize * this->YSize);
  double InvZFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muz * this->ZSize * this->ZSize);
  
  this->KineticElements = new double[this->Space->GetHilbertSpaceDimension ()];

  double FactorX = 0.0, FactorY = 0.0;
  int TotalIndex = 0;
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      FactorX = double((i - this->LowerImpulsionX) * (i - this->LowerImpulsionX)) * InvXFactor;
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  FactorY = double((j - this->LowerImpulsionY) * (j - this->LowerImpulsionY)) * InvYFactor + FactorX;
	  for (int k = 0; k < this->NbrStateZ; ++k)
	    {	      
	      this->KineticElements[TotalIndex] = FactorY + double((k - this->LowerImpulsionZ) * (k - this->LowerImpulsionZ)) * InvZFactor;	      
	      ++TotalIndex;
	    }
	}
    }

  int LengthX = (this->NbrStateX - 1) * 2 + 1; int LengthY = (this->NbrStateY - 1) * 2 + 1; int LengthZ = (this->NbrStateZ - 1) * 2 + 1;

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
      TmpRealWaveFunctionOverlapX = this->RealWaveFunctionOverlapX[m];
      TmpImaginaryWaveFunctionOverlapX = this->ImaginaryWaveFunctionOverlapX[m];	      	  
      for (int n = 0; n < LengthY; ++n)
	{	  
	  TmpReal[m][n] = new double [this->NbrCellZ];
	  TmpImaginary[m][n] = new double [this->NbrCellZ];
	  TmpRealWaveFunctionOverlapY = this->RealWaveFunctionOverlapY[n];
	  TmpImaginaryWaveFunctionOverlapY = this->ImaginaryWaveFunctionOverlapY[n];	  
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];		  
	  for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
	    {
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		{
		  TmpRe2 = TmpRealWaveFunctionOverlapY[CellY];
		  TmpIm2 = TmpImaginaryWaveFunctionOverlapY[CellY];
		  for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
		    {		      
		      TmpRe += this->InteractionFactors[CellZ][CellY][CellX] * (TmpRealWaveFunctionOverlapX[CellX] * TmpRe2 - TmpImaginaryWaveFunctionOverlapX[CellX] * TmpIm2);
		      TmpIm += this->InteractionFactors[CellZ][CellY][CellX] * (TmpRealWaveFunctionOverlapX[CellX] * TmpIm2 + TmpImaginaryWaveFunctionOverlapX[CellX] * TmpRe2);		      
		    }
		}
	      TmpRealPrecalculatedHamiltonian[CellZ] = TmpRe;  
	      TmpImaginaryPrecalculatedHamiltonian[CellZ] = TmpIm;  
	    }
	}
    }

  this->RealPrecalculatedHamiltonian = new double** [LengthX];
  this->ImaginaryPrecalculatedHamiltonian = new double** [LengthX];
  double* TmpRealWaveFunctionOverlapZ;
  double* TmpImaginaryWaveFunctionOverlapZ;
  for (int m = 0; m < LengthX; ++m)
    {
      this->RealPrecalculatedHamiltonian[m] = new double* [LengthY];      
      this->ImaginaryPrecalculatedHamiltonian[m] = new double* [LengthY]; 
      for (int n = 0; n < LengthY; ++n)
	{
	  this->RealPrecalculatedHamiltonian[m][n] = new double [LengthZ];      
	  this->ImaginaryPrecalculatedHamiltonian[m][n] = new double [LengthZ]; 
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];
	  for (int p = 0; p < LengthZ; ++p)
	    {
	      TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[p];
	      TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[p];
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  TmpRe += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ] - TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ]);
		  TmpIm += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ] + TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ]);
		}
	      this->RealPrecalculatedHamiltonian[m][n][p] = TmpRe;
	      this->ImaginaryPrecalculatedHamiltonian[m][n][p] = TmpIm;
	    }
	}
    }
  delete[] TmpReal; delete[] TmpImaginary;
}

bool NewPeriodicQuantumDots3DHamiltonian::EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray)
{
  double Diff = 0.0;
  double Tmp = 0.0;
  double Tmp1 = 1.0 / double (nbrStep);
  int Length = (nbrState - 1) * 2 + 1;
  realArray = new double* [Length];
  imaginaryArray = new double* [Length];  
  int Origin = nbrState - 1;
  for (int delta = 0; delta < Length; ++delta)
    {
      realArray[delta] = new double [nbrStep];
      imaginaryArray[delta] = new double [nbrStep];
      if (delta != Origin)
	{
	  Diff = 2.0 * M_PI * double (delta - Origin);
	  Tmp = Diff / nbrStep;	
	  Diff = 1.0 / Diff;	
	  for (int i = 0; i < nbrStep; ++i)
	    {
	      realArray[delta][i] = Diff * (sin(Tmp * (i + 1)) - sin(Tmp * i));
	      imaginaryArray[delta][i] = Diff * (cos(Tmp * (i + 1)) - cos(Tmp * i));
	    }
	}
      else
	for (int i = 0; i < nbrStep; ++i)
	  {
	    realArray[delta][i] = Tmp1;
	    imaginaryArray[delta][i] = 0.0;
	  }	
    }
  return true;
}
