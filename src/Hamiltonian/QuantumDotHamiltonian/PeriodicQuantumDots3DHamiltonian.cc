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
#include "Hamiltonian/QuantumDotHamiltonian/PeriodicQuantumDots3DHamiltonian.h"
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

PeriodicQuantumDots3DHamiltonian::PeriodicQuantumDots3DHamiltonian(Confined3DOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDPotential* PotentialInput, int memory)
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
  this->NbrStateY = this->Space->GetNbrStateY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->LeftNumber = PotentialInput->under;
  this->PreConstantRegionSize = PotentialInput->UnderSize;
  this->PreConstantRegionPotential = PotentialInput->UnderPotential;
  this->RightNumber = PotentialInput->above;
  this->PostConstantRegionSize =  PotentialInput->AboveSize;
  this->PostConstantRegionPotential =  PotentialInput->AbovePotential;
  this->InteractionFactors = PotentialInput->Potential;
  this->EvaluateInteractionFactors(memory);
}


// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

PeriodicQuantumDots3DHamiltonian::PeriodicQuantumDots3DHamiltonian(const PeriodicQuantumDots3DHamiltonian& hamiltonian)
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
  this->NbrStateY = this->Space->GetNbrStateY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->DiagonalElements = hamiltonian.DiagonalElements;
  this->NbrPrecalculatedDimension = hamiltonian.NbrPrecalculatedDimension;
  this->InteractionFactors = hamiltonian.InteractionFactors;
  this->RealWaveFunctionOverlapX = hamiltonian.RealWaveFunctionOverlapX;
  this->ImaginaryWaveFunctionOverlapX = hamiltonian.ImaginaryWaveFunctionOverlapX;
  this->ImaginaryWaveFunctionOverlapY = hamiltonian.ImaginaryWaveFunctionOverlapY;
  this->ImaginaryWaveFunctionOverlapZ = hamiltonian.ImaginaryWaveFunctionOverlapZ;
  this->LeftNumber = hamiltonian.LeftNumber;
  this->PreConstantRegionSize = hamiltonian.PreConstantRegionSize;
  this->PreConstantRegionPotential = hamiltonian.PreConstantRegionPotential;
  this->RightNumber = hamiltonian.RightNumber;
  this->PostConstantRegionSize = hamiltonian.PostConstantRegionSize;
  this->PostConstantRegionPotential = hamiltonian.PostConstantRegionPotential;
}

// destructor
//

PeriodicQuantumDots3DHamiltonian::~ PeriodicQuantumDots3DHamiltonian()
{
  delete[] this->DiagonalElements;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* PeriodicQuantumDots3DHamiltonian::Clone ()
{
  return new PeriodicQuantumDots3DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PeriodicQuantumDots3DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void PeriodicQuantumDots3DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->DiagonalElements[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicQuantumDots3DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
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

Complex PeriodicQuantumDots3DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
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

ComplexVector& PeriodicQuantumDots3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
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


ComplexVector& PeriodicQuantumDots3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Space->GetHilbertSpaceDimension();

  double TmpVSourceRe = 0.0, TmpVSourceIm = 0.0;
  double TmpRe = 0.0, TmpIm = 0.0;
  double TmpTotalRe = 0.0, TmpTotalIm = 0.0;
  int Index1 = firstComponent;
  int ReducedIndex1 = Index1 / this->NbrStateZ;
  int k1 = Index1 - ReducedIndex1 * this->NbrStateZ;
  int ReducedFirst = firstComponent / this->NbrStateZ;
  int kFirst = firstComponent - ReducedFirst * this->NbrStateZ;
  int ReducedLast = LastComponent / this->NbrStateZ;
  int kLast = LastComponent - ReducedLast * this->NbrStateZ;
  int k2 = 0;
  int CellZ;
  int ReducedDim = this->NbrStateX * this->NbrStateY;
  double* TmpRealWaveFunctionOverlapZ;
  double* TmpImaginaryWaveFunctionOverlapZ;
  double* TmpRealPrecalculatedHamiltonian;
  double* TmpImaginaryPrecalculatedHamiltonian;
  int ReducedIndex2 = 0;
  int Index2 = 0;  

  while (Index1 < LastComponent)
    {
      ReducedIndex2 = 0;
      Index2 = 0;     
      TmpVSourceRe = vSource.Re(Index1); 
      TmpVSourceIm = vSource.Im(Index1); 
      TmpTotalRe = 0.0; TmpTotalIm = 0.0;
      for (; ReducedIndex2 < ReducedFirst; ++ReducedIndex2)
	{
	  TmpRealPrecalculatedHamiltonian = this->RealPartial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	  TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPartial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	  for (k2 = 0; k2 < this->NbrStateZ; ++k2)
	    {
	      TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[k1][k2];
	      TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[k1][k2];
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  TmpRe += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ] - TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpImaginaryPrecalculatedHamiltonian[CellZ]);
		  TmpIm += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpImaginaryPrecalculatedHamiltonian[CellZ] + TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);		  
		}
	      TmpTotalRe += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
	      TmpTotalIm += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));
	      ++Index2;
	    }
	}
      TmpRealPrecalculatedHamiltonian = this->RealPartial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
      TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPartial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
      if (ReducedIndex2 == ReducedIndex1)
	for (k2 = 0; k2 < kFirst; ++k2)
	  {
	    TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[k1][k2];
	    TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[k1][k2];
	    TmpRe = 0.0; TmpIm = 0.0;
	    for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
	      {
		TmpRe += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		TmpIm += (TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
	      }
	    TmpRe += RealPartialZPrecalculatedHamiltonian[k1][k2];
	    TmpIm += ImaginaryPartialZPrecalculatedHamiltonian[k1][k2];
	    TmpTotalRe += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
	    TmpTotalIm += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));	    
	    ++Index2;
	  }
      else
	for (k2 = 0; k2 < kFirst; ++k2)
	  {
	    TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[k1][k2];
	    TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[k1][k2];
	    TmpRe = 0.0; TmpIm = 0.0;
	    for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
	      {
		TmpRe += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ] - TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpImaginaryPrecalculatedHamiltonian[CellZ]);
		TmpIm += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpImaginaryPrecalculatedHamiltonian[CellZ] + TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);		  
	      }
	    TmpTotalRe += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
	    TmpTotalIm += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));
	      ++Index2;
	  }
      // intersection
      if (ReducedIndex2 < ReducedIndex1)
	{
	  for (; k2 < this->NbrStateZ; ++k2)
	    {	   
	      TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[k1][k2];
	      TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[k1][k2];
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  TmpRe += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ] - TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpImaginaryPrecalculatedHamiltonian[CellZ]);
		  TmpIm += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpImaginaryPrecalculatedHamiltonian[CellZ] + TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);		  
		}
	      TmpTotalRe += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
	      TmpTotalIm += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));
	      vDestination.Re(Index2) += TmpRe * TmpVSourceRe + TmpIm * TmpVSourceIm;
	      vDestination.Im(Index2) += TmpRe * TmpVSourceIm - TmpIm * TmpVSourceRe;
	      ++Index2;
	    }
	  ++ReducedIndex2;
	  for (; ReducedIndex2 < ReducedIndex1; ++ReducedIndex2)
	    {
	      TmpRealPrecalculatedHamiltonian = this->RealPartial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	      TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPartial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	      for (k2 = 0; k2 < this->NbrStateZ; ++k2)
		{
		  TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[k1][k2];
		  TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[k1][k2];
		  TmpRe = 0.0; TmpIm = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		    {
		      TmpRe += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ] - TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpImaginaryPrecalculatedHamiltonian[CellZ]);
		      TmpIm += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpImaginaryPrecalculatedHamiltonian[CellZ] + TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);		  
		    }
		  TmpTotalRe += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
		  TmpTotalIm += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));
		  vDestination.Re(Index2) += TmpRe * TmpVSourceRe + TmpIm * TmpVSourceIm;
		  vDestination.Im(Index2) += TmpRe * TmpVSourceIm - TmpIm * TmpVSourceRe;		  
		  ++Index2;
		}
	    }
	  TmpRealPrecalculatedHamiltonian = this->RealPartial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex1];
	  // TmpImaginaryPrecalculatedHamiltonian = 0.0
	  for (k2 = 0; k2 < k1; ++k2)
	    {
	      TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[k1][k2];
	      TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[k1][k2];
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)	
		{	
		  TmpRe += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		  TmpIm += (TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		}
	      TmpRe += RealPartialZPrecalculatedHamiltonian[k1][k2];
	      TmpIm += ImaginaryPartialZPrecalculatedHamiltonian[k1][k2];
	      TmpTotalRe += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
	      TmpTotalIm += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));	     
	      vDestination.Re(Index2) += TmpRe * TmpVSourceRe + TmpIm * TmpVSourceIm;
	      vDestination.Im(Index2) += TmpRe * TmpVSourceIm - TmpIm * TmpVSourceRe;
	      ++Index2;
	    }  	  
	}
      else
	{
	  TmpRealPrecalculatedHamiltonian = this->RealPartial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex1];
	  // TmpImaginaryPrecalculatedHamiltonian = 0.0	  
	  for (k2 = kFirst; k2 < k1; ++k2)
	    {
	      TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[k1][k2];
	      TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[k1][k2];
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)	
		{	
		  TmpRe += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		  TmpIm += (TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		}
	      TmpRe += RealPartialZPrecalculatedHamiltonian[k1][k2];
	      TmpIm += ImaginaryPartialZPrecalculatedHamiltonian[k1][k2];
	      TmpTotalRe += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
	      TmpTotalIm += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));	     
	      vDestination.Re(Index2) += TmpRe * TmpVSourceRe + TmpIm * TmpVSourceIm;
	      vDestination.Im(Index2) += TmpRe * TmpVSourceIm - TmpIm * TmpVSourceRe;
	      ++Index2;
	    }
	}
      // diagonal
      TmpTotalRe += this->DiagonalElements[Index1] * vSource.Re(Index1);
      TmpTotalIm += this->DiagonalElements[Index1] * vSource.Im(Index1);
      // second part
      Index2 = LastComponent;
      ReducedIndex2 = ReducedLast;	  
      k2 = kLast;
      if (kLast != 0)
	{
	  TmpRealPrecalculatedHamiltonian = this->RealPartial2DPrecalculatedHamiltonian[ReducedIndex2][ReducedIndex1];
	  TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPartial2DPrecalculatedHamiltonian[ReducedIndex2][ReducedIndex1]; // inversed sign!!!
	  if (ReducedIndex1 == ReducedIndex2)
	    for (; k2 < this->NbrStateZ; ++k2)
	      {	
		TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[k2][k1];
		TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[k2][k1];
		TmpRe = 0.0; TmpIm = 0.0;
		for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)	
		  {	
		    TmpRe += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		    TmpIm -= (TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		  }
		TmpRe += RealPartialZPrecalculatedHamiltonian[k1][k2];
		TmpIm += ImaginaryPartialZPrecalculatedHamiltonian[k1][k2];
		TmpTotalRe += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
		TmpTotalIm += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));	          
		++Index2;
	      }
	  else
	    for (; k2 < this->NbrStateZ; ++k2)
	      {
		TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[k2][k1];
		TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[k2][k1];
		TmpRe = 0.0; TmpIm = 0.0;
		for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)	
		  {	
		    TmpRe += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		    TmpIm -= (TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		  }
		TmpTotalRe += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
		TmpTotalIm += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));
		++Index2;
	      }
	  ++ReducedIndex2;
	}
      for (; ReducedIndex2 < ReducedDim; ++ReducedIndex2)
	{
	  TmpRealPrecalculatedHamiltonian = this->RealPartial2DPrecalculatedHamiltonian[ReducedIndex2][ReducedIndex1];
	  TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPartial2DPrecalculatedHamiltonian[ReducedIndex2][ReducedIndex1];
	  for (k2 = 0; k2 < this->NbrStateZ; ++k2)
	    {
	      TmpRealWaveFunctionOverlapZ = this->RealWaveFunctionOverlapZ[k2][k1];
	      TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[k2][k1];
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)	
		{	
		  TmpRe += (TmpRealWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		  TmpIm -= (TmpImaginaryWaveFunctionOverlapZ[CellZ] * TmpRealPrecalculatedHamiltonian[CellZ]);
		}
	      TmpTotalRe += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
	      TmpTotalIm += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));
	      ++Index2;
	    }
	}
      vDestination.Re(Index1) += TmpTotalRe;
      vDestination.Im(Index1) += TmpTotalIm;
      ++k1;
      if (k1 == this->NbrStateZ)
	{
	  k1 = 0;
	  ++ReducedIndex1;
	}
      ++Index1;
    }
  return vDestination;            
}
  
void PeriodicQuantumDots3DHamiltonian::EvaluateInteractionFactors(int memory)
{
  int UsedMemory = 0;
  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellX, this->NbrStateX, UsedMemory, this->RealWaveFunctionOverlapX, this->ImaginaryWaveFunctionOverlapX))
    cout << "Error in evaluation of function overlap in X direction. Stop!" << endl;
  
  if (!this->EvaluateWaveFunctionOverlap(this->NbrCellY, this->NbrStateY, UsedMemory, this->RealWaveFunctionOverlapY, this->ImaginaryWaveFunctionOverlapY))
    cout << "Error in evaluation of function overlap in Y direction. Stop!" << endl;

  if (!this->EvaluateWaveFunctionZOverlap(UsedMemory, this->RealWaveFunctionOverlapZ, this->ImaginaryWaveFunctionOverlapZ))
    cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;

  double InvXFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Mux * this->XSize * this->XSize);
  double InvYFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muy * this->YSize * this->YSize);
  double InvZFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muz * this->ZSize * this->ZSize);
  
  this->DiagonalElements = new double[this->Space->GetHilbertSpaceDimension ()];

  double Factor = 0.0, FactorX = 0.0, FactorY = 0.0;
  double TmpElement = 0.0;
  int TotalIndex = 0;
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      FactorX = double(i * i) * InvXFactor;
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  FactorY = double(j * j) * InvYFactor + FactorX;
	  for (int k = 0; k < this->NbrStateZ; ++k)
	    {	      
	      TmpElement = FactorY + double(k * k) * InvZFactor + this->RealPartialZPrecalculatedHamiltonian[k][k];
	      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		    {
		      Factor = this->RealWaveFunctionOverlapZ[k][k][CellZ] * this->RealWaveFunctionOverlapY[j][j][CellY];		      
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)			
			TmpElement += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->RealWaveFunctionOverlapX[i][i][CellX];	
		    }
		}
	      this->DiagonalElements[TotalIndex] = TmpElement;	      
	      ++TotalIndex;
	    }
	}
    }

  int Dim = this->NbrStateX * this->NbrStateY;
  this->RealPartial2DPrecalculatedHamiltonian = new double** [Dim];
  this->ImaginaryPartial2DPrecalculatedHamiltonian = new double** [Dim]; 

  double TmpRe, TmpIm;
  double TmpRe2, TmpIm2;
  double* TmpRealWaveFunctionOverlapX;
  double* TmpImaginaryWaveFunctionOverlapX = new double[this->NbrCellX];
  double* TmpRealWaveFunctionOverlapY;
  double* TmpImaginaryWaveFunctionOverlapY = new double[this->NbrCellY];
  double* TmpRealPrecalculatedHamiltonian;
  double* TmpImaginaryPrecalculatedHamiltonian;
  int j1 = 0;
  int i1 = 0;
  for (int Index1 = 0; Index1 < Dim; ++Index1)
    {
      this->RealPartial2DPrecalculatedHamiltonian[Index1] = new double* [Index1 + 1];
      this->ImaginaryPartial2DPrecalculatedHamiltonian[Index1] = new double* [Index1 + 1];
      int Index2 = 0;
      int i2 = 0;
      int j2 = 0;
      while (Index2 <= Index1)
	{
	  TmpRealWaveFunctionOverlapX = this->RealWaveFunctionOverlapX[i1][i2];
	  TmpImaginaryWaveFunctionOverlapX = this->ImaginaryWaveFunctionOverlapX[i1][i2];	      	  
	  TmpRealWaveFunctionOverlapY = this->RealWaveFunctionOverlapY[j1][j2];
	  TmpImaginaryWaveFunctionOverlapY = this->ImaginaryWaveFunctionOverlapY[j1][j2];	      	    
	    
	  this->RealPartial2DPrecalculatedHamiltonian[Index1][Index2] = new double [this->NbrCellZ];
	  this->ImaginaryPartial2DPrecalculatedHamiltonian[Index1][Index2] = new double [this->NbrCellZ];
	  TmpRealPrecalculatedHamiltonian = this->RealPartial2DPrecalculatedHamiltonian[Index1][Index2];
	  TmpImaginaryPrecalculatedHamiltonian = this->ImaginaryPartial2DPrecalculatedHamiltonian[Index1][Index2];
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
	  ++Index2;
	  ++j2;
	  if (j2 == this->NbrStateY)
	    {
	      j2 = 0;
	      ++i2;
	    }
	}
      ++j1;
      if (j1 == this->NbrStateY)
	{
	  j1 = 0;
	  ++i1;
	}
    }
  
}

bool PeriodicQuantumDots3DHamiltonian::EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, int& memory, double*** &realArray, double*** &imaginaryArray)
{
  double Diff = 0.0;
  double Tmp = 0.0;
  double Tmp1 = 1.0 / double (nbrStep);

  realArray = new double** [nbrState];
  imaginaryArray = new double** [nbrState];  
  for (int m1 = 0; m1 < nbrState; ++m1)
    {
      realArray[m1] = new double* [nbrState];
      imaginaryArray[m1] = new double* [nbrState];
      for (int m2 = 0; m2 < m1; ++m2)
	{
	  Diff = 2.0 * M_PI * double (m2 - m1);
	  Tmp = Diff / nbrStep;
	  Diff = 1.0 / Diff;
	  realArray[m1][m2] = new double [nbrStep];
	  imaginaryArray[m1][m2] = new double [nbrStep];	
  	  realArray[m2][m1] = new double [nbrStep];
	  imaginaryArray[m2][m1] = new double [nbrStep];
	  for (int i = 0; i < nbrStep; ++i)
	    {
	      realArray[m1][m2][i] = Diff * (sin(Tmp * (i + 1)) - sin(Tmp * i));
	      imaginaryArray[m1][m2][i] = Diff * (cos(Tmp * i) - cos(Tmp * (i + 1)));
	      realArray[m2][m1][i] = realArray[m1][m2][i];
	      imaginaryArray[m2][m1][i] = -imaginaryArray[m1][m2][i];
	    }
	}
      realArray[m1][m1] = new double [nbrStep];
      imaginaryArray[m1][m1] = new double [nbrStep];
      for (int i = 0; i < nbrStep; ++i)
	{
	  realArray[m1][m1][i] = Tmp1;
	  imaginaryArray[m1][m1][i] = 0.0;
	}  
    }
  return true;
}


// evaluate wave function overlaps on a cell in the z direction
//
// memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
// return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)

bool PeriodicQuantumDots3DHamiltonian::EvaluateWaveFunctionZOverlap(int& memory, double*** &realArray, double*** &imaginaryArray)
{
  memory += 2 * ((((this->NbrStateZ * (this->NbrStateZ + 1)) >> 1) * ((this->NbrCellZ)  * sizeof(double) + sizeof(double*))) + this->NbrStateZ * sizeof(double**));

  realArray = new double** [this->NbrStateZ];
  imaginaryArray = new double** [this->NbrStateZ];

  this->RealPartialZPrecalculatedHamiltonian = new double* [this->NbrStateZ];
  this->ImaginaryPartialZPrecalculatedHamiltonian = new double* [this->NbrStateZ];

  double LeftRegionSize = 0.0;
  for (int k = 0; k < this->LeftNumber; ++k)
    LeftRegionSize += PreConstantRegionSize[k];
  double RightRegionSize = 0.0;
  for (int k = 0; k < this->RightNumber; ++k)
    RightRegionSize += PostConstantRegionSize[k];

  double Step = (this->ZSize - LeftRegionSize - RightRegionSize) / ((double) this->NbrCellZ) ;
  double PostShift = this->ZSize - RightRegionSize;
  double Diff;
  double Tmp1 = 0.0; double Tmp2 = 0.0; 
  double Tmp3 = 0.0; double Tmp4 = 0.0; 

  double* PreLength = new double[this->LeftNumber+ 1];
  PreLength[0] = 0.0;
  for (int k = 0; k < this->LeftNumber; ++k)
    PreLength[k + 1] = PreLength[k] + PreConstantRegionSize[k];
  double* PostLength = new double[this->RightNumber + 1];
  PostLength[0] = PostShift;
  for (int k = 0; k < this->RightNumber; ++k)
    PostLength[k + 1] = PostLength[k] + PostConstantRegionSize[k];
  double* Length = new double[this->NbrCellZ + 1];
  Length[0] = LeftRegionSize;
  for (int k = 0; k < this->NbrCellZ; ++k)
    Length[k + 1] = Length[k] + Step;

  for (int i = 0; i < this->NbrStateZ; ++i)
    {
      realArray[i] = new double* [this->NbrStateZ];
      imaginaryArray[i] = new double* [this->NbrStateZ];
      this->RealPartialZPrecalculatedHamiltonian[i] = new double [this->NbrStateZ];
      this->ImaginaryPartialZPrecalculatedHamiltonian[i] = new double [this->NbrStateZ];
      for (int j = 0; j < i; ++j)
	{
	  
	  Diff = 2.0 * M_PI * (double (i - j));
	  Tmp1 = 1.0 / Diff;
	  Tmp2 = Diff / this->ZSize;
	  //TmpLength1 = 0.0; TmpLength2 = 0.0;
	  Tmp3 = 0.0; Tmp4 = 0.0;
	  
	  for (int k = 0; k < this->LeftNumber; ++k)
	    {	      
	      Tmp3 += Tmp1 * (sin(Tmp2 * PreLength[k + 1]) - (sin(Tmp2 * PreLength[k]))) * (this->PreConstantRegionPotential[k]);
	      Tmp4 += Tmp1 * (cos(Tmp2 * PreLength[k + 1]) - (cos(Tmp2 * PreLength[k]))) * (this->PreConstantRegionPotential[k]);	      
	    }

	  for (int k = 0; k < this->RightNumber; ++k)
	    {	  
	      Tmp3 += Tmp1 * (sin(Tmp2 * PostLength[k + 1]) - (sin(Tmp2 * PostLength[k]))) * (this->PostConstantRegionPotential[k]);
	      Tmp4 += Tmp1 * (cos(Tmp2 * PostLength[k + 1]) - (cos(Tmp2 * PostLength[k]))) * (this->PostConstantRegionPotential[k]);
	    }

	  this->RealPartialZPrecalculatedHamiltonian[i][j] = Tmp3;	  
	  this->ImaginaryPartialZPrecalculatedHamiltonian[i][j] = Tmp4;
	  realArray[i][j] = new double [this->NbrCellZ];
	  imaginaryArray[i][j] = new double [this->NbrCellZ];
	  this->RealPartialZPrecalculatedHamiltonian[j][i] = Tmp3;	  
	  this->ImaginaryPartialZPrecalculatedHamiltonian[j][i] = -Tmp4;
	  realArray[j][i] = new double [this->NbrCellZ];
	  imaginaryArray[j][i] = new double [this->NbrCellZ];
	  for (int k = 0; k < this->NbrCellZ; ++k)
	    {
	      realArray[i][j][k] = Tmp1 * (sin(Tmp2 * Length[k + 1]) - (sin(Tmp2 * Length[k]))) ;
	      imaginaryArray[i][j][k] = Tmp1 * (cos(Tmp2 * Length[k + 1]) - (cos(Tmp2 * Length[k])));
	      realArray[j][i][k] = realArray[i][j][k] ;
	      imaginaryArray[j][i][k] = -imaginaryArray[i][j][k];
	    }
	}
      Tmp3 = 0.0; Tmp4 = 0.0;
      for (int k = 0; k < this->LeftNumber; ++k)	
	Tmp3 +=( this->PreConstantRegionPotential[k] * PreConstantRegionSize[k] / this->ZSize);
      for (int k = 0; k < this->RightNumber; ++k)
	Tmp3 += (this->PostConstantRegionPotential[k] * PostConstantRegionSize[k] / this->ZSize); 
      this->RealPartialZPrecalculatedHamiltonian[i][i] = Tmp3;	  
      this->ImaginaryPartialZPrecalculatedHamiltonian[i][i] = 0.0;
      realArray[i][i] = new double [this->NbrCellZ];
      imaginaryArray[i][i] = new double [this->NbrCellZ];
      for (int k = 0; k < this->NbrCellZ; ++k)
	{
	  realArray[i][i][k] = Step / this->ZSize;
	  imaginaryArray[i][i][k] = 0.0;
	}
    }
  return true;
}
