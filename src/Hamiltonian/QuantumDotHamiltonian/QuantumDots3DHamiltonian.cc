////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 26/02/2003                        //
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
#include "Hamiltonian/QuantumDots3DHamiltonian.h"
#include "Vector/RealVector.h"
#include "Complex.h"
#include "Potential/ThreeDPotential.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


#define HAMILTONIAN_FACTOR 37.60



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

QuantumDots3DHamiltonian::QuantumDots3DHamiltonian(Confined3DOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDPotential* PotentialInput, int memory)
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
  this->TotalNbrCells = this->NbrCellX * this->NbrCellY * this->NbrCellZ;
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

QuantumDots3DHamiltonian::QuantumDots3DHamiltonian(const QuantumDots3DHamiltonian& hamiltonian)
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
  this->TotalNbrCells = hamiltonian.TotalNbrCells; 
  this->NbrStateX = this->Space->GetNbrStateX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->DiagonalElements = hamiltonian.DiagonalElements;
  this->NbrPrecalculatedDimension = hamiltonian.NbrPrecalculatedDimension;
  this->InteractionFactors = hamiltonian.InteractionFactors;
  this->WaveFunctionOverlapX = hamiltonian.WaveFunctionOverlapX;
  this->WaveFunctionOverlapY = hamiltonian.WaveFunctionOverlapY;
  this->WaveFunctionOverlapZ = hamiltonian.WaveFunctionOverlapZ;
  this->FullPrecalculatedHamiltonian = hamiltonian.FullPrecalculatedHamiltonian;
  this->LeftNumber = hamiltonian.LeftNumber;
  this->PreConstantRegionSize = hamiltonian.PreConstantRegionSize;
  this->PreConstantRegionPotential = hamiltonian.PreConstantRegionPotential;
  this->RightNumber = hamiltonian.RightNumber;
  this->PostConstantRegionSize = hamiltonian.PostConstantRegionSize;
  this->PostConstantRegionPotential = hamiltonian.PostConstantRegionPotential;
}

// destructor
//

QuantumDots3DHamiltonian::~QuantumDots3DHamiltonian()
{
  delete[] this->DiagonalElements;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* QuantumDots3DHamiltonian::Clone ()
{
  return new QuantumDots3DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void QuantumDots3DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void QuantumDots3DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->DiagonalElements[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex QuantumDots3DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
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

Complex QuantumDots3DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& QuantumDots3DHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination)
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

RealVector& QuantumDots3DHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    vDestination[i] = 0.0;
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& QuantumDots3DHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
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

RealVector& QuantumDots3DHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
							  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Space->GetHilbertSpaceDimension();

  // no acceleration
  if (this->NbrPrecalculatedDimension == 0)
    {
      double Tmp;      
      double Factor;
      int i1 = 0;
      int j1 = 0;
      int k1 = 0;
      for (int Index1 = firstComponent; Index1 < LastComponent; ++Index1)
	{
	  vDestination[Index1] += this->DiagonalElements[Index1] * vSource[Index1];
	  int Index2 = 0;
	  int k2 = 0;
	  for (; k2 < k1; ++k2)
	    {
	      int j2 = 0;
	      for (; j2 <= j1; ++j2)
		{
		  int i2 = 0;
		  for (; i2 <= i1; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k1][k2][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k1][k2][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		}
	      for (; j2 < this->NbrStateY; ++j2)
		{
		  int i2 = 0;
		  for (; i2 <= i1; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor =  this->WaveFunctionOverlapZ[k1][k2][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor =  this->WaveFunctionOverlapZ[k1][k2][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		}
	    }
	  {
	    int j2 = 0;
	    for (; j2 < j1; ++j2)
	      {
		int i2 = 0;
		for (; i2 <= i1; ++i2)
		  {
		    Tmp = 0.0;
		    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			      {
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
			      }	 
			  }     
		      }
		    vDestination[Index1] += Tmp * vSource[Index2];
		    ++Index2;
		  }
		for (; i2 < this->NbrStateX; ++i2)
		  {
		    Tmp = 0.0;
		    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			      {
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
			      }	 
			  }     
		      }
		    vDestination[Index1] += Tmp * vSource[Index2];
		    ++Index2;
		  }
	      }
	    {
	      int i2 = 0;
	      for (; i2 < i1; ++i2)
		{
		  Tmp = 0.0;
		  for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		    {
		      for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			{
			  Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j1][j1][CellY];
			  for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			    {
			      Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
			    }	 
			}     
		    }
		  vDestination[Index1] += Tmp * vSource[Index2];
		  ++Index2;
		}
	      ++i2;
	      ++Index2;
	      for (; i2 < this->NbrStateX; ++i2)
		{
		    Tmp = 0.0;
		    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j1][j1][CellY];
			    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			      {
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
			      }	 
			  }     
		      }
		  vDestination[Index1] += Tmp * vSource[Index2];
		  ++Index2;
		}
	      
	    }
	    ++j2;
	    for (; j2 < this->NbrStateY; ++j2)
	      {
		int i2 = 0;
		for (; i2 <= i1; ++i2)
		  {
		    Tmp = 0.0;
		    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			      {
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
			      }	 
			  }     
		      }
		    vDestination[Index1] += Tmp * vSource[Index2];
		    ++Index2;
		  }
		for (; i2 < this->NbrStateX; ++i2)
		  {
		    Tmp = 0.0;
		    for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		      {
			for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			  {
			    Factor = this->WaveFunctionOverlapZ[k1][k1][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			    for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
			      {
				Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
			      }	 
			  }     
		      }
		    vDestination[Index1] += Tmp * vSource[Index2];
		    ++Index2;
		  }
	      }
	  }
	  ++k2;
	  for (; k2 < this->NbrStateZ; ++k2)
	    {
	      int j2 = 0;
	      for (; j2 <= j1; ++j2)
		{
		  int i2 = 0;
		  for (; i2 <= i1; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k2][k1][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k2][k1][CellZ] * this->WaveFunctionOverlapY[j1][j2][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		}
	      for (; j2 < this->NbrStateY; ++j2)
		{
		  int i2 = 0;
		  for (; i2 <= i1; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k2][k1][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i1][i2][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		  for (; i2 < this->NbrStateX; ++i2)
		    {
		      Tmp = 0.0;
		      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
			{
			  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
			    {
			      Factor = this->WaveFunctionOverlapZ[k2][k1][CellZ] * this->WaveFunctionOverlapY[j2][j1][CellY];
			      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)
				{
				  Tmp += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i2][i1][CellX];
				}	 
			    }     
			}
		      vDestination[Index1] += Tmp * vSource[Index2];
		      ++Index2;
		    }
		}
	      ++i1;
	      if (i1 == this->NbrStateX)
		{
		  i1 = 0;
		  ++j1;
		  if (j1 == this->NbrStateY)
		    {
		      j1 = 0;
		      ++k1;		      
		    }
		}
	    }	      
	}
      return vDestination;
    }
/*
  if (this->NbrPrecalculatedDimension == 0x010)
    {
      double Tmp = 0.0;
      double TmpTotal = 0.0;
      int Index1 = firstComponent;
      int ReducedIndex1 = Index1 / this->NbrStateZ;
      int k1 = Index1 - ReducedIndex1 * this->NbrStateZ;
      int k2 = 0;
      int CellZ;
      int ReducedDim = this->NbrStateX * this->NbrStateY;
      double* TmpWaveFunctionOverlapZ;
      double* TmpPrecalculatedHamiltonian;
      int ReducedIndex2 = 0;
      int Index2 = 0;

      while (Index1 < LastComponent)
	{
	  ReducedIndex2 = 0;
	  Index2 = 0;
	  Tmp = 0.0;
	  TmpTotal = 0.0;
	  for (; ReducedIndex2 < ReducedIndex1; ++ReducedIndex2)
	    {
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	      for (k2 = 0; k2 < this->NbrStateZ; ++k2)
		{	     
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];		    
		  TmpTotal += Tmp * vSource[Index2];
		  ++Index2;
		}
	    }
	  TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	  for (k2 = 0; k2 < k1; ++k2)
	    {
	      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
	      Tmp = 0.0;
	      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		
		Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];		
	      Tmp += this->PartialZPrecalculatedHamiltonian[k1][k2];	      
	      TmpTotal += Tmp * vSource[Index2];
	      ++Index2;
	    }
       
	  TmpTotal += this->DiagonalElements[Index1] * vSource[Index1];	      
	  ++k2;
	  ++Index2;
	  for (; k2 < this->NbrStateZ; ++k2)
	    {
	      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
	      Tmp = 0.0;
	      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		
		Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];		
	      Tmp += this->PartialZPrecalculatedHamiltonian[k1][k2];	     
	      TmpTotal += Tmp * vSource[Index2];
	      ++Index2;
	    }
	  ++ReducedIndex2;
	  for (; ReducedIndex2 < ReducedDim; ++ReducedIndex2)
	    {
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	      for (k2 = 0; k2 < this->NbrStateZ; ++k2)
		{		  
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];		    		
		  TmpTotal += Tmp * vSource[Index2];
		  ++Index2;
		}
	    }
	  vDestination[Index1] += TmpTotal;
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
*/
  if (this->NbrPrecalculatedDimension == 0x001)
    {
      double TmpVSource = 0.0;
      double Tmp = 0.0;
      double TmpTotal = 0.0;
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
      double* TmpWaveFunctionOverlapZ;
      double* TmpPrecalculatedHamiltonian;

      while (Index1 < LastComponent)
	{
	  int ReducedIndex2 = 0;
	  int Index2 = 0;
	  Tmp = 0.0;
	  TmpTotal = 0.0;
	  TmpVSource = vSource[Index1];	  
	  for (; ReducedIndex2 < ReducedFirst; ++ReducedIndex2)
	    {
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	      for (k2 = 0; k2 < this->NbrStateZ; ++k2)
		{	     
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		  TmpTotal += Tmp * vSource[Index2];
		  ++Index2;
		}
	    }
	  TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
	  if (ReducedIndex2 == ReducedIndex1)	    
	    for (k2 = 0; k2 < kFirst; ++k2)
	      {
		TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		Tmp = 0.0;
		for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		
		  Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ]; 	  
		Tmp += PartialZPrecalculatedHamiltonian[k1][k2];
		TmpTotal += Tmp * vSource[Index2];
		++Index2;
	      }
	  else
	    for (k2 = 0; k2 < kFirst; ++k2)
	      {
		TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		Tmp = 0.0;
		for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		
		  Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ]; 	  
		TmpTotal += Tmp * vSource[Index2];
		++Index2;
	      }		      
	  // intersection
	  if (ReducedIndex2 < ReducedIndex1)
	    {
	      for (; k2 < this->NbrStateZ; ++k2)
		{	      
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		  TmpTotal += Tmp * vSource[Index2];
		  vDestination[Index2] += Tmp * TmpVSource;
		  ++Index2;			      
		}
	      ++ReducedIndex2;
	      for (; ReducedIndex2 < ReducedIndex1; ++ReducedIndex2)
		{
		  TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex2];
		  for (k2 = 0; k2 < this->NbrStateZ; ++k2)
		    {	  
		      TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		      Tmp = 0.0;
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		      TmpTotal += Tmp * vSource[Index2];
		      vDestination[Index2] += Tmp * TmpVSource;
		      ++Index2;
		    }
		}
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex1];
	      for (k2 = 0; k2 < k1; ++k2)
		{
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];		
		  Tmp += this->PartialZPrecalculatedHamiltonian[k1][k2];	      
		  TmpTotal += Tmp * vSource[Index2];
		  vDestination[Index2] += Tmp * TmpVSource;
		  ++Index2;
		}  	  
	    }
	  else
	    {
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex1][ReducedIndex1];
	      for (k2 = kFirst; k2 < k1; ++k2)
		{
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];		
		  Tmp += this->PartialZPrecalculatedHamiltonian[k1][k2];	      
		  TmpTotal += Tmp * vSource[Index2];
		  vDestination[Index2] += Tmp * TmpVSource;
		  ++Index2;
		}     
	    }
	  // diagonal
	  TmpTotal += this->DiagonalElements[Index1] * vSource[Index1];
	  // second part	  
	  Index2 = LastComponent;
	  ReducedIndex2 = ReducedLast;	  
	  k2 = kLast;	 		
	  if (kLast != 0)
	    {
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex2][ReducedIndex1];
	      if (ReducedIndex1 == ReducedIndex2)
		for (; k2 < this->NbrStateZ; ++k2)
		  {		  	       
		    TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		    Tmp = 0.0;
		    for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		
		      Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	    
		    Tmp += PartialZPrecalculatedHamiltonian[k1][k2];
		    TmpTotal += Tmp * vSource[Index2];	    
		    ++Index2;
		  }	
	      else
		for (; k2 < this->NbrStateZ; ++k2)
		  {		  	       
		    TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		    Tmp = 0.0;
		    for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		
		      Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];	    
		    TmpTotal += Tmp * vSource[Index2];	    
		    ++Index2;
		  }
	      ++ReducedIndex2;	    
	    }
	  for (; ReducedIndex2 < ReducedDim; ++ReducedIndex2)
	    {
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[ReducedIndex2][ReducedIndex1];
	      for (k2 = 0; k2 < this->NbrStateZ; ++k2)
		{
		  TmpWaveFunctionOverlapZ = this->WaveFunctionOverlapZ[k1][k2];
		  Tmp = 0.0;
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
		    Tmp += TmpWaveFunctionOverlapZ[CellZ] * TmpPrecalculatedHamiltonian[CellZ];
		  TmpTotal += Tmp * vSource[Index2];
		  ++Index2;
		}
	    }	  
	  vDestination[Index1] += TmpTotal;
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
  if (this->NbrPrecalculatedDimension == 0x100)
    {
      double Tmp = 0.0; int Index2 = 0;
      for (int Index1 = firstComponent; Index1 < LastComponent; ++Index1)
	{
	  Tmp = (this->DiagonalElements[Index1] * vSource[Index1]);	  
	  Index2 = 0;
	  for (; Index2 < Index1; ++Index2)
	    Tmp += this->FullPrecalculatedHamiltonian[Index1][Index2] * vSource[Index2];	  
	  ++Index2;
	  for (; Index2 < Dim; ++Index2)
	    Tmp += this->FullPrecalculatedHamiltonian[Index2][Index1] * vSource[Index2];
	  vDestination[Index1] += Tmp;
	}
      return vDestination;
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& QuantumDots3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& QuantumDots3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							  int firstComponent, int nbrComponent)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& QuantumDots3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
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

ComplexVector& QuantumDots3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							     int firstComponent, int nbrComponent)
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> QuantumDots3DHamiltonian::LeftInteractionOperators()  
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators 
//
// return value = list of right interaction operators

List<Matrix*> QuantumDots3DHamiltonian::RightInteractionOperators()  
{
  List<Matrix*> TmpList;
  return TmpList;
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, QuantumDots3DHamiltonian& H)
{
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, QuantumDots3DHamiltonian& H)
{
  return Str;
}

// evaluate all interaction factors
//   
// memory = amount of memory available to store precalculated values

void QuantumDots3DHamiltonian::EvaluateInteractionFactors(int memory)
{
  int UsedMemory = 0;
  this->WaveFunctionOverlapX = this->EvaluateWaveFunctionOverlap (this->XSize, this->NbrCellX, this->NbrStateX, UsedMemory);
  this->WaveFunctionOverlapY = this->EvaluateWaveFunctionOverlap (this->YSize, this->NbrCellY, this->NbrStateY, UsedMemory);
  this->WaveFunctionOverlapZ = this->EvaluateWaveFunctionZOverlap (UsedMemory);

  UsedMemory -= this->Space->GetHilbertSpaceDimension () * sizeof(double);
  double InvXFactor = HAMILTONIAN_FACTOR / (this->Mux * this->XSize * this->XSize);
  double InvYFactor = HAMILTONIAN_FACTOR / (this->Muy * this->YSize * this->YSize);
  double InvZFactor = HAMILTONIAN_FACTOR / (this->Muz * this->ZSize * this->ZSize);
  int TotalIndex = 0;
  double FactorX = 0.0;
  double FactorY = 0.0;
  double Factor;
  double TmpElement = 0.0;
  int IncNbrCellZ = this->NbrCellZ + 1;
  this->DiagonalElements = new double[ this->Space->GetHilbertSpaceDimension ()];
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      FactorX = ((double) ((i + 1) * (i + 1))) * InvXFactor;
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  FactorY = ((double) ((j + 1) * (j + 1))) * InvYFactor + FactorX;
	  for (int k = 0; k < this->NbrStateZ; ++k)
	    {
	      TmpElement = FactorY + ((double) ((k + 1) * (k + 1))) * InvZFactor  + (this->PartialZPrecalculatedHamiltonian[k][k]);
	      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		    {
		      Factor = this->WaveFunctionOverlapZ[k][k][CellZ] * this->WaveFunctionOverlapY[j][j][CellY];
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)			
			TmpElement += this->InteractionFactors[CellZ][CellY][CellX] * Factor * this->WaveFunctionOverlapX[i][i][CellX];			
		    }
		}
	      this->DiagonalElements[TotalIndex] = TmpElement;
	      ++TotalIndex;
	    }
	}
    }
  memory -= UsedMemory;
  if (memory < 0)
    {
      cout << "Calculation from scratch is done" << endl;
      this->NbrPrecalculatedDimension = 0;
      return;
    }

  if ( (double)memory > ((double(this->NbrStateX * this->NbrStateX) * double(this->NbrStateY * this->NbrStateY) * double(this->NbrStateZ * this->NbrStateZ)) * (double(sizeof(double)))))
    {      
      cout << "Calculation from fully calculated Hamiltonian" << endl << endl;
      this->NbrPrecalculatedDimension = 0x100;
      int Dim = this->NbrStateX * this->NbrStateY * this->NbrStateZ;
      this->FullPrecalculatedHamiltonian = new double* [Dim];      
      double**** Inter1 = new double*** [this->NbrStateX];   
      double tmp = 0.0;
      
      for (int m1 = 0; m1 < this->NbrStateX; ++m1)
	{
	  Inter1[m1] = new double** [this->NbrStateX];
	  for (int m2 = 0; m2 < m1; ++m2)
	    {
	      Inter1[m1][m2] = new double* [this->NbrCellZ];
	      Inter1[m2][m1] = new double* [this->NbrCellZ];
	      for (int k = 0; k < this->NbrCellZ; ++k)
		{
		  Inter1[m1][m2][k] = new double [this->NbrCellY];
		  Inter1[m2][m1][k] = new double [this->NbrCellY];
		  for (int j = 0; j < this->NbrCellY; ++j)
		    {		     
		      tmp = 0.0;
		      for (int i = 0; i < this->NbrCellX; ++i)			
			tmp += this->InteractionFactors[k][j][i] * this->WaveFunctionOverlapX[m1][m2][i];
		      Inter1[m1][m2][k][j] = tmp;		     
		      Inter1[m2][m1][k][j] = tmp;
		    }
		}	      
	    }
	  Inter1[m1][m1] = new double* [this->NbrCellZ];
	  for (int k = 0; k < this->NbrCellZ; ++k)
	    {
	      Inter1[m1][m1][k] = new double [this->NbrCellY];
	      for (int j = 0; j < this->NbrCellY; ++j)
		{		     
		  tmp = 0.0;
		  for (int i = 0; i < this->NbrCellX; ++i)			
		    tmp += this->InteractionFactors[k][j][i] * this->WaveFunctionOverlapX[m1][m1][i];			
		  Inter1[m1][m1][k][j] = tmp;		     
		}
	    }	      
	}
      double***** Inter2 = new double**** [this->NbrStateY];
      for (int n1 = 0; n1 < this->NbrStateY; ++n1)
	{
	  Inter2[n1] = new double*** [this->NbrStateY];
	  for (int n2 = 0; n2 < n1; ++n2)
	    {
	      Inter2[n1][n2] = new double** [this->NbrStateX];
	      Inter2[n2][n1] = new double** [this->NbrStateX];
	      for (int m1 = 0; m1 < this->NbrStateX; ++m1)
		{
		  Inter2[n1][n2][m1] = new double* [this->NbrStateX];
		  Inter2[n2][n1][m1] = new double* [this->NbrStateX];
		  for(int m2 = 0; m2 < this->NbrStateX; ++m2)
		    {
		      Inter2[n1][n2][m1][m2] = new double [this->NbrCellZ];
		      Inter2[n2][n1][m1][m2] = new double [this->NbrCellZ];
		      for (int k = 0; k < this->NbrCellZ; ++k)
			{
			  tmp = 0.0;
			  for (int j = 0; j < this->NbrCellY; ++j)
			    tmp += this->WaveFunctionOverlapY[n1][n2][j] * Inter1[m1][m2][k][j];
			  Inter2[n1][n2][m1][m2][k] = tmp;			  
			  Inter2[n2][n1][m1][m2][k] = tmp;
			}
		    }
		}
	    }
	  Inter2[n1][n1] = new double** [this->NbrStateX];
	  for (int m1 = 0; m1 < this->NbrStateX; ++m1)
	    {
	      Inter2[n1][n1][m1] = new double* [this->NbrStateX];
	      for(int m2 = 0; m2 < this->NbrStateX; ++m2)
		{
		  Inter2[n1][n1][m1][m2] = new double [this->NbrCellZ];
		  for (int k = 0; k < this->NbrCellZ; ++k)
		    {
		      tmp = 0.0;
		      for (int j = 0; j < this->NbrCellY; ++j)
			tmp += this->WaveFunctionOverlapY[n1][n1][j] * Inter1[m1][m2][k][j];
		      Inter2[n1][n1][m1][m2][k] = tmp;			  
		    }
		}
	    }
	} 
      
      int ReducedIndex1 = 0, ReducedIndex2 = 0, m1 = 0, m2 = 0, n1 = 0, n2 = 0, p1 = 0, p2 = 0;
   
      for (int Index1 = 0; Index1 < Dim; ++Index1) 
	{
	  this->FullPrecalculatedHamiltonian[Index1] = new double [Index1];
	  ReducedIndex1 = Index1 / this->NbrStateZ;
	  p1 = Index1 - ReducedIndex1 * this->NbrStateZ;
	  m1 = ReducedIndex1 / this->NbrStateY;
	  n1 = ReducedIndex1 - m1 * this->NbrStateY;	  
	  for (int Index2 = 0; Index2 < Index1; ++Index2) 
	    {	    
	      ReducedIndex2 = Index2 / this->NbrStateZ;
	      p2 = Index2 - ReducedIndex2 * this->NbrStateZ;
	      m2 = ReducedIndex2 / this->NbrStateY;
	      n2 = ReducedIndex2 - m2 * this->NbrStateY;	  	      
	      tmp = 0.0;
	      for (int k = 0; k < this->NbrCellZ; ++k)
		tmp += this->WaveFunctionOverlapZ[p1][p2][k] * Inter2[n1][n2][m1][m2][k]; 
	      if (ReducedIndex1 == ReducedIndex2)
		tmp += this->PartialZPrecalculatedHamiltonian[p1][p2];
	      this->FullPrecalculatedHamiltonian[Index1][Index2] = tmp;
	    }
	    
	}
      delete[] Inter1; delete[] Inter2;
      return; 
    }
/*
  if (((double) memory) > (((double) (this->NbrStateX * this->NbrStateY)) * 
			    ((double) (this->NbrStateX * this->NbrStateY))) * 
			   (((double) sizeof(double*)) + ((double) this->NbrCellZ) * ((double) sizeof(double))))
    {
      cout << "Calculation with all precalculated Hamiltonian is done" << endl;
      this->NbrPrecalculatedDimension = 0x010;
      int Dim = this->NbrStateX * this->NbrStateY;
      this->Partial2DPrecalculatedHamiltonian = new double** [Dim];
      this->Partial2DDiagonalPrecalculatedHamiltonian = new double* [Dim];
      double Tmp;
      double Tmp2;
      double* TmpWaveFunctionOverlapX;
      double* TmpWaveFunctionOverlapY;
      double* TmpPrecalculatedHamiltonian1;
      double* TmpPrecalculatedHamiltonian2;
      int j1 = 0;
      int i1 = 0;
      for (int Index1 = 0; Index1 < Dim; ++Index1)
	{
	  this->Partial2DPrecalculatedHamiltonian[Index1] = new double* [Dim];
	  int Index2 = 0;
	  int i2 = 0;
	  int j2 = 0;
	  while (Index2 < Index1)
	    {
	      if (i2 <= i1)
		TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[i1][i2];
	      else
		TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[i2][i1];
	      if (j2 <= j1)
		TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[j1][j2];
	      else
		TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[j2][j1];
	      this->Partial2DPrecalculatedHamiltonian[Index1][Index2] = new double [this->NbrCellZ];
	      TmpPrecalculatedHamiltonian1 = this->Partial2DPrecalculatedHamiltonian[Index1][Index2];
	      this->Partial2DPrecalculatedHamiltonian[Index2][Index1] = new double [this->NbrCellZ];
	      TmpPrecalculatedHamiltonian2 = this->Partial2DPrecalculatedHamiltonian[Index2][Index1];	       

	      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  Tmp = 0.0;
		  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		    {
		      Tmp2 = TmpWaveFunctionOverlapY[CellY];
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)			
			Tmp += this->InteractionFactors[CellZ][CellY][CellX] * TmpWaveFunctionOverlapX[CellX] * Tmp2;
		    } 
		  TmpPrecalculatedHamiltonian1[CellZ] = Tmp;  		 
		  TmpPrecalculatedHamiltonian2[CellZ] = Tmp;
		}    
	      ++Index2;
	      ++j2;
	      if (j2 == this->NbrStateY)
		{
		  j2 = 0;
		  ++i2;
		}
	    }

	  if (i2 <= i1)
	    TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[i1][i2];
	  else
	    TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[i2][i1];
	  if (j2 <= j1)
	    TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[j1][j2];
	  else
	    TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[j2][j1];
	  this->Partial2DPrecalculatedHamiltonian[Index1][Index2] = new double [this->NbrCellZ];
	  TmpPrecalculatedHamiltonian1 = this->Partial2DPrecalculatedHamiltonian[Index1][Index1];
	  
	  for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
	    {
	      Tmp = 0.0;
	      for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		{
		  Tmp2 = TmpWaveFunctionOverlapY[CellY];
		  for (int CellX = 0; CellX < this->NbrCellX; ++CellX)			
		    Tmp += this->InteractionFactors[CellZ][CellY][CellX] * TmpWaveFunctionOverlapX[CellX] * Tmp2;
		} 
	      TmpPrecalculatedHamiltonian1[CellZ] = Tmp;  
	    } 
	  
	  ++j1;
	  if (j1 == this->NbrStateY)
	    {
	      j1 = 0;
	      ++i1;
	    }
	}   
      return;
    }
*/
  if (((double) memory) > ((0.5 * ((double) (this->NbrStateX * this->NbrStateY)) * 
			    ((double) (this->NbrStateX * this->NbrStateY))) * 
			   (((double) sizeof(double*)) + ((double) this->NbrCellZ) * ((double) sizeof(double)))))
    {
      cout << "Calculation with symmetric precalculated Hamiltonian is done" << endl;
      this->NbrPrecalculatedDimension = 0x001;
      int Dim = this->NbrStateX * this->NbrStateY;
      this->Partial2DPrecalculatedHamiltonian = new double** [Dim];
      this->Partial2DDiagonalPrecalculatedHamiltonian = new double* [Dim];
      double Tmp;
      double Tmp2;
      double* TmpWaveFunctionOverlapX;
      double* TmpWaveFunctionOverlapY;
      double* TmpPrecalculatedHamiltonian;
      int j1 = 0;
      int i1 = 0;
      for (int Index1 = 0; Index1 < Dim; ++Index1)
	{
	  this->Partial2DPrecalculatedHamiltonian[Index1] = new double* [Index1 + 1];
	  int Index2 = 0;
	  int i2 = 0;
	  int j2 = 0;
	  while (Index2 <= Index1)
	    {
	      if (i2 <= i1)
		TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[i1][i2];
	      else
		TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[i2][i1];
	      if (j2 <= j1)
		TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[j1][j2];
	      else
		TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[j2][j1];
	      this->Partial2DPrecalculatedHamiltonian[Index1][Index2] = new double [this->NbrCellZ];
	      TmpPrecalculatedHamiltonian = this->Partial2DPrecalculatedHamiltonian[Index1][Index2];
	      for (int CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		{
		  Tmp = 0.0;
		  for (int CellY = 0; CellY < this->NbrCellY; ++CellY)
		    {
		      Tmp2 = TmpWaveFunctionOverlapY[CellY];
		      for (int CellX = 0; CellX < this->NbrCellX; ++CellX)			
			Tmp += this->InteractionFactors[CellZ][CellY][CellX] * TmpWaveFunctionOverlapX[CellX] * Tmp2;
		    } 
		  TmpPrecalculatedHamiltonian[CellZ] = Tmp;  
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
      return;
    }

}

// evaluate wave function overlaps on a cell in a given direction
//
// size = system length in the choosen direction
// nbrStep = number of subdivision in the choosen direction
// nbrState = number of state in the choosen direction
// memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
// return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)

double*** QuantumDots3DHamiltonian::EvaluateWaveFunctionOverlap(double size, int nbrStep, int nbrState, int& memory)
{
  memory += ((((nbrState * (nbrState + 1)) >> 1) * (nbrStep * sizeof(double) + sizeof(double*)))
	     + nbrState * sizeof(double**));
  double*** TmpArray = new double** [nbrState];
  double StepInc = 1.0 / ((double) nbrStep);
  double Tmp;
  double Diff;
  for (int i = 0; i < nbrState; ++i)
    {
      TmpArray[i] = new double* [i + 1];
      for (int j = 0; j < i; ++j)
	{
	  TmpArray[i][j] = new double [nbrStep];
	  for (int k = 0; k < nbrStep; ++k)
	    {
	      Diff = (double) (i - j);
	      Tmp = M_PI * Diff * StepInc;	      
	      TmpArray[i][j][k] = M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / Diff;
	      Diff = (double) (i + j + 2);
	      Tmp = M_PI * Diff * StepInc;	      
	      TmpArray[i][j][k] -= M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / Diff;
	    }
	}
      TmpArray[i][i] = new double [nbrStep];
      for (int k = 0; k < nbrStep; ++k)
	{
	  Tmp = M_PI * (double) (2 * i + 2) * StepInc;	      
	  TmpArray[i][i][k] = StepInc - M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / ((double) (2 * i + 2));
	}     
    }
  return TmpArray;
}

// evaluate wave function overlaps on a cell in the z direction
//
// memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
// return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)

double*** QuantumDots3DHamiltonian::EvaluateWaveFunctionZOverlap(int& memory)
{
  memory += ((((this->NbrStateZ * (this->NbrStateZ + 1)) >> 1) * ((this->NbrCellZ)  * sizeof(double) + sizeof(double*)))
	     + this->NbrStateZ * sizeof(double**));
  double*** TmpArray = new double** [this->NbrStateZ];
  this->PartialZPrecalculatedHamiltonian = new double*[this->NbrStateZ];
  double LeftRegionSize = 0.0;
  for (int k = 0; k < this->LeftNumber; ++k)
    LeftRegionSize += PreConstantRegionSize[k];
  double RightRegionSize = 0.0;
  for (int k = 0; k < this->RightNumber; ++k)
    RightRegionSize += PostConstantRegionSize[k];

  double StepInc = (this->ZSize - LeftRegionSize - RightRegionSize) / ((((double) this->NbrCellZ)) * this->ZSize);
  double Shift = LeftRegionSize / this->ZSize;
  double PostShift = this->ZSize - RightRegionSize;
  double Tmp;
  double TmpShift;
  double Diff;
  double Tmp2;
  double TmpShift2;
  double Diff2;
  double TmpLength1; double TmpLength2; double Tmp3; double Tmp4; double Tmp5;
  for (int i = 0; i < this->NbrStateZ; ++i)
    {
      TmpArray[i] = new double* [this->NbrStateZ];
      this->PartialZPrecalculatedHamiltonian[i] = new double [this->NbrStateZ];
      for (int j = 0; j < i; ++j)
	{
	  Diff = (double) (i - j);
	  Tmp = M_PI * Diff;
	  TmpShift = Tmp * Shift;
	  TmpArray[i][j] = new double [this->NbrCellZ];
	  TmpArray[j][i] = new double [this->NbrCellZ];
	  Diff2 = (double) (i + j + 2);
	  Tmp2 = M_PI * Diff2;
	  TmpShift2 = Tmp2 * Shift;
	  Diff = M_1_PI / Diff;
	  Diff2 = M_1_PI / Diff2;
	  TmpLength1 = 0.0; TmpLength2 = 0.0;
	  Tmp3 = Tmp/this->ZSize;
	  Tmp4 = Tmp2/this->ZSize;
	  Tmp5 = 0.0;
	  for (int k = 0; k < this->LeftNumber; ++k)
	    {
	      TmpLength1 += PreConstantRegionSize[k];
	      Tmp5 +=  ((sin(Tmp3 * TmpLength1) - sin(Tmp3 * TmpLength2)) * Diff
		                  - (sin(Tmp4 * TmpLength1) - sin(Tmp4 * TmpLength2)) * Diff2) * (this->PreConstantRegionPotential[k]);
	      TmpLength2 += PreConstantRegionSize[k];
	    }

	  Tmp *= StepInc;
	  Tmp2 *= StepInc;
	  for (int k = 0; k < this->NbrCellZ; ++k)
	    {
	      TmpArray[i][j][k] =  ((sin ((Tmp * ((double) (k + 1))) + TmpShift) - sin ((Tmp * ((double) k)) + TmpShift)) * Diff
					     - (sin ((Tmp2 * ((double) (k + 1))) + TmpShift2) - sin ((Tmp2 * ((double) k)) + TmpShift2)) * Diff2);
	      TmpArray[j][i][k] =  TmpArray[i][j][k];
	    }
	  TmpLength1 = PostShift;
	  TmpLength2 = PostShift;
	  for (int k = 0; k < this->RightNumber; ++k)
	    {
	      TmpLength1 += PostConstantRegionSize[k];
	      Tmp5 +=  ((sin(Tmp3 * TmpLength1) - sin(Tmp3 * TmpLength2)) * Diff
		     - (sin(Tmp4 * TmpLength1) - sin(Tmp4 * TmpLength2)) * Diff2) * (this->PostConstantRegionPotential[k]);
	      TmpLength2 += PostConstantRegionSize[k];
	    }
	  this->PartialZPrecalculatedHamiltonian[i][j] = Tmp5;
	  this->PartialZPrecalculatedHamiltonian[j][i] = Tmp5;
	}

      TmpArray[i][i] = new double [this->NbrCellZ];
      Diff = (double) (2 * i + 2);
      Tmp = M_PI * Diff;
      TmpShift = Tmp * Shift;
      Diff = M_1_PI / Diff;
      Tmp3 = Tmp/this->ZSize;
      TmpLength1 = 0.0; TmpLength2 = 0.0; Tmp5 = 0.0;
      for (int k = 0; k < this->LeftNumber; ++k)
	{
	  TmpLength1 += PreConstantRegionSize[k];
	  Tmp5 += ((PreConstantRegionSize[k]/this->ZSize) - (sin(Tmp3 * TmpLength1) - sin(Tmp3 * TmpLength2)) * Diff) * (this->PreConstantRegionPotential[k]);
	  TmpLength2 += PreConstantRegionSize[k];
	}
      TmpLength1 = PostShift; TmpLength2 = PostShift;
      for (int k = 0; k < this->RightNumber; ++k)
	{
	  TmpLength1 += PostConstantRegionSize[k];
	  Tmp5 += ((PostConstantRegionSize[k]/this->ZSize) - (sin(Tmp3 * TmpLength1) - sin(Tmp3 * TmpLength2)) * Diff) * (this->PostConstantRegionPotential[k]);
	  TmpLength2 += PostConstantRegionSize[k];
	}
      this->PartialZPrecalculatedHamiltonian[i][i] = Tmp5;

      Tmp *= StepInc;
      for (int k = 0; k < this->NbrCellZ; ++k)
	{
	  TmpArray[i][i][k] = StepInc - (sin ((Tmp * ((double) (k + 1))) + TmpShift) - sin ((Tmp * ((double) k)) + TmpShift)) * Diff;
	}
    }
  return TmpArray;
}

