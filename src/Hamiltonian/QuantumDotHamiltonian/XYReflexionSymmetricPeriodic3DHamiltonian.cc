////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2004 Duc-Phuong Nguyen                    //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 25/03/2004                        //
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
#include "Hamiltonian/QuantumDotHamiltonian/XYReflexionSymmetricPeriodic3DHamiltonian.h"
#include "Complex.h"
#include "Vector/ComplexVector.h"
#include "Tools/QuantumDot/Potential/ThreeDConstantCellPotential.h"
#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>

using std::ostream;
using std::cout;
using std::endl;


#define PERIODIC_HAMILTONIAN_FACTOR 150.4


// constructor from data
//
// space = Hilbert space
// pairX = whether basis is pair in X direction, if not impair
// pairY = whether basis is pair in Y direction, if not impair
// xSize = the sample length in X direction
// ySize = the sample length in Y direction
// zSize = the sample length in Z direction
// mux = effective mass in X direction
// muy = effective mass in Y direction
// muz = effective mass in Z direction
// nbrCellX = number of steps in X direction
// nbrCellY = number of steps in Y direction
// nbrCellZ = number of steps in Z direction
// PotentielInput = pointer to a 3D potential with constant value in a cell

XYReflexionSymmetricPeriodic3DHamiltonian::XYReflexionSymmetricPeriodic3DHamiltonian(XYReflexionSymmetricPeriodic3DOneParticle* space, bool pairX, bool pairY, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDConstantCellPotential* PotentialInput)
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
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();

  this->InteractionFactors = new double** [this->NbrCellZ];      
  int CenterX = this->NbrCellX / 2; int CenterY = this->NbrCellY / 2;
  for (int k = 0; k < this->NbrCellZ; ++k)
    {
      this->InteractionFactors[k] = new double* [this->NbrCellY];
      for (int j = 0; j < this->NbrCellY; ++j)
	{
	  this->InteractionFactors[k][j] = new double [this->NbrCellX];
	  for (int i = 0; i < this->NbrCellX; ++i)
	    {
	      this->InteractionFactors[k][j][i] = PotentialInput->GetPotential(i, j, k);	
	      
	      if ((i > CenterX) && (j > CenterY))
	        if ((this->InteractionFactors[k][j][i] != this->InteractionFactors[k][this->NbrCellY - j][i]) || (this->InteractionFactors[k][j][i] != this->InteractionFactors[k][j][this->NbrCellX - i]) ||
		    (this->InteractionFactors[k][j][i] != this->InteractionFactors[k][this->NbrCellY - j][this->NbrCellX - i]))
		  {
		    cout << "The potential is not reflexion symmetric in X or Y direction. Exit now!" << endl;
		    exit(0);
		  }
	      
	    }
	}	        
    }
  cout << "Hamiltonian dimension: " << this->Space->GetHilbertSpaceDimension () << endl;
  cout << "Evaluation of Hamiltionian elements ..." << endl;
  this->EvaluateInteractionFactors(pairX, pairY);
  cout << "Evaluation finished ..." << endl;
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

XYReflexionSymmetricPeriodic3DHamiltonian::XYReflexionSymmetricPeriodic3DHamiltonian(const XYReflexionSymmetricPeriodic3DHamiltonian& hamiltonian)
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
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();
  this->KineticElements = hamiltonian.KineticElements;
  this->InteractionFactors = hamiltonian.InteractionFactors;
  this->WaveFunctionOverlapX = hamiltonian.WaveFunctionOverlapX;
  this->WaveFunctionOverlapY = hamiltonian.WaveFunctionOverlapY;
  this->RealWaveFunctionOverlapZ = hamiltonian.RealWaveFunctionOverlapZ;
  this->ImaginaryWaveFunctionOverlapZ = hamiltonian.ImaginaryWaveFunctionOverlapZ;
  this->RealHamiltonian =  hamiltonian.RealHamiltonian;
  this->ImaginaryHamiltonian = hamiltonian.ImaginaryHamiltonian;
}

// destructor
//

XYReflexionSymmetricPeriodic3DHamiltonian::~ XYReflexionSymmetricPeriodic3DHamiltonian()
{  
  delete[] this->KineticElements;
  delete[] this->InteractionFactors;
  delete[] this->WaveFunctionOverlapX;
  delete[] this->WaveFunctionOverlapY;
  delete[] this->RealWaveFunctionOverlapZ;
  delete[] this->ImaginaryWaveFunctionOverlapZ;
  delete[] this->RealHamiltonian;
  delete[] this->ImaginaryHamiltonian;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* XYReflexionSymmetricPeriodic3DHamiltonian::Clone ()
{
  return new XYReflexionSymmetricPeriodic3DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void XYReflexionSymmetricPeriodic3DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void XYReflexionSymmetricPeriodic3DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->KineticElements[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex XYReflexionSymmetricPeriodic3DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
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

Complex XYReflexionSymmetricPeriodic3DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& XYReflexionSymmetricPeriodic3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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

ComplexVector& XYReflexionSymmetricPeriodic3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
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

ComplexVector& XYReflexionSymmetricPeriodic3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  int tmpDim = this->NbrStateY * this->NbrStateZ;
  int** RI = new int* [this->NbrStateX];
  int m1 = 0, n1 = 0, p1 = 0, m2 = 0, n2 = 0, p2 = 0;
  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {
      RI[m1] = new int [this->NbrStateY];
      for (n1 = 0; n1 < this->NbrStateY; ++n1)
	RI[m1][n1] = m1 * tmpDim + n1 * this->NbrStateZ;
    }
  int OriginZ = this->NbrStateZ - 1;
  int Index1 = 0, Index2 = 0;
  double TmpRe = 0.0, TmpIm = 0.0;
  int IndexZ = 0, tmpIndexZ = 0;
  double* TmpRealHamiltonian; double* TmpImaginaryHamiltonian;

  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {
      for (n1 = 0; n1 < this->NbrStateY; ++n1)
	{
	  Index1 = RI[m1][n1]; // Index1 = (m1 * N + n1) * H
	  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
	    {
	      for (m2 = 0; m2 < this->NbrStateX; ++m2)
	        {
	          for (n2 = 0; n2 < this->NbrStateY; ++n2)
		    {
		      Index2 = RI[m2][n2]; // Index2 = (m2 * N + n2) * H
		      for (p2 = 0; p2 < this->NbrStateZ; ++p2)
			{
			  IndexZ = -p1 + p2 + OriginZ;
			  vDestination.Re(Index1) += (this->RealHamiltonian[m1][n1][m2][n2][IndexZ] * vSource.Re(Index2) - this->ImaginaryHamiltonian[m1][n1][m2][n2][IndexZ] * vSource.Im(Index2));
			  vDestination.Im(Index1) += (this->RealHamiltonian[m1][n1][m2][n2][IndexZ] * vSource.Im(Index2) + this->ImaginaryHamiltonian[m1][n1][m2][n2][IndexZ] * vSource.Re(Index2));
			  ++Index2;
			}
		    }
		}
	      vDestination.Re(Index1) += (this->KineticElements[Index1] * vSource.Re(Index1));
	      vDestination.Im(Index1) += (this->KineticElements[Index1] * vSource.Im(Index1));
	      ++Index1;
	    }
	}
    }
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

ComplexVector& XYReflexionSymmetricPeriodic3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  if ((firstComponent == 0) && (nbrComponent == this->Space->GetHilbertSpaceDimension()))
    return this->LowLevelAddMultiply(vSource, vDestination);
  else
    {
      int lastComponent = firstComponent + nbrComponent;
      int m1, m2, n1, n2, p1, p2;
      m1 = firstComponent / (this->NbrStateY * this->NbrStateZ);
      int tmpIndex = firstComponent - m1 * this->NbrStateY * this->NbrStateZ;
      n1 = tmpIndex / this->NbrStateZ;
      p1 = tmpIndex - n1 * this->NbrStateZ;
      int Index1 = firstComponent;
      /*
      int ReducedIndex1 = Index1 / this->NbrStateZ;
      p1 = Index1 - ReducedIndex1 * this->NbrStateZ;
      m1 = ReducedIndex1 / this->NbrStateY;
      n1 = ReducedIndex1 - m1 * this->NbrStateY;*/
      double TmpRe = 0.0, TmpIm = 0.0;      
      int Index2 = 0; int IndexZ = 0; int tmpIndexZ = 0;
      double* TmpRealHamiltonian; double* TmpImaginaryHamiltonian;

      for (; Index1 < lastComponent; ++Index1)
	{
	  Index2 = 0;
	  TmpRe = 0.0; TmpIm = 0.0;
	  for (m2 = 0; m2 < m1; ++m2)
	    {
	      for (n2 = 0; n2 < n1; ++n2)
		{
		  TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m2][n2];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m2][n2];
		  //for (p2 = 0; p2 < p1; ++p2)
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {		      
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  tmpIndexZ = this->NbrStateZ - p1;
		  //for (p2 = p1; p2 < this->NbrStateZ; ++p2)
		  for (IndexZ = 0; IndexZ < tmpIndexZ; ++IndexZ)
		    {		      
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		}
	      for (n2 = n1; n2 < this->NbrStateY; ++n2)
		{
		  TmpRealHamiltonian = this->RealHamiltonian[m1][n2][m2][n1];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n2][m2][n1];
		  //for (p2 = 0; p2 < p1; ++p2)
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {		      
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  tmpIndexZ = this->NbrStateZ - p1;
		  //for (p2 = p1; p2 < this->NbrStateZ; ++p2)
		  for (IndexZ = 0; IndexZ < tmpIndexZ; ++IndexZ)
		    {		      
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		}
	    }
	  for (m2 = m1; m2 < this->NbrStateX; ++m2)
	    {
	      for (n2 = 0; n2 < n1; ++n2)
		{
		  TmpRealHamiltonian = this->RealHamiltonian[m2][n1][m1][n2];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n1][m1][n2];
		  //for (p2 = 0; p2 < p1; ++p2)
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {		      
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  tmpIndexZ = this->NbrStateZ - p1;
		  //for (p2 = p1; p2 < this->NbrStateZ; ++p2)
		  for (IndexZ = 0; IndexZ < tmpIndexZ; ++IndexZ)
		    {		      
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		}
	      for (n2 = n1; n2 < this->NbrStateY; ++n2)
		{
		  TmpRealHamiltonian = this->RealHamiltonian[m2][n2][m1][n1];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n2][m1][n1];
		  //for (p2 = 0; p2 < p1; ++p2)
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {		      
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  tmpIndexZ = this->NbrStateZ - p1;
		  //for (p2 = p1; p2 < this->NbrStateZ; ++p2)
		  for (IndexZ = 0; IndexZ < tmpIndexZ; ++IndexZ)
		    {		      
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		}
	    }

	  TmpRe += KineticElements[Index1] * vSource.Re(Index1);
	  TmpIm += KineticElements[Index1] * vSource.Im(Index1);
	  
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
}

// evaluate all interaction factors
// 
// pairX = whether basis is pair in X direction, if not impair
// pairY = whether basis is pair in Y direction, if not impair

void XYReflexionSymmetricPeriodic3DHamiltonian::EvaluateInteractionFactors(bool pairX, bool pairY)
{
  if (pairX)    
    this->WaveFunctionOverlapX = this->EvaluateCosinusWaveFunctionOverlap(this->XSize, this->NbrCellX, this->NbrStateX);
  else
    this->WaveFunctionOverlapX = this->EvaluateSinusWaveFunctionOverlap(this->XSize, this->NbrCellX, this->NbrStateX);
  if (pairY)
    this->WaveFunctionOverlapY = this->EvaluateCosinusWaveFunctionOverlap(this->YSize, this->NbrCellY, this->NbrStateY);
  else
    this->WaveFunctionOverlapY = this->EvaluateSinusWaveFunctionOverlap(this->YSize, this->NbrCellY, this->NbrStateY);

  if (!this->EvaluatePlaneWaveFunctionOverlap(this->NbrCellZ, this->NbrStateZ, this->RealWaveFunctionOverlapZ, this->ImaginaryWaveFunctionOverlapZ))
    {
      cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;
      exit(0);
    }

  double InvXFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Mux * this->XSize * this->XSize);
  double InvYFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muy * this->YSize * this->YSize);
  double InvZFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muz * this->ZSize * this->ZSize);
  
  this->KineticElements = new double[this->Space->GetHilbertSpaceDimension()];
  // this->NbrStateX * this->NbrStateY * this->NbrStateZ this->Space->GetHilbertSpaceDimension()
  double FactorX = 0.0, FactorY = 0.0;
  int TotalIndex = 0;
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      if (pairX)
	FactorX = double(i * i) * InvXFactor;
      else
	FactorX = double((i + 1) * (i + 1)) * InvXFactor;
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  if (pairY)
	    FactorY = double(j * j) * InvYFactor + FactorX;
	  else
	    FactorY = double((j + 1) * (j + 1)) * InvYFactor + FactorX;
	  for (int k = 0; k < this->NbrStateZ; ++k)
	    {
	      this->KineticElements[TotalIndex] = FactorY + double((k + this->LowerImpulsionZ) * (k + this->LowerImpulsionZ)) * InvZFactor;
	      ++TotalIndex;
	    }
	}
    }

  double* TmpRealHamiltonian;
  double* TmpImaginaryHamiltonian;
  double* TmpWaveFunctionOverlapX;
  double* TmpWaveFunctionOverlapY;
  double* TmpRealWaveFunctionOverlapZ;
  double* TmpImaginaryWaveFunctionOverlapZ;
  int m1 = 0, m2 = 0, n1 = 0, n2 = 0, p = 0;
  double TmpRe = 0.0, TmpIm = 0.0, TmpXY = 0.0, TmpX = 0.0, TmpY = 0.0;
  int CellX = 0, CellY = 0, CellZ = 0;

  this->RealHamiltonian = new double**** [this->NbrStateX];
  this->ImaginaryHamiltonian = new double**** [this->NbrStateX];
  double* TmpInter = new double [this->NbrCellZ];
  int LengthZ = (this->NbrStateZ - 1) * 2 + 1;
  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {	  
      this->RealHamiltonian[m1] = new double*** [this->NbrStateY];	
      this->ImaginaryHamiltonian[m1] = new double*** [this->NbrStateY];	
      for (n1 = 0; n1 < this->NbrStateY; ++n1)
	{
	  this->RealHamiltonian[m1][n1] = new double** [this->NbrStateX];
	  this->ImaginaryHamiltonian[m1][n1] = new double** [this->NbrStateX];	      	      
	  for (m2 = 0; m2 < this->NbrStateX; ++m2)
	    {	      
	      this->RealHamiltonian[m1][n1][m2] = new double* [this->NbrStateY];	
 	      this->ImaginaryHamiltonian[m1][n1][m2] = new double* [this->NbrStateY];
	      TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[m1][m2]; 
	      for (n2 = 0; n2 < this->NbrStateY; ++n2)
		{
		  TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[n1][n2];	  
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		    {
		      TmpXY = 0.0;
		      for (CellY = 0; CellY < this->NbrCellY; ++CellY)
			{
			  TmpY = TmpWaveFunctionOverlapY[CellY];
			  TmpX = 0.0;
			  for (CellX = 0; CellX < this->NbrCellX; ++CellX)			     
			    TmpX += this->InteractionFactors[CellZ][CellY][CellX] * TmpWaveFunctionOverlapX[CellX];
			  TmpXY += (TmpX * TmpY);
			} 		      
		      TmpInter[CellZ] = TmpXY;
		    }

		  this->RealHamiltonian[m1][n1][m2][n2] = new double [LengthZ];
		  this->ImaginaryHamiltonian[m1][n1][m2][n2] = new double [LengthZ];
		  TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m2][n2];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m2][n2];	
		  
		  for (p = 0; p < LengthZ; ++p)
		    {
		      TmpRealWaveFunctionOverlapZ = RealWaveFunctionOverlapZ[p];
		      TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[p];
		      
		      TmpRe = 0.0; TmpIm = 0.0;
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			{			  			  
			  TmpRe += (TmpInter[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ]);
			  TmpIm += (TmpInter[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ]);				    
			}
		      TmpRealHamiltonian[p] = TmpRe;
		      TmpImaginaryHamiltonian[p] = TmpIm;
		    }
		}		  	      
	    }
	}
    }  
  cout << "Hamitonian elements: " << endl;
  cout << this->RealHamiltonian[1][3][2][4][this->NbrStateZ - 1] << '\t' << this->ImaginaryHamiltonian[1][3][2][4][this->NbrStateZ - 1] << endl;
  cout << this->RealHamiltonian[1][3][2][4][2] << '\t' << this->ImaginaryHamiltonian[1][3][2][4][2] << endl;
  cout << this->RealHamiltonian[2][3][1][4][2] << '\t' << this->ImaginaryHamiltonian[2][3][1][4][2] << endl;
  cout << this->RealHamiltonian[2][4][1][3][2] << '\t' << this->ImaginaryHamiltonian[2][4][1][3][2] << endl;
  
  delete[] TmpInter;
}

// evaluate sinus wave function overlaps on a cell in a given direction
//
// size = system length in the choosen direction
// nbrStep = number of subdivision in the choosen direction
// nbrState = number of state in the choosen direction
// memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
// return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)

double*** XYReflexionSymmetricPeriodic3DHamiltonian::EvaluateSinusWaveFunctionOverlap(double size, int nbrStep, int nbrState)
{
  double*** TmpArray = new double** [nbrState];
  double StepInc = 1.0 / ((double) nbrStep);
  double Tmp;
  double Diff;
  double p = 0.0;
  for (int i = 0; i < nbrState; ++i)
    {
      TmpArray[i] = new double* [nbrState];
      for (int j = 0; j < i; ++j)
	{
	  TmpArray[i][j] = new double [nbrStep];
	  p = -double(nbrStep) / 2.0;
	  for (int k = 0; k < nbrStep; ++k)
	    {
	      Diff = (double) 2 * (i - j);
	      Tmp = M_PI * Diff * StepInc;
	      TmpArray[i][j][k] = M_1_PI * (sin (Tmp * ((double) (p + 1))) - sin (Tmp * ((double) (p)))) / Diff;
	      Diff = (double) 2 * (i + j + 2);
	      Tmp = M_PI * Diff * StepInc;
	      TmpArray[i][j][k] -= M_1_PI * (sin (Tmp * ((double) (p + 1))) - sin (Tmp * ((double) (p)))) / Diff;
	      p += 1.0;
	    }
	}
      TmpArray[i][i] = new double [nbrStep];
      p = -double(nbrStep) / 2.0;
      for (int k = 0; k < nbrStep; ++k)
	{
	  Tmp = M_PI * (double) (4 * i + 4) * StepInc;
	  TmpArray[i][i][k] = StepInc - M_1_PI * (sin (Tmp * ((double) (p + 1))) - sin (Tmp * ((double) (p)))) / ((double) (4 * i + 4));
	  p += 1.0;
	}
      for (int j = i + 1; j < nbrState; ++j)
	{
	  TmpArray[i][j] = new double [nbrStep];
	  p = -double(nbrStep) / 2.0;
	  for (int k = 0; k < nbrStep; ++k)
	    {
	      Diff = (double) 2 * (i - j);
	      Tmp = M_PI * Diff * StepInc;
	      TmpArray[i][j][k] = M_1_PI * (sin (Tmp * ((double) (p + 1))) - sin (Tmp * ((double) (p)))) / Diff;
	      Diff = (double) 2 * (i + j + 2);
	      Tmp = M_PI * Diff * StepInc;
	      TmpArray[i][j][k] -= M_1_PI * (sin (Tmp * ((double) (p + 1))) - sin (Tmp * ((double) (p)))) / Diff;
	      p += 1.0;
	    }
	}
    }
  return TmpArray;
}

// evaluate cosinus wave function overlaps on a cell in a given direction
//
// size = system length in the choosen direction
// nbrStep = number of subdivision in the choosen direction
// nbrState = number of state in the choosen direction
// memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
// return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)

double*** XYReflexionSymmetricPeriodic3DHamiltonian::EvaluateCosinusWaveFunctionOverlap(double size, int nbrStep, int nbrState)
{
  double*** TmpArray = new double** [nbrState];
  double StepInc = 1.0 / ((double) nbrStep);
  double Tmp;
  double Diff;
  double p = 0.0;
  for (int i = 0; i < nbrState; ++i)
    {
      TmpArray[i] = new double* [nbrState];
      for (int j = 0; j < i; ++j)
	{
	  TmpArray[i][j] = new double [nbrStep];
	  p = -double(nbrStep) / 2.0;
	  for (int k = 0; k < nbrStep; ++k)
	    {	      
	      Diff = (double) 2 * (i - j);
	      Tmp = M_PI * Diff * StepInc;
	      TmpArray[i][j][k] = M_1_PI * (sin (Tmp * ((double) (p + 1))) - sin (Tmp * ((double) (p)))) / Diff;
	      Diff = (double) 2 * (i + j);
	      Tmp = M_PI * Diff * StepInc;
	      TmpArray[i][j][k] += M_1_PI * (sin (Tmp * ((double) (p + 1))) - sin (Tmp * ((double) (p)))) / Diff;
	      p += 1.0;
	    }
	}
      TmpArray[i][i] = new double [nbrStep];
      p = -double(nbrStep) / 2.0;
      for (int k = 0; k < nbrStep; ++k)
	{
	  Tmp = M_PI * (double) (4 * i) * StepInc;	      
	  TmpArray[i][i][k] = StepInc + M_1_PI * (sin (Tmp * ((double) (p + 1))) - sin (Tmp * ((double) (p)))) / ((double) (4 * i + 4));
	  p += 1.0;
	}
      for (int j = i + 1; j < nbrState; ++j)
	{
	  TmpArray[i][j] = new double [nbrStep];
	  p = -double(nbrStep) / 2.0;
	  for (int k = 0; k < nbrStep; ++k)
	    {
	      Diff = (double) 2 * (i - j);
	      Tmp = M_PI * Diff * StepInc;
	      TmpArray[i][j][k] = M_1_PI * (sin (Tmp * ((double) (p + 1))) - sin (Tmp * ((double) (p)))) / Diff;
	      Diff = (double) 2 * (i + j);
	      Tmp = M_PI * Diff * StepInc;
	      TmpArray[i][j][k] += M_1_PI * (sin (Tmp * ((double) (p + 1))) - sin (Tmp * ((double) (p)))) / Diff;
	      p += 1.0;
	    }
	}
    }
  return TmpArray;
}

// evaluate the plane wave function overlap
//
// nbrStep = number of steps in the given direction
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool XYReflexionSymmetricPeriodic3DHamiltonian::EvaluatePlaneWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray)
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
