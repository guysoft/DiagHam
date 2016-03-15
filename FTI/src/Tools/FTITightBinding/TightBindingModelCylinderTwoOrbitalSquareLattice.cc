////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//  class of tight binding model for the two-orbital model on square lattice  //
//          with periodic boundary conditions only along the x direction      //
//                                                                            //
//                        last modification : 27/10/2013                      //
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
#include "Tools/FTITightBinding/TightBindingModelCylinderTwoOrbitalSquareLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "MathTools/BinomialCoefficients.h"
#include "Matrix/HermitianMatrix.h"

#include <iostream>


using std::cout;
using std::endl;


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = imag part of the inter-orbital hopping amplitude between nearest neighbors along the x direction
// t2 = the inter-orbital hopping amplitude between nearest neighbors along the y direction
// t3 = the intra-orbital hopping amplitude between nearest neighbors
// foldingFactor = folding factor for the momenta along sigma_x and sigma_y
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelCylinderTwoOrbitalSquareLattice::TightBindingModelCylinderTwoOrbitalSquareLattice(int nbrSiteX, int nbrSiteY, 
												   int t1, int t2, int t3, int foldingFactor, 
												   double mus, double gammaX,
												   AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
    this->NbrSiteX = nbrSiteX;
    this->NbrSiteY = nbrSiteY;
    this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
    this->NNHoppingInterX = t1;
    this->NNHoppingInterY = t2;
    this->NNHoppingIntra = t3;
    this->FoldingFactor = foldingFactor;
    this->MuS = mus;
    this->GammaX = gammaX;
    this->NbrBands = 2 * this->NbrSiteY;
    this->NbrStatePerBand = this->NbrSiteX;

    this->EmbeddingX = RealVector(this->NbrBands, true);

    this->Architecture = architecture;

    if (storeOneBodyMatrices == true)
        this->OneBodyBasis = new ComplexMatrix[this->NbrStatePerBand];
    else
        this->OneBodyBasis = 0;
    this->EnergyBandStructure = new double*[this->NbrBands];
    for (int i = 0; i < this->NbrBands; ++i)
        this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    this->ComputeBandStructure();
}

// destructor
//

TightBindingModelCylinderTwoOrbitalSquareLattice::~TightBindingModelCylinderTwoOrbitalSquareLattice()
{
}

// compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelCylinderTwoOrbitalSquareLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      double x = this->KxFactor * (kx + this->GammaX);
      int Index = kx;
      if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	{
	  TmpOneBodyHamiltonian.ClearMatrix();
	  double B1 = 2.0 * this->NNHoppingInterX * sin(this->FoldingFactor * x);
	  double d3 = this->MuS - 2.0 * this->NNHoppingIntra * cos(x);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, + d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, - d3);
	  for (int y = 1; y < this->NbrSiteY; ++y)
	    {
	      TmpOneBodyHamiltonian.SetMatrixElement(2 * y, 2 * y, + d3);
	      TmpOneBodyHamiltonian.SetMatrixElement(2 * y, 2 * y + 1, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(2 * y + 1, 2 * y + 1, - d3);
	      TmpOneBodyHamiltonian.SetMatrixElement(2 * (y - 1), 2 * y, - this->NNHoppingIntra);
	      TmpOneBodyHamiltonian.SetMatrixElement(2 * (y - 1) + 1, 2 * y + 1,  this->NNHoppingIntra);
	      TmpOneBodyHamiltonian.SetMatrixElement(2 * (y - 1), 2 * y + 1, this->NNHoppingInterY);
	      TmpOneBodyHamiltonian.SetMatrixElement(2 * (y - 1) + 1, 2 * y, -this->NNHoppingInterY);
	    }
//     turn on periodic boundary conditions
// 	  if (this->NbrSiteY > 2)
// 	    {
// 	      TmpOneBodyHamiltonian.SetMatrixElement(2 * (this->NbrSiteY - 1), 0, - this->NNHoppingIntra);
// 	      TmpOneBodyHamiltonian.SetMatrixElement(2 * (this->NbrSiteY - 1) + 1, 1,  this->NNHoppingIntra);
// 	      TmpOneBodyHamiltonian.SetMatrixElement(2 * (this->NbrSiteY - 1), 1, this->NNHoppingInterY);
// 	      TmpOneBodyHamiltonian.SetMatrixElement(2 * (this->NbrSiteY - 1) + 1, 0, -this->NNHoppingInterY);
// 	    }
// 	  else
// 	    {
// 	      TmpOneBodyHamiltonian.AddToMatrixElement(2 * (this->NbrSiteY - 1), 0, -this->NNHoppingIntra);
// 	      TmpOneBodyHamiltonian.AddToMatrixElement(2 * (this->NbrSiteY - 1) + 1, 1, this->NNHoppingIntra);
// 	      TmpOneBodyHamiltonian.AddToMatrixElement(2 * (this->NbrSiteY - 1), 1, this->NNHoppingInterY);
// 	      TmpOneBodyHamiltonian.AddToMatrixElement(2 * (this->NbrSiteY - 1) + 1, 0, this->NNHoppingInterY);
// 	    }
	  if (this->OneBodyBasis != 0)
	    {
	      ComplexMatrix TmpMatrix(this->NbrBands, this->NbrBands, true);
	      TmpMatrix.SetToIdentity();
	      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	      TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	      TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
	      this->OneBodyBasis[Index] = TmpMatrix;
	      for (int i = 0; i < this->NbrBands; ++i)
		this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
	    }
	  else
	    {
	      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	      TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
	      TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
	      for (int i = 0; i < this->NbrBands; ++i)
		this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
	    }
	}
    }
}


// compute the many-body real space entanglement spectrum of a full band
// 
// outputFile = name of the output file where the spectrum has to be stored
// minEnergy = lowest energy of the full band
// maxEnergy = highest energy of the full band
// nbrSiteYA = number of site to keep for the A part along the y direction    

void TightBindingModelCylinderTwoOrbitalSquareLattice::ComputeManyBodyRealSpaceEntanglementSpectrum(char* outputFile, double minEnergy, double maxEnergy, int nbrSiteYA)
{
  if (this->HaveOneBodyBasis() == false)
    {
      cout << "error, the tight binding model does not provide the one body basis" << endl;
    }
  int* MinIndex = new int[this->NbrSiteX];
  int* MaxIndex = new int[this->NbrSiteX];
  long TotalNbrStates = 0l;
  double** Weights = new double* [this->NbrSiteX];
  int* NbrKeptWeights = new int [this->NbrSiteX];
  int** KeptWeights = new int* [this->NbrSiteX];
  cout << "nbrSiteYA " << nbrSiteYA << endl;
  int TwiceNbrSiteYA = 2 * nbrSiteYA;
  for (int Index = 0; Index < this->NbrSiteX; ++Index)
    {
      MinIndex[Index] = this->NbrBands;
      MaxIndex[Index] = -1;
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  if ((this->EnergyBandStructure[i][Index] >= minEnergy) && (this->EnergyBandStructure[i][Index] <= maxEnergy))
	    {
	      ++TotalNbrStates;
	      if (MinIndex[Index] > i)
		{
		  MinIndex[Index] = i;
		}
	      if (MaxIndex[Index] < i)
		{
		  MaxIndex[Index] = i;
		}
	    }
	}
      NbrKeptWeights[Index] = 0;
      if (MaxIndex[Index] >= MinIndex[Index])
	{
	  Weights[Index] = new double [MaxIndex[Index] - MinIndex[Index] + 1];
	  KeptWeights[Index] = new int [MaxIndex[Index] - MinIndex[Index] + 1];
	  for (int i = MinIndex[Index]; i <= MaxIndex[Index]; ++i)
	    {
	      double Tmp = 0.0;
	      for (int j = 0; j < TwiceNbrSiteYA; ++j)
		{
		  Tmp += SqrNorm(this->OneBodyBasis[Index][i][j]);
		}
// 	      Complex Tmp = 0.0;
// 	      for (int j = 0; j < TwiceNbrSiteYA; ++j)
// 		{
// 		  Tmp += Conj(this->OneBodyBasis[Index][i][j]);
// 		}
// 	      Complex Tmp2 = 0.0;
// 	      for (int j = TwiceNbrSiteYA; j < this->NbrBands; ++j)
// 		{
// 		  Tmp2 += Conj(this->OneBodyBasis[Index][i][j]);
// 		}
//  	      Weights[Index][i - MinIndex[Index]] = Norm(Tmp) / sqrt(SqrNorm(Tmp) + SqrNorm(Tmp2));
	      Weights[Index][i - MinIndex[Index]] = sqrt(Tmp);
// 	      if (((1.0 - Weights[Index][i - MinIndex[Index]]) != 0.0) && ((1.0 - Weights[Index][i - MinIndex[Index]]) != 1.0))
 		{
		  KeptWeights[Index][NbrKeptWeights[Index]] = i - MinIndex[Index];
		  NbrKeptWeights[Index]++;
		}
	    }
	}
      else
	{
	  Weights[Index] = 0;
	}
    }
  cout << "nbr of states in the Slater determinant = " <<   TotalNbrStates << endl;
  cout << "nbr of states with a non zero weight : " << endl; 
  ofstream File;
  File.open(outputFile);
  File << "# kx  weight_A " << endl;
  for (int Index = 0; Index < this->NbrSiteX; ++Index)
    {
      cout << "total nbr of states in the kx=" << Index << " sector : " << NbrKeptWeights[Index] << endl;
      for (int i = 0; i < NbrKeptWeights[Index]; ++i)
	{
	  File << Index << " " << Weights[Index][KeptWeights[Index][i]] << endl;
	}
      delete[] KeptWeights[Index];
      delete[] Weights[Index];
    }  
  File.close();
  delete[] NbrKeptWeights;
  delete[] KeptWeights;
  delete[] Weights;
}

// compute the one-body real space entanglement spectrum of a full band
// 
// outputFile = name of the output file where the spectrum has to be stored
// minEnergy = lowest energy of the full band
// maxEnergy = highest energy of the full band
// nbrSiteYA = number of site to keep for the A part along the y direction    

void TightBindingModelCylinderTwoOrbitalSquareLattice::ComputeOneBodyRealSpaceEntanglementSpectrum(char* outputFile, double minEnergy, double maxEnergy, int nbrSiteYA)
{
  if (this->HaveOneBodyBasis() == false)
    {
      cout << "error, the tight binding model does not provide the one body basis" << endl;
    }
  int* MinIndex = new int[this->NbrSiteX];
  int* MaxIndex = new int[this->NbrSiteX];
  long TotalNbrStates = 0l;
  double** Weights = new double* [this->NbrSiteX];
  int* NbrKeptWeights = new int [this->NbrSiteX];
  int** KeptWeights = new int* [this->NbrSiteX];
  cout << "nbrSiteYA " << nbrSiteYA << endl;
  int NbrOrbitalPerUnitCell =  this->NbrBands / this->NbrSiteY;
  int PropagatorDimension = NbrOrbitalPerUnitCell * nbrSiteYA;
  ofstream File;
  File.open(outputFile);
  File << "# kx  weight_A " << endl;
  for (int TmpKx = 0; TmpKx < this->NbrSiteX; ++TmpKx)
    {
      int MinIndex = this->NbrBands;
      int MaxIndex = -1;
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  if ((this->EnergyBandStructure[i][TmpKx] >= minEnergy) && (this->EnergyBandStructure[i][TmpKx] <= maxEnergy))
	    {
	      ++TotalNbrStates;
	      if (MinIndex > i)
		{
		  MinIndex = i;
		}
	      if (MaxIndex < i)
		{
		  MaxIndex = i;
		}
	    }
	}
      if ((MinIndex != this->NbrBands) && (MaxIndex >= 0))
	{
	  cout << "nbr of states at kx=" << TmpKx << " : " <<  (MaxIndex - MinIndex + 1)  << endl;
	  HermitianMatrix Propagator (PropagatorDimension, true);
	  for (; MinIndex <= MaxIndex; ++MinIndex)
	    {
	      for (int i = 0; i < PropagatorDimension; ++i)
		{
		  Complex Tmp = Conj(this->OneBodyBasis[TmpKx][MinIndex][i]);
		  Propagator.AddToMatrixElement(i, i, Tmp * this->OneBodyBasis[TmpKx][MinIndex][i]);		  
		  for (int j = i + 1; j < PropagatorDimension; ++j)
		    {
		      Propagator.AddToMatrixElement(i, j, Tmp * this->OneBodyBasis[TmpKx][MinIndex][j]);		  
		    }
		}
	    }
	  RealDiagonalMatrix TmpDiag (Propagator.GetNbrRow(), true);
#ifdef __LAPACK__
	  Propagator.LapackDiagonalize(TmpDiag);
#else
	  Propagator.Diagonalize(TmpDiag);
#endif		  
	  for (int i = 0; i < PropagatorDimension; ++i)	 
	    File << TmpKx << " " << TmpDiag[i] << endl;
	}
    }
  File.close();
}

// evaluate the mixed two point correlation function in a given region, assuming translation invariance along the x direction
//
// maxY = length along the y direction of the region 
// kx = momentum along the x direction
// occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
// nbrOccupiedMomenta = number of occupied momenta
// bandIndex = index of the band to consider
// return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)

HermitianMatrix TightBindingModelCylinderTwoOrbitalSquareLattice::EvaluateFullMixedTwoPointCorrelationFunctionWithKy(int maxY, int kx, int* occupiedMomenta, int nbrOccupiedMomenta, int bandIndex)
{
  int TmpNbrOrbitalPerUnitCell = this->NbrBands / this->NbrSiteY;
  int TotalNbrSites = maxY * TmpNbrOrbitalPerUnitCell;
  int TmpMomentumX;
  int TmpMomentumY;
  HermitianMatrix EntanglementHamiltonian(TotalNbrSites, true);
  int TmpPos = -1;
  for (int i = 0; (i < nbrOccupiedMomenta) && (TmpPos == -1); ++i)
    {
      if (occupiedMomenta[i] == kx)
	{
	  TmpPos = i;
	}
    }
  if (TmpPos == -1)
    {
      return EntanglementHamiltonian;
    }


  return EntanglementHamiltonian;
}
  
