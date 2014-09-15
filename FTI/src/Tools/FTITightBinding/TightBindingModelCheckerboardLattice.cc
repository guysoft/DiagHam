////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of tight binding model for the checkerboard lattice        //
//                                                                            //
//                        last modification : 08/05/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// t2p = hoping amplitude between second next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelCheckerboardLattice::TightBindingModelCheckerboardLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double t2p, double mus, 
									   double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->SecondNextNNHoping = t2p;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 2;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;

  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
    }
  else
    {
      this->OneBodyBasis = 0;
    }
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    }
  this->ComputeBandStructure();
}

// destructor
//

TightBindingModelCheckerboardLattice::~TightBindingModelCheckerboardLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelCheckerboardLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double KX;
  double KY;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int Index = this->GetLinearizedMomentumIndex(kx, ky);
	  if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	    {
	      Complex B1 = 4.0 * this->NNHoping * Complex (cos (1.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (1.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)) * cos(M_PI * 0.25), 
							   sin (1.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * sin (1.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)) * sin(M_PI * 0.25));
	      double d1 = 4.0 * this->SecondNextNNHoping * cos (2.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (2.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY));
	      double d3 =  this->MuS + (2.0 * this->NextNNHoping * (cos (2.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX))
								    - cos (2.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY))));

	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	      
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
}

// get the tight binding hamiltonian in real space 
// 
// return value = tight binding hamiltonian

HermitianMatrix TightBindingModelCheckerboardLattice::GetRealSpaceTightBindingHamiltonian()
{
  int* NbrConnectedOrbitals = new int [this->NbrBands];
  int** OrbitalIndices = new int* [this->NbrBands];
  int** SpatialIndices = new int* [this->NbrBands];
  Complex** HoppingAmplitudes = new Complex* [this->NbrBands];
  if (this->SecondNextNNHoping != 0.0)
    {
      NbrConnectedOrbitals[0] = 13; 
      NbrConnectedOrbitals[1] = 13;      
    } 
  else
    {
      NbrConnectedOrbitals[0] = 9; 
      NbrConnectedOrbitals[1] = 9;
    }
  for (int i = 0; i < this->NbrBands; ++i)
    {
      OrbitalIndices[i] = new int[NbrConnectedOrbitals[i]];
      SpatialIndices[i] = new int[2 * NbrConnectedOrbitals[i]];
      HoppingAmplitudes[i] = new Complex[NbrConnectedOrbitals[i]];
    }

  int TmpIndex = 0;
  OrbitalIndices[0][TmpIndex] = 0;
  SpatialIndices[0][TmpIndex * 2] = 0;
  SpatialIndices[0][(TmpIndex * 2) +1] = 0;
  HoppingAmplitudes[0][TmpIndex] = this->MuS;
  OrbitalIndices[1][TmpIndex] = 1;
  SpatialIndices[1][TmpIndex * 2] = 0;
  SpatialIndices[1][(TmpIndex * 2) +1] = 0;
  HoppingAmplitudes[1][TmpIndex] = -this->MuS;
  ++TmpIndex;

  OrbitalIndices[0][TmpIndex] = 1;
  SpatialIndices[0][TmpIndex * 2] = 0;
  SpatialIndices[0][(TmpIndex * 2) +1] = 0;
  HoppingAmplitudes[0][TmpIndex] = this->NNHoping * Phase (M_PI * 0.25);
  OrbitalIndices[1][TmpIndex] = 0;
  SpatialIndices[1][TmpIndex * 2] = 0;
  SpatialIndices[1][(TmpIndex * 2) +1] = 0;
  HoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (-M_PI * 0.25);
  ++TmpIndex;
  OrbitalIndices[0][TmpIndex] = 1;
  SpatialIndices[0][TmpIndex * 2] = -1;
  SpatialIndices[0][(TmpIndex * 2) +1] = -1;
  HoppingAmplitudes[0][TmpIndex] = this->NNHoping * Phase (M_PI * 0.25);
  OrbitalIndices[1][TmpIndex] = 0;
  SpatialIndices[1][TmpIndex * 2] = 1;
  SpatialIndices[1][(TmpIndex * 2) +1] = 1;
  HoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (-M_PI * 0.25);
  ++TmpIndex;
  OrbitalIndices[0][TmpIndex] = 1;
  SpatialIndices[0][TmpIndex * 2] = -1;
  SpatialIndices[0][(TmpIndex * 2) +1] = 0;
  HoppingAmplitudes[0][TmpIndex] = this->NNHoping * Phase (-M_PI * 0.25);
  OrbitalIndices[1][TmpIndex] = 0;
  SpatialIndices[1][TmpIndex * 2] = 0;
  SpatialIndices[1][(TmpIndex * 2) +1] = 1;
  HoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (M_PI * 0.25);
  ++TmpIndex;
  OrbitalIndices[0][TmpIndex] = 1;
  SpatialIndices[0][TmpIndex * 2] = 0;
  SpatialIndices[0][(TmpIndex * 2) +1] = -1;
  HoppingAmplitudes[0][TmpIndex] = this->NNHoping * Phase (-M_PI * 0.25);
  OrbitalIndices[1][TmpIndex] = 0;
  SpatialIndices[1][TmpIndex * 2] = 1;
  SpatialIndices[1][(TmpIndex * 2) +1] = 0;
  HoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (M_PI * 0.25);
  ++TmpIndex;

  OrbitalIndices[0][TmpIndex] = 0;
  SpatialIndices[0][TmpIndex * 2] = 1;
  SpatialIndices[0][(TmpIndex * 2) +1] = 0;
  HoppingAmplitudes[0][TmpIndex] = this->NextNNHoping;
  OrbitalIndices[1][TmpIndex] = 1;
  SpatialIndices[1][TmpIndex * 2] = 1;
  SpatialIndices[1][(TmpIndex * 2) +1] = 0;
  HoppingAmplitudes[1][TmpIndex] = -this->NextNNHoping;
  ++TmpIndex;
  OrbitalIndices[0][TmpIndex] = 0;
  SpatialIndices[0][TmpIndex * 2] = -1;
  SpatialIndices[0][(TmpIndex * 2) +1] = 0;
  HoppingAmplitudes[0][TmpIndex] = this->NextNNHoping;
  OrbitalIndices[1][TmpIndex] = 1;
  SpatialIndices[1][TmpIndex * 2] = -1;
  SpatialIndices[1][(TmpIndex * 2) +1] = 0;
  HoppingAmplitudes[1][TmpIndex] = -this->NextNNHoping;
  ++TmpIndex;
  OrbitalIndices[0][TmpIndex] = 0;
  SpatialIndices[0][TmpIndex * 2] = 0;
  SpatialIndices[0][(TmpIndex * 2) +1] = 1;
  HoppingAmplitudes[0][TmpIndex] = -this->NextNNHoping;
  OrbitalIndices[1][TmpIndex] = 1;
  SpatialIndices[1][TmpIndex * 2] = 0;
  SpatialIndices[1][(TmpIndex * 2) +1] = 1;
  HoppingAmplitudes[1][TmpIndex] = this->NextNNHoping;
  ++TmpIndex;
  OrbitalIndices[0][TmpIndex] = 0;
  SpatialIndices[0][TmpIndex * 2] = 0;
  SpatialIndices[0][(TmpIndex * 2) +1] = -1;
  HoppingAmplitudes[0][TmpIndex] = -this->NextNNHoping;
  OrbitalIndices[1][TmpIndex] = 1;
  SpatialIndices[1][TmpIndex * 2] = 0;
  SpatialIndices[1][(TmpIndex * 2) +1] = -1;
  HoppingAmplitudes[1][TmpIndex] = this->NextNNHoping;
  ++TmpIndex;

  if (this->SecondNextNNHoping != 0.0)
    {
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  OrbitalIndices[i][TmpIndex] = i;
	  SpatialIndices[i][TmpIndex * 2] = 1;
	  SpatialIndices[i][(TmpIndex * 2) + 1] = 1;
	  HoppingAmplitudes[i][TmpIndex] = this->SecondNextNNHoping;
	}
      ++TmpIndex;
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  OrbitalIndices[i][TmpIndex] = i;
	  SpatialIndices[i][TmpIndex * 2] = -1;
	  SpatialIndices[i][(TmpIndex * 2) +1] = 1;
	  HoppingAmplitudes[i][TmpIndex] = this->SecondNextNNHoping;
	}
      ++TmpIndex;
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  OrbitalIndices[i][TmpIndex] = i;
	  SpatialIndices[i][TmpIndex * 2] = 1;
	  SpatialIndices[i][(TmpIndex * 2) +1] = -1;
	  HoppingAmplitudes[i][TmpIndex] = this->SecondNextNNHoping;
	}
      ++TmpIndex;
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  OrbitalIndices[i][TmpIndex] = i;
	  SpatialIndices[i][TmpIndex * 2] = -1;
	  SpatialIndices[i][(TmpIndex * 2) +1] = -1;
	  HoppingAmplitudes[i][TmpIndex] = this->SecondNextNNHoping;
	}
      ++TmpIndex;
    }


  HermitianMatrix TmpMatrix = this->BuildTightBindingHamiltonianRealSpace(NbrConnectedOrbitals, OrbitalIndices, SpatialIndices, HoppingAmplitudes);
  for (int i = 0; i < this->NbrBands; ++i)
    {
      delete[] HoppingAmplitudes[i];
      delete[] SpatialIndices[i];
      delete[] OrbitalIndices[i];
    }
  delete[] HoppingAmplitudes;
  delete[] SpatialIndices;
  delete[] OrbitalIndices;
  delete[] NbrConnectedOrbitals;
  return TmpMatrix;
}

