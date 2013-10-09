////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of tight binding model for the kagome lattice          //
//                                                                            //
//                        last modification : 07/05/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelKagomeLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A1 sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelKagomeLattice::TightBindingModelKagomeLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double lambda1, double lambda2, double mus, 
							       double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->NNSpinOrbit = lambda1;
  this->NextNNSpinOrbit = lambda2;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 3;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->ProjectedMomenta = new double* [this->NbrStatePerBand];
  for (int i = 0; i < this->NbrStatePerBand; ++i)
    this->ProjectedMomenta[i] = new double [2];
  
  this->ComputeAllProjectedMomenta();
  

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


// constructor for a tilted lattice
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nx1 = first coordinate of the first spanning vector of the tilted lattice
// ny1 = second coordinate of the first spanning vector of the tilted lattice
// nx2 = first coordinate of the second spanning vector of the tilted lattice
// ny2 = second coordinate of the second spanning vector of the tilted lattice
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A1 sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelKagomeLattice::TightBindingModelKagomeLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double lambda1, double lambda2, double mus, 
							       double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->Nx1 = nx1;
  this->Ny1 = ny1;
  this->Nx2 = nx2;
  this->Ny2 = ny2;
  this->Offset = offset;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->NNSpinOrbit = lambda1;
  this->NextNNSpinOrbit = lambda2;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 3;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->ProjectedMomenta = new double* [this->NbrStatePerBand];
  for (int i = 0; i < this->NbrStatePerBand; ++i)
    this->ProjectedMomenta[i] = new double [2];
  
  this->ComputeAllProjectedMomenta();
  
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

TightBindingModelKagomeLattice::~TightBindingModelKagomeLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelKagomeLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	      KX = this->GetProjectedMomentum(kx, ky, 0);
	      KY = this->GetProjectedMomentum(kx, ky, 1);
	      Complex HAB (-2.0 * this->NNHopping, -2.0 * this->NNSpinOrbit);
	      HAB *= cos (KX);
	      Complex HAC(-2.0 * this->NNHopping, 2.0 * this->NNSpinOrbit);
	      HAC *= cos (KY);
	      Complex HBC(-2.0 * this->NNHopping, -2.0 * this->NNSpinOrbit);
	      HBC *= cos(KX - KY);
	      
	      Complex HAB2 (-2.0 * this->NextNNHopping, 2.0 * this->NextNNSpinOrbit);
	      HAB2 *= cos (KX - 2.0 * KY);
	      Complex HAC2 (-2.0 * this->NextNNHopping, -2.0 * this->NextNNSpinOrbit);
	      HAC2 *= cos (2.0 * KX - KY);
	      Complex HBC2 (-2.0 * this->NextNNHopping, 2.0  *  this->NextNNSpinOrbit);
	      HBC2 *= cos (KX + KY);
	      
	      HAB += HAB2;
	      HAC += HAC2;
	      HBC += HBC2;

	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, this->MuS);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HAB);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HAC);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HBC);
	      
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

