////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                     Class author Cecile Repellin                           //
//                                                                            //
//               class of tight binding model for the kagome lattice          //
//                    with tilted boundary conditions                         //
//                      last modification : 14/05/2013                        //
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
#include "Tools/FTITightBinding/TightBindingModelKagomeLatticeTilted.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
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

TightBindingModelKagomeLatticeTilted::TightBindingModelKagomeLatticeTilted(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double lambda1, double lambda2, double mus, 
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

TightBindingModelKagomeLatticeTilted::~TightBindingModelKagomeLatticeTilted()
{
  delete[] this->ProjectedMomenta;
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelKagomeLatticeTilted::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	      Complex HAB (-1.0*this->NNHopping, -1.0*this->NNSpinOrbit);
	      HAB *= 1 + Phase(KX);
	      Complex HAC(-1.0*this->NNHopping, this->NNSpinOrbit);
	      HAC *= 1 + Phase(KY);
	      Complex HBC(-1.0*this->NNHopping, -1.0*this->NNSpinOrbit);
	      HBC *= 1 + Phase(KY - KX);
	      
// 	      Complex HAB2 (-2.0 * this->NextNNHopping, 2.0 * this->NextNNSpinOrbit);
// 	      HAB2 *= cos (KX - 2.0 * KY);
// 	      Complex HAC2 (-2.0 * this->NextNNHopping, -2.0 * this->NextNNSpinOrbit);
// 	      HAC2 *= cos (2.0 * KX - KY);
// 	      Complex HBC2 (-2.0 * this->NextNNHopping, 2.0  *  this->NextNNSpinOrbit);
// 	      HBC2 *= cos (KX + KY);
// 	      
// 	      HAB += HAB2;
// 	      HAC += HAC2;
// 	      HBC += HBC2;

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

//computes all the values of the projected momentum and stores them in a double array
//
void TightBindingModelKagomeLatticeTilted::ComputeAllProjectedMomenta()
{
 double projectedMomentum1;
 double projectedMomentum2;
 for (int kx = 0; kx < this->NbrSiteX; ++kx)
 {
   for (int ky = 0; ky < this->NbrSiteY; ++ky)
   {
     int kx_trans = kx + this->Offset*ky;
     int ky_trans = ky;
     projectedMomentum1 = 2.0 * M_PI * ((double) kx_trans * (double) this->Ny2 - (double) ky_trans * (double) this->Ny1) / ((double) (this->NbrSiteX * this->NbrSiteY));
     projectedMomentum2 = 2.0 * M_PI * ((double) kx_trans * (double) (-this->Nx2) + (double) ky_trans * (double)this->Nx1) / ((double) (this->NbrSiteX * this->NbrSiteY));
     this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky)][0] = projectedMomentum1;
     this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky)][1] = projectedMomentum2;
   }
 }
}