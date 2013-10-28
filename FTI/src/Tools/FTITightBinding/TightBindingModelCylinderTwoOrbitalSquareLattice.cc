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
    this->Inversion = ComplexMatrix(this->NbrBands, this->NbrBands, true);

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
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      double x = this->KxFactor * (kx + this->GammaX);
      double y = 0.0;
      int Index = kx;

      if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	{
	  Complex B1 = 2.0 * Complex(this->NNHoppingInterX * sin(this->FoldingFactor * x), - this->NNHoppingInterY * sin(this->FoldingFactor * y));
	  double d3 = this->MuS - 2.0 * this->NNHoppingIntra * (cos(x) + cos(y));
          
	  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, + d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, - d3);
	  
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
