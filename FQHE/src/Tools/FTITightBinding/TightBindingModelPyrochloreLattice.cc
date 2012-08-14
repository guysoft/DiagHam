////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of tight binding model for the 3D pyrochlore lattice         //
//                                                                            //
//                        last modification : 10/08/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelPyrochloreLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// lambdaNN = spin orbit coupling to neareast neighbor sites
// lambdaNNN = spin orbit coupling to next neareast neighbor sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// gammaZ = boundary condition twisting angle along y
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelPyrochloreLattice::TightBindingModelPyrochloreLattice(int nbrSiteX, int nbrSiteY, int nbrSiteZ,
								       double lambdaNN, double lambdaNNN,
								       double gammaX, double gammaY, double gammaZ, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->KzFactor = 2.0 * M_PI / ((double) this->NbrSiteZ);
  this->NNSpinOrbit = lambdaNN;
  this->NextNNSpinOrbit = lambdaNNN;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->GammaZ = gammaZ;
  this->NbrBands = 8;
  this->NbrStatePerBand = 2 * this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;

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

TightBindingModelPyrochloreLattice::~TightBindingModelPyrochloreLattice()
{
}

// compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelPyrochloreLattice::ComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double KX;
  double KY;
  double KZ;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  for (int kz = 0; kz < this->NbrSiteZ; ++kz)
	    {
	      int Index = this->GetLinearizedMomentumIndex(kx, ky, kz);
	      if ((Index >= minStateIndex) && (Index < MaxStateIndex))
		{
		  KX = this->KxFactor * (((double) kx) + this->GammaX);
		  KY = this->KyFactor * (((double) ky) + this->GammaY);
		  KZ = this->KzFactor * (((double) kz) + this->GammaZ);
		  
		  Complex HHopping12 = 1.0 + Phase(KX-KZ);
		  Complex HHopping13 = 1.0 + Phase(KX);
		  Complex HHopping14 = 1.0 + Phase(KX-KY);
		  Complex HHopping23 = 1.0 + Phase(KZ);
		  Complex HHopping24 = 1.0 + Phase(KZ-KY);
		  Complex HHopping34 = 1.0 + Phase(-KY);

		  
		  Complex HSONN12 = this->NNSpinOrbit * I() * (1.0+Phase(KX-KZ));
		  Complex HSONN13 = this->NNSpinOrbit * I() * (1.0+Phase(KX));
		  Complex HSONN14 = this->NNSpinOrbit * I() * (1.0+Phase(KX-KY));
		  Complex HSONN23 = this->NNSpinOrbit * I() * (1.0+Phase(KZ));
		  Complex HSONN24 = this->NNSpinOrbit * I() * (Phase(-KY+KZ)+1.0);
		  Complex HSONN34 = this->NNSpinOrbit * I() * (Phase(-KY)+1.0);
		  
		  Complex HSONN56 = -this->NNSpinOrbit * I() * (1.0+Phase(KX-KZ));
		  Complex HSONN57 = -this->NNSpinOrbit * I() * (1.0+Phase(KX));
		  Complex HSONN58 = -this->NNSpinOrbit * I() * (1.0+Phase(KX-KY));
		  Complex HSONN67 = -this->NNSpinOrbit * I() * (1.0+Phase(KZ));
		  Complex HSONN68 = -this->NNSpinOrbit * I() * (Phase(-KY+KZ)+1.0);
		  Complex HSONN78 = -this->NNSpinOrbit * I() * (Phase(-KY)+1.0);


		  Complex HSONNN12 = this->NextNNSpinOrbit * I() * (-0.5 * Phase(-KZ) + 0.5 * Phase(KY - KZ) - 0.5 * Phase(KX - KY) + 0.5 * Phase(KX));
		  Complex HSONNN13 = this->NextNNSpinOrbit * I() * (Phase(KY) - Phase(KX - KY));
		  Complex HSONNN14 = this->NextNNSpinOrbit * I() * (- Phase(-KY) + Phase(KX));
		  Complex HSONNN23 = this->NextNNSpinOrbit * I() * (Phase(-KY + KZ) - Phase(KY));
		  Complex HSONNN24 = this->NextNNSpinOrbit * I() * (Phase(-KY) - Phase(KZ));
		  Complex HSONNN34 = this->NextNNSpinOrbit * I() * (-0.5 * Phase(-KX) - 0.5 * Phase(-KY + KZ) + 0.5 * Phase(-KZ) + 0.5 * Phase(KX - KY));
		  
		  Complex HSONNN56 = -this->NextNNSpinOrbit * I() * (-0.5 * Phase(-KZ) + 0.5 * Phase(KY - KZ) - 0.5 * Phase(KX - KY) + 0.5 * Phase(KX));
		  Complex HSONNN57 = -this->NextNNSpinOrbit * I() * (Phase(KY) - Phase(KX - KY));
		  Complex HSONNN58 = -this->NextNNSpinOrbit * I() * (- Phase(-KY) + Phase(KX));
		  Complex HSONNN67 = -this->NextNNSpinOrbit * I() * (Phase(-KY + KZ) - Phase(KY));
		  Complex HSONNN68 = -this->NextNNSpinOrbit * I() * (Phase(-KY) - Phase(KZ));
		  Complex HSONNN78 = -this->NextNNSpinOrbit * I() * (-0.5 * Phase(-KX) - 0.5 * Phase(-KY + KZ) + 0.5 * Phase(-KZ) + 0.5 * Phase(KX - KY));
		  
		  Complex HSONNN16 = this->NextNNSpinOrbit * ((I() * (Phase(-KZ) - Phase(KX)))   - (Phase(KY - KZ) - Phase(KX - KY)));
		  Complex HSONNN17 = this->NextNNSpinOrbit * ((I() * (- Phase(KZ) + Phase(KX - KZ)))     - (-0.5 * Phase(KZ) + 0.5 * Phase(KY) - 0.5 * Phase(KX - KY) + 0.5 * Phase(KX - KZ)));
		  Complex HSONNN18 = this->NextNNSpinOrbit * ((I() * (0.5 * Phase(-KY) - 0.5 * Phase(-KY + KZ) + 0.5 * Phase(KX - KZ) - 0.5 * Phase(KX))) - (- Phase(-KY + KZ) + Phase(KX - KZ)));
		  Complex HSONNN25 = this->NextNNSpinOrbit * ((I() * (- Phase(-KX + KY) + Phase(-KY + KZ)))      - (- Phase(-KX) + Phase(KZ)));
		  Complex HSONNN27 = this->NextNNSpinOrbit * ((I() * (-0.5 * Phase(-KX + KZ) + 0.5 * Phase(-KY + KZ) - 0.5 * Phase(KY) + 0.5 * Phase(KX)))        - (- Phase(-KX + KZ) + Phase(KX)));
		  Complex HSONNN28 = this->NextNNSpinOrbit * ((I() * (- Phase(-KX + KZ) + Phase(KX - KY)))       - (-0.5 * Phase(-KX + KZ) - 0.5 * Phase(-KY) + 0.5 * Phase(KZ) + 0.5 * Phase(KX - KY)));
		  Complex HSONNN35 = this->NextNNSpinOrbit * ((I() * (Phase(-KX + KY) - Phase(-KY)))     - (0.5 * Phase(-KX + KZ) - 0.5 * Phase(-KX + KY) + 0.5 * Phase(-KY) - 0.5 * Phase(-KZ)));
		  Complex HSONNN36 = this->NextNNSpinOrbit * ((I() * (0.5 * Phase(-KX) - 0.5 * Phase(-KY) + 0.5 * Phase(KY - KZ) - 0.5 * Phase(KX - KZ))) - (Phase(-KY) - Phase(KY - KZ)));
		  Complex HSONNN38 = this->NextNNSpinOrbit * ((I() * (Phase(-KX) - Phase(KX - KY)))      - (Phase(-KY + KZ) - Phase(-KZ)));
		  Complex HSONNN45 = this->NextNNSpinOrbit * ((I() * (-0.5 * Phase(-KX) + 0.5 * Phase(-KX + KZ) - 0.5 * Phase(KY - KZ) + 0.5 * Phase(KY)))        - (Phase(-KX) - Phase(KY)));
		  Complex HSONNN46 = this->NextNNSpinOrbit * ((I() * (- Phase(-KZ) + Phase(KY))) - (0.5 * Phase(-KX + KY) + 0.5 * Phase(-KZ) - 0.5 * Phase(KY) - 0.5 * Phase(KX - KZ)));
		  Complex HSONNN47 = this->NextNNSpinOrbit * ((I() * (Phase(KZ) - Phase(KY - KZ)))       - (Phase(-KX + KY) - Phase(KX)));

		  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
		  
		  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HHopping12 + HSONNN12);
		  TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HHopping13 + HSONNN13);
		  TmpOneBodyHamiltonian.SetMatrixElement(0, 3, HHopping14 + HSONNN14);
		  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HHopping23 + HSONNN23);
		  TmpOneBodyHamiltonian.SetMatrixElement(1, 3, HHopping24 + HSONNN24);
		  TmpOneBodyHamiltonian.SetMatrixElement(2, 3, HHopping34 + HSONNN34);
		  
		  TmpOneBodyHamiltonian.SetMatrixElement(4, 5, HHopping12 + HSONNN56);
		  TmpOneBodyHamiltonian.SetMatrixElement(4, 6, HHopping13 + HSONNN57);
		  TmpOneBodyHamiltonian.SetMatrixElement(4, 7, HHopping14 + HSONNN58);
		  TmpOneBodyHamiltonian.SetMatrixElement(5, 6, HHopping23 + HSONNN67);
		  TmpOneBodyHamiltonian.SetMatrixElement(5, 7, HHopping24 + HSONNN68);
		  TmpOneBodyHamiltonian.SetMatrixElement(6, 7, HHopping34 + HSONNN78);

		  TmpOneBodyHamiltonian.SetMatrixElement(0, 5, HSONNN16);
		  TmpOneBodyHamiltonian.SetMatrixElement(0, 6, HSONNN17);
		  TmpOneBodyHamiltonian.SetMatrixElement(0, 7, HSONNN18);
		  TmpOneBodyHamiltonian.SetMatrixElement(1, 4, HSONNN25);
		  TmpOneBodyHamiltonian.SetMatrixElement(1, 6, HSONNN27);
		  TmpOneBodyHamiltonian.SetMatrixElement(1, 7, HSONNN28);
		  TmpOneBodyHamiltonian.SetMatrixElement(2, 4, HSONNN35);
		  TmpOneBodyHamiltonian.SetMatrixElement(2, 5, HSONNN36);
		  TmpOneBodyHamiltonian.SetMatrixElement(2, 7, HSONNN38);
		  TmpOneBodyHamiltonian.SetMatrixElement(3, 4, HSONNN45);
		  TmpOneBodyHamiltonian.SetMatrixElement(3, 5, HSONNN46);
		  TmpOneBodyHamiltonian.SetMatrixElement(3, 6, HSONNN47);

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
}
