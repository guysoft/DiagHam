////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                         Class author : Cecile Repellin                     //
//                                                                            //
//           class of tight binding model for trilayer grahene-hBN            //
//                with a relative twist angle and Moire coupling              //
//                      with valley and spin conservation                     //
//               (must be TR conjugated to obtain other valley index)         //
//                   last modification : 04/06/2018                           //
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
#include "Tools/FTITightBinding/TightBindingModelTrilayerGrapheneTwistBN.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

using std::cout;
using std::endl;


// basic constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mixingTerm12 = mixing term coupling the two copies of the kagome lattice (sites 1 and 2)
// mixingTerm13 = mixing term coupling the two copies of the kagome lattice (sites 1 and 3)
// mixingTerm23 = mixing term coupling the two copies of the kagome lattice (sites 2 and 3)
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelTrilayerGrapheneTwistBN::TightBindingModelTrilayerGrapheneTwistBN(int nbrSiteX, int nbrSiteY, int nbrPointsX, int nbrPointsY, double uVoltage, double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
   
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  
  
  this->NNHoppingIntra = 2676;
  this->GammaOne = 340;
  this->TrigonalWarping = 0.051*this->NNHoppingIntra;
  this->UVoltage = uVoltage;
    
  this->NbrPointSuperlatticeX = nbrPointsX;
  this->NbrPointSuperlatticeY = nbrPointsY;
  
    
  this->KxFactor = 1.0 / ((double) this->NbrSiteX);
  this->KyFactor = 1.0 / ((double) this->NbrSiteY);
  
  
  this->GammaX = gammaX;
  this->GammaY = gammaY; 
  
  this->LatticeConstant = 58.8;
  this->Gx1 = 0.0;
  this->Gy1 = -4.0 * M_PI / (sqrt(3.0) * this->LatticeConstant);
  this->Gx2 = 2.0*M_PI / this->LatticeConstant;
  this->Gy2 = this->Gy1 / 2.0;
  
  this->NbrBands = 6*nbrPointsX*nbrPointsY;
  this->BandIndex = this->NbrBands/2 - 1;
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
  
  
//   for (int kx = 0; kx < this->NbrSiteX; ++kx)
//     for (int ky = 0; ky < this->NbrSiteY; ++ky)      
//     {
//         int TmpIndex = 
//           cout << kx << " " << ky << " " << this->EnergyBandStructure[][] << endl;
//     }
}

// destructor
//

TightBindingModelTrilayerGrapheneTwistBN::~TightBindingModelTrilayerGrapheneTwistBN()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelTrilayerGrapheneTwistBN::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	  int Index = (kx * this->NbrSiteY) + ky;
      
	  HermitianMatrix TmpOneBodyHamiltonian = this->ComputeBlochHamiltonian(((double) (kx + this->GammaX)), ((double)(ky + + this->GammaY)));
// 	  	  cout << TmpOneBodyHamiltonian << endl;
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
//             cout << kx << " " << ky << " " << (this->OneBodyBasis[Index]) << endl;
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
          {
		this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
          }
          
	    }
	}
    }
}



// compute the Bloch hamiltonian at a point of the Brillouin zone
//
// kx = momentum along the x axis
// ky = momentum along the x axis
// return value = Bloch hamiltonian

HermitianMatrix TightBindingModelTrilayerGrapheneTwistBN::ComputeBlochHamiltonian(double kx, double ky)
{
  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
    
  double KxLayer;
  double KyLayer;
  int TmpIndex;
  int TmpIndex1;
  
  
  double KX = this->Gx1 * this->KxFactor * (kx + this->GammaX) + this->Gx2 * this->KyFactor * (ky + this->GammaY);
  double KY = this->Gy1 * this->KxFactor * (kx + this->GammaX) + this->Gy2 * this->KyFactor * (ky + this->GammaY);
  
  double C0 = -10.13;
  double Cz = -9.01;
  double Cab = 11.34;
  
  double phi_0 = 86.53/ 180.0 * M_PI;
  double phi_z = 8.43/ 180.0 * M_PI;
  double phi_ab = 19.6/ 180.0 * M_PI;
   
    
  for (int kxLayer = 0; kxLayer < this->NbrPointSuperlatticeX; ++kxLayer)
  {
      for (int kyLayer = 0; kyLayer < this->NbrPointSuperlatticeY; ++kyLayer)
      {
          KxLayer = KX + (kxLayer - this->NbrPointSuperlatticeX/2)*this->Gx1 + (kyLayer - this->NbrPointSuperlatticeY/2)*this->Gx2;
          KyLayer = KY + (kxLayer - this->NbrPointSuperlatticeX/2)*this->Gy1 + (kyLayer - this->NbrPointSuperlatticeY/2)*this->Gy2;
          
          
          Complex TmpComplexMomentum = Complex (KxLayer, -KyLayer);
          TmpIndex = this->GetLinearizedMomentumIndexSuperlattice(kxLayer, kyLayer);
          
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex, this->UVoltage);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex + 1, this->UVoltage);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 2, 6*TmpIndex + 2, this->UVoltage / 2.0);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 3, 6*TmpIndex + 3, this->UVoltage / 2.0);
          
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex + 2, this->GammaOne);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 3, 6*TmpIndex + 4, this->GammaOne);
          
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex + 1, this->NNHoppingIntra * TmpComplexMomentum);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 2, 6*TmpIndex + 3, this->NNHoppingIntra * TmpComplexMomentum);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 4, 6*TmpIndex + 5, this->NNHoppingIntra * TmpComplexMomentum);
          
          TmpIndex1 = this->GetLinearizedMomentumIndexSuperlatticeSafe(kxLayer + 1, kyLayer);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1, C0*Phase(phi_0)+Cz*Phase(phi_z));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1 + 1, Cab*Phase(2*M_PI/3.0 - phi_ab));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1, Cab*Phase(2*M_PI/3.0 - phi_ab));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1 + 1, C0*Phase(phi_0)-Cz*Phase(phi_z));
                  
          TmpIndex1 = this->GetLinearizedMomentumIndexSuperlatticeSafe(kxLayer, kyLayer + 1);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1, C0*Phase(-phi_0)+Cz*Phase(-phi_z));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1 + 1, Cab*Phase(phi_ab));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1, Cab*Phase(2*M_PI/3.0 + phi_ab));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1 + 1, C0*Phase(-phi_0)-Cz*Phase(-phi_z));
                    
          TmpIndex1 = this->GetLinearizedMomentumIndexSuperlatticeSafe(kxLayer - 1, kyLayer + 1);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1, C0*Phase(phi_0)+Cz*Phase(phi_z));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1 + 1, Cab*Phase(-phi_ab));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1, Cab*Phase(-(2*M_PI/3.0 + phi_ab)));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1 + 1, C0*Phase(phi_0)-Cz*Phase(phi_z));
          
          TmpIndex1 = this->GetLinearizedMomentumIndexSuperlatticeSafe(kxLayer - 1, kyLayer);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1, C0*Phase(-phi_0)+Cz*Phase(-phi_z));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1 + 1, Cab*Phase(-(2*M_PI/3.0-phi_ab)));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1, Cab*Phase(-(2*M_PI/3.0-phi_ab)));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1 + 1, C0*Phase(-phi_0)-Cz*Phase(-phi_z));
          
          TmpIndex1 = this->GetLinearizedMomentumIndexSuperlatticeSafe(kxLayer, kyLayer - 1);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1, C0*Phase(phi_0)+Cz*Phase(phi_z));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1 + 1, Cab*Phase(-(2*M_PI/3.0 + phi_ab)));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1, Cab*Phase(-phi_ab));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1 + 1, C0*Phase(phi_0)-Cz*Phase(phi_z));
  
          TmpIndex1 = this->GetLinearizedMomentumIndexSuperlatticeSafe(kxLayer + 1, kyLayer - 1);
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1, C0*Phase(-phi_0)+Cz*Phase(-phi_z));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex, 6*TmpIndex1 + 1, Cab*Phase(2*M_PI/3.0 + phi_ab));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1, Cab*Phase(phi_ab));
          TmpOneBodyHamiltonian.SetMatrixElement(6*TmpIndex + 1, 6*TmpIndex1 + 1, C0*Phase(-phi_0)-Cz*Phase(-phi_z));
          
      }
  }
  
  return TmpOneBodyHamiltonian;
}

// compute the form factor for the density operator 
// 
// kx = momentum along x of annihilation operator
// ky = momentum along y of creation operator
// qx = momentum transfer along x direction
// qy = momentum transfer along y direction
// valleyIndex = valley index of density operator

Complex TightBindingModelTrilayerGrapheneTwistBN::ComputeDensityFormFactor(int kx, int ky, int qx, int qy, int valleyIndex)
{
    Complex Tmp = 0.0;
    int momentumIndex1;
    int momentumIndex2;
    
    int TmpSumX = kx + qx;
    int TmpSumY = ky + qy;
    
    
    if (valleyIndex == 1)
    {
        momentumIndex1 = this->GetLinearizedMomentumIndexSafe(kx, ky);
        momentumIndex2 = this->GetLinearizedMomentumIndexSafe(kx + qx, ky + qy);
    }
    else
    {
        int new_kx = -kx - qx;
        int new_ky = -ky - qy;
        while (new_kx < 0)
            new_kx += this->NbrSiteX;
        while (new_ky < 0)
            new_ky += this->NbrSiteY;
        new_kx %= this->NbrSiteX;
        new_ky %= this->NbrSiteY;
        momentumIndex1 = this->GetLinearizedMomentumIndexSafe(new_kx, new_ky);
        momentumIndex2 = this->GetLinearizedMomentumIndexSafe(new_kx + qx, new_ky + qy);
        
        TmpSumX = new_kx + qx;
        TmpSumY = new_ky + qy;        
    }
    
    int BZIndexX = (abs(TmpSumX) / this->NbrSiteX);
    int BZIndexY = (abs(TmpSumY) / this->NbrSiteY);
      
    if ((TmpSumX) < 0)
    {
        if ((abs(TmpSumX) % this->NbrSiteX) == 0)
            BZIndexX = -BZIndexX;
        else
            BZIndexX = -BZIndexX - 1;
    }
    if ((TmpSumY) < 0)
    {
        if (((abs(TmpSumY)) % this->NbrSiteY) == 0)
            BZIndexY = -BZIndexY;
        else
            BZIndexY = -BZIndexY - 1;
    }
    
    int TmpIndex1;
    int TmpIndex2;
    for (int i = 0; i < this->NbrPointSuperlatticeX; ++i)
    {
        for (int j = 0; j < this->NbrPointSuperlatticeY; ++j)
        {
            TmpIndex1 = 6 * this->GetLinearizedMomentumIndexSuperlattice(i,j);
            TmpIndex2 = 6 * this->GetLinearizedMomentumIndexSuperlatticeSafe(i + BZIndexX, j + BZIndexY);
            
            for (int k = 0; k < 6; ++k)
                Tmp += Conj(this->OneBodyBasis[momentumIndex1][this->BandIndex][TmpIndex1 + k]) * this->OneBodyBasis[momentumIndex2][this->BandIndex][TmpIndex2 +k];            
            
        }
    }
    
    return Tmp;
}
