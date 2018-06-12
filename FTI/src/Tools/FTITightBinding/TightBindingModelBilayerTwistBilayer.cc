////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                         Class author : Cecile Repellin                     //
//                                                                            //
// class of tight binding model for two bilayers with a relative twist angle  //
//                      with valley and spin conservation                     //
//               (must be TR conjugated to obtain other valley index)         //
//                   last modification : 15/05/2018                           //
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
#include "Tools/FTITightBinding/TightBindingModelBilayerTwistBilayer.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

using std::cout;
using std::endl;


//default constructor
//

TightBindingModelBilayerTwistBilayer::TightBindingModelBilayerTwistBilayer()
{
    this->GammaX = 0;
    this->GammaY = 0;
    
    this->NNHoppingIntra = 0;
    this->NNHoppingInter = 0;
    this->GammaOne = 0;
    this->TrigonalWarping = 0;
    this->UVoltage = 0;
    this->TwistingAngle = 0;
}

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

TightBindingModelBilayerTwistBilayer::TightBindingModelBilayerTwistBilayer(int nbrSiteX, int nbrSiteY, int nbrPointsX, int nbrPointsY, double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
   
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NNHoppingIntra = 2400;
  this->NNHoppingInter = 114;
  this->GammaOne = 340;
  this->TrigonalWarping = 0.051*this->NNHoppingIntra;
  this->UVoltage = 20;
  this->TwistingAngle = 0.7 / 180.0 * M_PI;
    
  this->EmbeddingX = RealVector(2*nbrPointsX, true);
  this->EmbeddingY = RealVector(2*nbrPointsY, true);
  
  this->NbrPointSuperlatticeX = nbrPointsX;
  this->NbrPointSuperlatticeY = nbrPointsY;
  
  this->BandIndex = 2 * nbrPointsX * nbrPointsY - 1;
  
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  
  this->KxFactor = 1.0 / ((double) this->NbrSiteX);
  this->KyFactor = 1.0 / ((double) this->NbrSiteY);
  
  
  double latticeConstant = 1.0 / (2.0 * sin(this->TwistingAngle / 2.0));
  this->Gx1 = 2.0 * M_PI / (sqrt(3.0) * latticeConstant);
  this->Gy1 = 2.0 * M_PI / latticeConstant;
  this->Gx2 = -this->Gx1;
  this->Gy2 = this->Gy1;
  
  this->NbrBands = 4*nbrPointsX*nbrPointsY;
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
  
//   cout << this->GetOneBodyMatrix(0) << endl;
}

// destructor
//

TightBindingModelBilayerTwistBilayer::~TightBindingModelBilayerTwistBilayer()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelBilayerTwistBilayer::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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

HermitianMatrix TightBindingModelBilayerTwistBilayer::ComputeBlochHamiltonian(double kx, double ky)
{
  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
    
  double KxLayer;
  double KyLayer;
  double KxLayerA;
  double KyLayerA;
  double KxLayerB;
  double KyLayerB;
  int TmpIndex;
  int TmpIndex1;
  int TmpIndex2;
  
  int num_layer1 = 2;
  int num_layer2 = 2;
  
  
  double KX = this->Gx1 * this->KxFactor * kx + this->Gx2 * this->KyFactor * ky;
  double KY = this->Gy1 * this->KxFactor * kx + this->Gy2 * this->KyFactor * ky;
	  
    
  for (int kxLayer = 0; kxLayer < NbrPointSuperlatticeX; ++kxLayer)
  {
      for (int kyLayer = 0; kyLayer < NbrPointSuperlatticeY; ++kyLayer)
      {
          KxLayer = KX + (kxLayer - this->NbrPointSuperlatticeX/2)*this->Gx1 + (kyLayer - this->NbrPointSuperlatticeY/2)*this->Gx2;
          KyLayer = KY + (kxLayer - this->NbrPointSuperlatticeX/2)*this->Gy1 + (kyLayer - this->NbrPointSuperlatticeY/2)*this->Gy2;
          KxLayerA = KxLayer + this->Gx1;
          KyLayerA = KyLayer + this->Gy1 / 3.0;
          KxLayerB = KxLayer + this->Gx1;
          KyLayerB = KyLayer - this->Gy1 / 3.0;
          Complex TmpComplexMomentumA = Complex(KxLayerA, -KyLayerA);
          Complex TmpComplexMomentumB = Complex(KxLayerB, -KyLayerB);
          
          TmpIndex = this->GetLinearizedMomentumIndexSuperlattice(kxLayer, kyLayer);
          TmpIndex1 = this->GetLinearizedMomentumIndexSuperlatticeSafe(kxLayer+1, kyLayer);
          TmpIndex2 = this->GetLinearizedMomentumIndexSuperlatticeSafe(kxLayer, kyLayer + 1);
          
          Complex TmpA = this->TrigonalWarping * TmpComplexMomentumA * Phase(-this->TwistingAngle / 2.0) - pow(this->NNHoppingIntra, num_layer1) * pow(Conj(TmpComplexMomentumA), num_layer1) / pow(this->GammaOne, num_layer1 - 1) * Phase(-num_layer1 * this->TwistingAngle / 2);
                   
          Complex TmpB = this->TrigonalWarping * TmpComplexMomentumB * Phase(this->TwistingAngle / 2.0) - pow(this->NNHoppingIntra, num_layer2) * pow(Conj(TmpComplexMomentumB), num_layer2) / pow(this->GammaOne, num_layer2 - 1) * Phase(num_layer2 * this->TwistingAngle / 2);
          
          
          // First layer HAA
          TmpOneBodyHamiltonian.SetMatrixElement(2*TmpIndex, 2*TmpIndex, this->UVoltage);
          TmpOneBodyHamiltonian.SetMatrixElement(2*TmpIndex + 1, 2*TmpIndex + 1, 2.0 / 3.0 * this->UVoltage);
          TmpOneBodyHamiltonian.SetMatrixElement(2*TmpIndex, 2*TmpIndex + 1, TmpA);
          
          // Second layer HBB
          TmpOneBodyHamiltonian.SetMatrixElement(this->NbrBands/2 + 2*TmpIndex, this->NbrBands/2 + 2*TmpIndex, this->UVoltage / 3.0);
          TmpOneBodyHamiltonian.SetMatrixElement(this->NbrBands/2 + 2*TmpIndex, this->NbrBands/2 + 2*TmpIndex + 1, TmpB);
          
          // Hybridization terms HAB and HBA
          TmpOneBodyHamiltonian.SetMatrixElement(2*TmpIndex + 1, this->NbrBands/2 + 2*TmpIndex, this->NNHoppingInter);
          TmpOneBodyHamiltonian.SetMatrixElement(2*TmpIndex + 1, this->NbrBands/2 + 2*TmpIndex1, this->NNHoppingInter * Phase(-2.0*M_PI / 3.0));
          TmpOneBodyHamiltonian.SetMatrixElement(2*TmpIndex + 1, this->NbrBands/2 + 2*TmpIndex2, this->NNHoppingInter * Phase(2.0*M_PI / 3.0));
          
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

Complex TightBindingModelBilayerTwistBilayer::ComputeDensityFormFactor(int kx, int ky, int qx, int qy, int valleyIndex)
{
    Complex Tmp = 0.0;
    int momentumIndex1;
    int momentumIndex2;
    int shiftIndex = (this->EmbeddingX.GetVectorDimension() * this->EmbeddingY.GetVectorDimension()) / 2;
    if (valleyIndex == 1)
    {
        momentumIndex1 = this->GetLinearizedMomentumIndex(kx, ky);
        momentumIndex2 = this->GetLinearizedMomentumIndexSafe(kx + qx, ky + qy);
    }
    else
    {
        momentumIndex1 = this->GetLinearizedMomentumIndexSafe(-kx, -ky);
        momentumIndex2 = this->GetLinearizedMomentumIndexSafe(-kx-qx, -ky-qy);
    }
    int BZIndexX = (kx + qx) / this->NbrSiteX;
    int BZIndexY = (ky + qy) / this->NbrSiteY;
    int TmpIndex1;
    int TmpIndex2;
    for (int i = 0; i < this->NbrPointSuperlatticeX; ++i)
    {
        for (int j = 0; j < this->NbrPointSuperlatticeY; ++j)
        {
            TmpIndex1 = 2 * this->GetLinearizedMomentumIndexSuperlattice(i,j);
            TmpIndex2 = 2 * this->GetLinearizedMomentumIndexSuperlatticeSafe(i + BZIndexX, j + BZIndexY);
            
            Tmp += Conj(this->OneBodyBasis[momentumIndex1][this->BandIndex][TmpIndex1]) * this->OneBodyBasis[momentumIndex2][this->BandIndex][TmpIndex2];
            Tmp += Conj(this->OneBodyBasis[momentumIndex1][this->BandIndex][TmpIndex1 + 1]) * this->OneBodyBasis[momentumIndex2][this->BandIndex][TmpIndex2 + 1];
            Tmp += Conj(this->OneBodyBasis[momentumIndex1][this->BandIndex][shiftIndex + TmpIndex1]) * this->OneBodyBasis[momentumIndex2][this->BandIndex][shiftIndex + TmpIndex2];
            Tmp += Conj(this->OneBodyBasis[momentumIndex1][this->BandIndex][shiftIndex + TmpIndex1 + 1]) * this->OneBodyBasis[momentumIndex2][this->BandIndex][shiftIndex + TmpIndex2 + 1];
            
        }
    }
    
    if (valleyIndex == 1)
        return Tmp;
    else
        return Conj(Tmp);
}



// get the high symmetry points 
//
// pointNames = name of each high symmetry point
// pointCoordinates = coordinates in the first Brillouin zone of the high symmetry points
// return value = number of high symmetry points

int TightBindingModelBilayerTwistBilayer::GetHighSymmetryPoints(char**& pointNames, double**& pointCoordinates)
{
  int NbrHighSymmetryPoints = 3;
  pointNames = new char*[NbrHighSymmetryPoints];
  pointCoordinates = new double*[NbrHighSymmetryPoints];

  pointNames[0] = new char[16];
  sprintf (pointNames[0], "Gamma");
  pointCoordinates[0] = new double[2];
  pointCoordinates[0][0] = 0.0; 
  pointCoordinates[0][1] = 0.0; 

  pointNames[1] = new char[16];
  sprintf (pointNames[1], "K");
  pointCoordinates[1] = new double[2];
  pointCoordinates[1][0] = 4.0 * M_PI / 3.0; 
  pointCoordinates[1][1] = 2.0 * M_PI / 3.0; 

  pointNames[2] = new char[16];
  sprintf (pointNames[2], "M");
  pointCoordinates[2] = new double[2];
  pointCoordinates[2][0] = M_PI; 
  pointCoordinates[2][1] = 0.0; 

  return NbrHighSymmetryPoints;
}

// compute the distance between two points in the first Brillouin zone, changing the coordinates the second one by a reciprocal lattice vector if needed
//
// kx1 = momentum of the first point along the x axis
// ky1 = momentum of the first point along the y axis
// kx2 = reference on the momentum of the second point along the x axis
// ky2 = reference on the momentum of the second point along the y axis
// return value = distance between the two points

double TightBindingModelBilayerTwistBilayer::GetDistanceReciprocalSpace(double kx1, double ky1, double& kx2, double& ky2)
{
  double AngleFactor = 2.0 * cos(2.0 * M_PI / 3.0);
  double DiffKx = kx1 - kx2;
  double DiffKy = ky1 - ky2;
  double MinDistance = sqrt ((DiffKx * DiffKx) + (DiffKy * DiffKy) + (AngleFactor * DiffKx * DiffKy));
  double MinKx2 = kx2;
  double MinKy2 = ky2;
  for (int i = -1; i <= 1; ++i)
    {
      double TmpKx2  = kx2 + (2.0 * ((double) i) * M_PI);
      for (int j = -1; j <= 1; ++j)
	{
	  double TmpKy2  = ky2 + (2.0 * ((double) j) * M_PI);	  
	  double DiffKx = kx1 - TmpKx2;
	  double DiffKy = ky1 - TmpKy2;
	  double TmpDistance =  sqrt ((DiffKx * DiffKx) + (DiffKy * DiffKy) + (AngleFactor * DiffKx * DiffKy));
	  if (TmpDistance < MinDistance)
	    {
	      MinDistance = TmpDistance;
	      MinKx2 = TmpKx2;
	      MinKy2 = TmpKy2;
	    }
	}
    }
  kx2 = MinKx2;
  ky2 = MinKy2;
  return MinDistance;
}




// compute the distance between two points in the same unit cell, changing the coordinates the second one by a lattice vector if needed
//
// x1 = x coordinate of the first point
// y1 = y coordinate of the first point
// x2 = reference to the x coordinate of the second point
// y2 = referenceto the y coordinate of the second point
// return value = distance between the two points

double TightBindingModelBilayerTwistBilayer::GetDistanceRealSpace(double x1, double y1, double& x2, double& y2)
{
  double AngleFactor = 2.0 * cos(2.0 * M_PI / 3.0);
  double DiffX = x1 - x2;
  double DiffY = y1 - y2;
  double MinDistance = sqrt ((DiffX * DiffX) + (DiffY * DiffY) + (AngleFactor * DiffX * DiffY));
  double MinX2 = x2;
  double MinY2 = y2;
  for (int i = -1; i <= 1; ++i)
    {
      double TmpX2 = x2 + ((double) i);
      for (int j = -1; j <= 1; ++j)
	{
	  double TmpY2 = y2 + ((double) j);
	  DiffX = x1 - TmpX2;
	  DiffY = y1 - TmpY2;
	  double TmpDistance =  sqrt ((DiffX * DiffX) + (DiffY * DiffY) + (AngleFactor * DiffX * DiffY));
	  if (TmpDistance < MinDistance)
	    {
	      MinDistance = TmpDistance;
	      MinX2 = TmpX2;
	      MinY2 = TmpY2;
	    }
	}
    }
  x2 = MinX2;
  y2 = MinY2;
  return MinDistance;
}
