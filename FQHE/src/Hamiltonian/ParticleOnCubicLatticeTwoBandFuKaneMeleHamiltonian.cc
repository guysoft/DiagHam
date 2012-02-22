////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 3d topological insulator based on the Fu-Kane-Mele     //
//                       model and restricted to two bands                    //
//                                                                            //
//                        last modification : 18/07/2011                      //
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
#include "Hamiltonian/ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// uPotential = strength of the repulsive two body neareast neighbor interaction
// vPotential = strength of the repulsive two body on site interaction
// nnHopingDistortion111 = distortion of nearest neighbor hoping amplitude in the (111) direction
// spinOrbitCoupling = amplitude of the spin orbit coupling
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// gammaZ = boundary condition twisting angle along z
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, int nbrSiteZ, double uPotential, double vPotential, double nnHopingDistortion111, double spinOrbitCoupling, double gammaX, double gammaY, double gammaZ, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteYZ = this->NbrSiteY * this->NbrSiteZ;
  this->LzMax = nbrSiteX * nbrSiteY * nbrSiteZ - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->KzFactor = 2.0 * M_PI / ((double) this->NbrSiteZ);
  this->HamiltonianShift = 0.0;
  this->NNHopingDistortion111 = nnHopingDistortion111;
  this->SpinOrbitCoupling = spinOrbitCoupling;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->GammaZ = gammaZ;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1 << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1 << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::~ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian()
{
}
  

// compute the matrix element for the two body interaction between two sites A and B with up spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kz2 = momentum along z for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// kz4 = momentum along z for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return 0.0;
  Complex Tmp = 1.0 ;
//   Tmp += Phase (0.5 * ((((double) (kx1 - kx3)) * this->KxFactor) + (((double) (ky1 - ky3)) * this->KyFactor)));
//   Tmp += Phase (0.5 * ((((double) (kx1 - kx3)) * this->KxFactor) + (((double) (kz1 - kz3)) * this->KzFactor)));
//   Tmp += Phase (0.5 * ((((double) (ky1 - ky3)) * this->KyFactor) + (((double) (kz1 - kz3)) * this->KzFactor)));
//   Tmp *= Phase (-0.25 * ((((double) (kx4 - kx2)) * this->KxFactor) + (((double) (ky4 - ky2)) * this->KyFactor)
// 			+ (((double) (kz4 - kz2)) * this->KzFactor)));
  Tmp += Phase ((((double) (kx1 - kx3)) * this->KxFactor));
  Tmp += Phase ((((double) (kz1 - kz3)) * this->KzFactor));
  Tmp += Phase ((((double) (ky1 - ky3)) * this->KyFactor));
  Tmp *= Phase (0.25 * ((((double) (kx4 - kx2)) * this->KxFactor) + (((double) (ky4 - ky2)) * this->KyFactor)
		       + (((double) (kz4 - kz2)) * this->KzFactor)));
//   Tmp *= Phase (0.25 * ((((double) (kx4 - kx2)) * this->KxFactor) + (((double) (ky4 - ky2)) * this->KyFactor)
// 			+ (((double) (kz4 - kz2)) * this->KzFactor)));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and B with down spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kz1 = momentum along z for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kz3 = momentum along z for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kz1 = momentum along z for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kz2 = momentum along z for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kz3 = momentum along z for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// kz4 = momentum along z for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementADownBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites A with opposite spins 
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on A site with spin down
// ky2 = momentum along y for the creation operator on A site with spin down
// kz2 = momentum along z for the creation operator on A site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on A site with spin down
// ky4 = momentum along y for the annihilation operator on A site with spin down
// kz4 = momentum along z for the annihilation operator on A site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  Complex Tmp = 1.0;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B with opposite spins 
//
// kx1 = momentum along x for the creation operator on B site with spin up
// ky1 = momentum along y for the creation operator on B site with spin up
// kz1 = momentum along z for the creation operator on B site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on B site with spin up
// ky3 = momentum along y for the annihilation operator on B site with spin up
// kz3 = momentum along z for the annihilation operator on B site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return 0.0;
//   Complex Tmp = Phase (0.25  * ((((double) (kx4 + kx3 - kx1 - kx2)) * this->KxFactor) + (((double) (ky4 + ky3 - ky1 - ky2)) * this->KyFactor) + (((double) (kz4 + kz3 - kz1 - kz2)) * this->KzFactor)));
  Complex Tmp = Phase (0.25  * ((((double) (kx4 + kx3 - kx1 - kx2)) * this->KxFactor) + 
				(((double) (ky4 + ky3 - ky1 - ky2)) * this->KyFactor) + 
				(((double) (kz4 + kz3 - kz1 - kz2)) * this->KzFactor)));
  return Tmp;
}

// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices

void ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      for (int kz = 0; kz < this->NbrSiteZ; ++kz)
	{
	  HermitianMatrix TmpOneBodyHamiltonian(4, true);
	  int Index = ((kx * this->NbrSiteY) + ky) * this->NbrSiteZ + kz;
	  Complex B1 = 1.0 + this->NNHopingDistortion111 + Phase(((double) kx) * this->KxFactor) + Phase(((double) ky) * this->KyFactor) + Phase(((double) kz) * this->KzFactor);
	  // 	  Complex TmpPhase = Phase (-0.25 * ((((double) kx) * this->KxFactor) +(((double) ky) * this->KyFactor) + (((double) kz) * this->KzFactor)));
	  // 	  Complex TmpPhase = Phase (-0.5 * ((((double) kx) * this->KxFactor) +(((double) ky) * this->KyFactor) + (((double) kz) * this->KzFactor)));
	  // 	  B1 *= TmpPhase;


	  //	  Complex B1 = 1.0 + this->NNHopingDistortion111 + Phase(0.5 * (((double) ky) * this->KyFactor) + 0.5 * (((double) kz) * this->KzFactor))   + Phase(0.5 * (((double) kx) * this->KxFactor) + 0.5 * (((double) kz) * this->KzFactor))  + Phase(0.5 * (((double) kx) * this->KxFactor) + 0.5 * (((double) ky) * this->KyFactor)) ;
// 	  double d3 = this->SpinOrbitCoupling * (sin (0.5 * (((double) kx) * this->KxFactor) + 0.5 * (((double) kz) * this->KzFactor))
// 						 - sin (0.5 * (((double) kx) * this->KxFactor) + 0.5 * (((double) ky) * this->KyFactor))
// 						 - sin (0.5 * (((double) kx) * this->KxFactor) - 0.5 * (((double) ky) * this->KyFactor))
// 						 + sin (0.5 * (((double) kx) * this->KxFactor) - 0.5 * (((double) kz) * this->KzFactor)));
// 	  double d4 = this->SpinOrbitCoupling * (sin (0.5 * (((double) kx) * this->KxFactor) + 0.5 * (((double) kz) * this->KzFactor))
// 						 - sin (0.5 * (((double) ky) * this->KyFactor) + 0.5 * (((double) kz) * this->KzFactor))
// 						 - sin (0.5 * (((double) ky) * this->KyFactor) - 0.5 * (((double) kz) * this->KzFactor))
// 						 + sin (0.5 * (((double) ky) * this->KyFactor) - 0.5 * (((double) kx) * this->KxFactor)));
// 	  double d5 = this->SpinOrbitCoupling * (sin (0.5 * (((double) ky) * this->KyFactor) + 0.5 * (((double) kx) * this->KxFactor))
// 						 - sin (0.5 * (((double) kx) * this->KxFactor) + 0.5 * (((double) kz) * this->KzFactor))
// 						 - sin (0.5 * (((double) kz) * this->KzFactor) - 0.5 * (((double) kx) * this->KxFactor))
// 						 + sin (0.5 * (((double) kz) * this->KzFactor) - 0.5 * (((double) ky) * this->KyFactor)));
	  double d3 = this->SpinOrbitCoupling * (sin ((((double) ky) * this->KyFactor))
						 - sin ((((double) kz) * this->KzFactor))
						 - sin ((((double) ky) * this->KyFactor) - (((double) kx) * this->KxFactor))
						 + sin ((((double) kz) * this->KzFactor) - (((double) kx) * this->KxFactor)));
	  double d4 = this->SpinOrbitCoupling * (sin ((((double) kz) * this->KzFactor))
						 - sin ((((double) kx) * this->KxFactor))
						 - sin ((((double) kz) * this->KzFactor) - (((double) ky) * this->KyFactor))
						 + sin ((((double) kx) * this->KxFactor) - (((double) ky) * this->KyFactor)));
	  double d5 = this->SpinOrbitCoupling * (sin ((((double) kx) * this->KxFactor))
						 - sin ((((double) ky) * this->KyFactor))
						 - sin ((((double) kx) * this->KxFactor) - (((double) kz) * this->KzFactor))
						 + sin ((((double) ky) * this->KyFactor) - (((double) kz) * this->KzFactor)));
	  Complex B2 = d3 - I() * d4;
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d5);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -d5);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 2, -d5);
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 3, d5);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 3, B1);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 2, B2);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 3, -B2);

// 	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, 0.0);
// 	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, 0.0);
// 	  TmpOneBodyHamiltonian.SetMatrixElement(2, 2, 1.0);
// 	  TmpOneBodyHamiltonian.SetMatrixElement(3, 3, 1.0);


	  //	  cout << TmpOneBodyHamiltonian << endl;

	  ComplexMatrix TmpMatrix(4, 4, true);
	  TmpMatrix[0][0] = 1.0;
	  TmpMatrix[1][1] = 1.0;
	  TmpMatrix[2][2] = 1.0;
	  TmpMatrix[3][3] = 1.0;
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	  oneBodyBasis[Index] = TmpMatrix;	
	  if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactorsupup[Index] = TmpDiag(0, 0);
	      this->OneBodyInteractionFactorsdowndown[Index] = TmpDiag(1, 1);
	    }
	  cout << "Index = " << Index << endl;
	  cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << " " << TmpDiag(3, 3) << endl;
	  //	  cout << endl << TmpMatrix << endl;
	}
}
