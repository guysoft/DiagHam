////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Gunnar Möller                         //
//                                                                            //
//  class of tight binding model for the square lattice with homogeneous flux //
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
#include "Tools/FTITightBinding/TightBindingModelOFLNOrbitalTriangularLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>
using std::cout;
using std::endl;

//#define DEBUG_OUTPUT


// default constructor
//
// nbrCellsX = number of unit cells in the x direction
// nbrCellsY = number of unit cella in the y direction
// unitCellX = number of sites in unit cell in x direction
// unitCellY = number of sites in unit cell in y direction
// nbrFlux = number of flux quanta per unit cell
// axis = direction of Landau gauge within cell ('x' or 'y')
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
TightBindingModelOFLNOrbitalTriangularLattice::TightBindingModelOFLNOrbitalTriangularLattice(double laserStrength,  int nbrInternalDegree,  int nbrSiteX,  int nbrSiteY, double gammaX, double gammaY, AbstractArchitecture* architecture, int cutOFF, bool storeOneBodyMatrices)
{
  this->NbrStep = cutOFF+1;
  this->NbrInternalDegree =  nbrInternalDegree;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LaserStrength = laserStrength; 
  this->InvMomentum= 1.0/((2.0 * M_PI)*(2.0 * M_PI));
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = this->NbrInternalDegree*this->NbrStep*this->NbrStep;
  this->NbrStatePerBand =  this->NbrSiteX *  this->NbrSiteY;
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

TightBindingModelOFLNOrbitalTriangularLattice::~TightBindingModelOFLNOrbitalTriangularLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelOFLNOrbitalTriangularLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double K1;
  double K2;
  for (int kx = 0; kx < this->NbrSiteX ; ++kx)
    {
      for (int ky = 0; ky <  this->NbrSiteY; ++ky)
	{
	  int Index = this->GetLinearizedMomentumIndex(kx, ky);
	  if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	    {
	      K1 = this->KxFactor*(((double) kx) + this->GammaX);
	      K2 = this->KyFactor*(((double) ky) + this->GammaY);
	      
	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      cout<<"this->NbrBands"<<this->NbrBands<<endl;
	      for (int p = 0; p< this->NbrStep ; p++)
		{
		  
		  for (int t = 0; t< this->NbrStep; t++)
		    {
		      for (int Spin = 0; Spin < this->NbrInternalDegree; Spin++)
			{
			  double MomentaY = K2-M_PI*(this->NbrStep-1) + 2*M_PI*t;
			  double MomentaX = this->NbrInternalDegree*(-M_PI*(this->NbrStep -1) + 2*M_PI*p + K1) + 2*M_PI*Spin;
			  int IntermediateIndex = this->GetIntermediateLinearizedIndices(p, t, Spin);
			  TmpOneBodyHamiltonian.AddToMatrixElement(IntermediateIndex, IntermediateIndex,this->InvMomentum*(MomentaX*MomentaX+MomentaY*MomentaY-MomentaX*MomentaY ));
			  
			  
			  int IntermediateIndex1 = this->GetIntermediateLinearizedIndices(p, t+1, Spin);
			  if (IntermediateIndex1 < IntermediateIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(IntermediateIndex1, IntermediateIndex,-LaserStrength*Phase(2*M_PI*Spin/((double) this->NbrInternalDegree ))); 
			  
			  IntermediateIndex1=this->GetIntermediateLinearizedIndices(p, t-1, Spin);
			  if (IntermediateIndex1 < IntermediateIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(IntermediateIndex1, IntermediateIndex,-LaserStrength*Phase(-2*M_PI*Spin/((double) this->NbrInternalDegree ))); 		 
			  
			  //
			  
			  IntermediateIndex1 = this->GetIntermediateLinearizedIndices(p, t, Spin-1);
			  if (IntermediateIndex1 < IntermediateIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(IntermediateIndex1, IntermediateIndex,-LaserStrength);
			  
			  IntermediateIndex1 = this->GetIntermediateLinearizedIndices(p, t-1, Spin-1);
			  if (IntermediateIndex1 < IntermediateIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(IntermediateIndex1, IntermediateIndex,-LaserStrength*Phase(-M_PI*(2*Spin-1)/((double) this->NbrInternalDegree ))); 


			  IntermediateIndex1 = this->GetIntermediateLinearizedIndices(p, t, Spin+1);
			  if (IntermediateIndex1 < IntermediateIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(IntermediateIndex1, IntermediateIndex,-LaserStrength);
			  
			  IntermediateIndex1 = this->GetIntermediateLinearizedIndices(p, t+1, Spin+1);
			  if (IntermediateIndex1 < IntermediateIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(IntermediateIndex1, IntermediateIndex,-LaserStrength*Phase(M_PI*(2*Spin+1)/((double) this->NbrInternalDegree ))); 


			}
		    }
		}

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


// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool TightBindingModelOFLNOrbitalTriangularLattice::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  int LimitOut = 100;
  if (this->NbrBands < LimitOut)
    LimitOut = this->NbrBands;
  File << "# kx    ky";
  //for (int i = 0; i < LimitOut ; ++i)
  //  File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky <  this->NbrSiteY; ++ky)
	{
	  int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(kx, ky);

	  for (int i = 0; i < LimitOut; ++i)
	    File << kx << " " << ky << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex]<<endl;
	  File << endl; 
	}
      File << endl;
    }
  File.close();
  return true;
}


