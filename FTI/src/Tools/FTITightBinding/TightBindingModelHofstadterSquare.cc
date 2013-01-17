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
#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
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
TightBindingModelHofstadterSquare::TightBindingModelHofstadterSquare(int nbrCellX, int nbrCellY, int unitCellX, int unitCellY, int nbrFlux, char axis,
								     double gammaX, double gammaY, 
								     AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->UnitCellX = unitCellX;
  this->UnitCellY = unitCellY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = UnitCellX*UnitCellY;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->LandauGaugeAxis=axis;
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
  this->SetNbrFluxQuanta(nbrFlux);

  this->ComputeBandStructure();
  
}

// destructor
//

TightBindingModelHofstadterSquare::~TightBindingModelHofstadterSquare()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelHofstadterSquare::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double K1;
  double K2;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int Index = this->GetLinearizedMomentumIndex(kx, ky);
	  if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	    {
	      K1 = this->KxFactor*(((double) kx) + this->GammaX);
	      K2 = this->KyFactor*(((double) ky) + this->GammaY);

	      // construct magnetic unit cell:
	      
	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      
	      Complex TranslationPhase;
	      switch (this->LandauGaugeAxis)
		{
		case 'y': {
		  for (int i=0; i<UnitCellX; ++i)
		    {
		      Complex Phase=Polar(1.0,2.0*M_PI*this->FluxDensity*(double)i);
		      for (int j=0; j<UnitCellY; ++j)
			{
			  int InitialIndex = this->EncodeSublatticeIndex(i, j, K1, K2, TranslationPhase); // TranlationPhase always one, so can be discarded
			  int FinalIndex = this->EncodeSublatticeIndex(i+1, j, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif

			  FinalIndex = this->EncodeSublatticeIndex(i-1, j, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif
			  
			  FinalIndex = this->EncodeSublatticeIndex(i, j+1, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Conj(Phase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase*Conj(Phase)<<endl;
#endif

			  FinalIndex = this->EncodeSublatticeIndex(i, j-1, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Phase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase*Phase<<endl;
#endif

			}
		    }
		  break;
		}


		case 'x': {
		  for (int j=0; j<UnitCellY; ++j)
		    {
		      Complex Phase=Polar(1.0,-2.0*M_PI*this->FluxDensity*(double)j);
		      for (int i=0; i<UnitCellX; ++i)
			{

			  int InitialIndex = this->EncodeSublatticeIndex(i, j, K1, K2, TranslationPhase); // TranlationPhase always one, so can be discarded
			  int FinalIndex = this->EncodeSublatticeIndex(i+1, j, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -Conj(Phase)*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-Conj(Phase)*TranslationPhase<<endl;
#endif

			  FinalIndex = this->EncodeSublatticeIndex(i-1, j, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -Phase*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-Phase*TranslationPhase<<endl;
#endif
			  
			  FinalIndex = this->EncodeSublatticeIndex(i, j+1, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif

			  FinalIndex = this ->EncodeSublatticeIndex(i, j-1, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif

			}
		    }
		  break;
		}
		default:
		  cout << "Invalid Landau quantization axis encountered in ParticleOnLatticeDeltaHamiltonian."<<endl;
		  exit(1);
		  break;
		}

#ifdef DEBUG_OUTPUT
	      cout << TmpOneBodyHamiltonian<< endl;
#endif
	    

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



// nbrFluxQuanta = number of quanta of flux piercing the unit cell
void TightBindingModelHofstadterSquare::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->FluxDensity = ((double)NbrFluxQuanta)/this->NbrBands;
  switch (this->LandauGaugeAxis)
    {
    case 'x':
      this->LxTranslationPhase = 1.0;  // no phase for translating in the y-direction in Landau gauge...
      this->LyTranslationPhase = Polar(1.0, 2.0*M_PI*FluxDensity*this->UnitCellY);
      break;
    case 'y':
      this->LxTranslationPhase = Polar(1.0, -2.0*M_PI*FluxDensity*this->UnitCellX);
      this->LyTranslationPhase = 1.0;  // no phase for translating in the y-direction in Landau gauge...
      break;
    default:
      cout << "Unknown Quantization axis! Exiting TightBindingModelHofstadterSquare..."<<endl;
      exit(1);
      break;
    }
}
      
 
// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// KX = current momentum in x-direction
// KY = current momentum in y-direction
// translationPhase = phase factor associated with any crossings of unit cell boundary
//
int TightBindingModelHofstadterSquare::EncodeSublatticeIndex(int posx, int posy, double KX, double KY, Complex &translationPhase)
{
  //cout << "Encoding " << posx<<", "<<posy<<": ";
  int numXTranslations=0, numYTranslations=0;  
  while (posx<0)
    {
      posx+=this->UnitCellX;
      ++numXTranslations;      
    }
  while (posx>=this->UnitCellX)
    {
      posx-=this->UnitCellX;
      --numXTranslations;
    }
  while (posy<0)
    {
      posy+=this->UnitCellY;
      ++numYTranslations;      
    }
  while (posy>=this->UnitCellY)
    {
      posy-=this->UnitCellY;
      --numYTranslations;
    }
  int rst = posx + this->UnitCellX*posy;
  // determine phase for shifting site to the simulation cell:
  Complex tmpPhase(1.0,0.0);
  Complex tmpPhase2;
  translationPhase=tmpPhase;
  if (numXTranslations>0)
    tmpPhase2=LxTranslationPhase;
  else
    tmpPhase2=Conj(LxTranslationPhase);
  for (int i=0; i<abs(numXTranslations); ++i)
    tmpPhase*=tmpPhase2;
  //cout<<" tmpPhaseX="<<tmpPhase;
  for (int y=1;y<=posy; ++y)
    translationPhase*=tmpPhase;
  translationPhase*=Polar(1.0, KX*numXTranslations);
  tmpPhase=1.0;
  if (numYTranslations>0)
    tmpPhase2=LyTranslationPhase;
  else
    tmpPhase2=Conj(LyTranslationPhase);
  for (int i=0; i<abs(numYTranslations); ++i)
    tmpPhase*=tmpPhase2;
  //cout<<" tmpPhaseY="<<tmpPhase;
  for (int x=1;x<=posx; ++x)
    translationPhase*=tmpPhase;
  translationPhase*=Polar(1.0, KY*numYTranslations);
  //cout << "tX="<<numXTranslations<< ", tY="<<numYTranslations<<", translationPhase= " <<translationPhase<<endl;
  return rst;
}
