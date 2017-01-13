////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the Checkerboard lattice       //
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


#ifndef TIGHTBINDINGMODELHOFSTADTERSQUARE_H
#define TIGHTBINDINGMODELHOFSTADTERSQUARE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelHofstadterSquare : public Abstract2DTightBindingModel
{

 protected:


  // axis of Landau gauge:
  char LandauGaugeAxis;
  // number of sites in cell in x-direction
  int UnitCellX;
  // number of sites in cell in y-direction
  int UnitCellY;

  // number of flux quanta in cell (cancelled by opposite flux)
  int NbrFluxQuanta;

  // auxiliary variables:
  // flux density:
  double FluxDensity;
  // magnetic translation phases;
  Complex LxTranslationPhase;
  Complex LyTranslationPhase;

  // flag indicating whether natural embedding is used.
  bool UsingNaturalEmbedding;

 public:

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
  TightBindingModelHofstadterSquare(int nbrCellX, int nbrCellY, int unitCellX, int unitCellY, int nbrFlux, char axis,
				       double gammaX, double gammaY, 
				       AbstractArchitecture* architecture, bool storeOneBodyMatrices = true,  bool useEmbedding = false);

  // destructor
  //
  ~TightBindingModelHofstadterSquare();


  // get the position of a sublattice site
  //
  // position = reference on a vector where the answer is supplied
  // sublatticeIndex = index of the sub-lattice position
  virtual void GetSublatticeVector(RealVector &position, int sublatticeIndex);


  // convert absolute coordinates into lattice coordinates and sublattice index
  //
  // position = coordinates of the site to be identified
  // tx, ty = translations of the site in units of lattice vectors
  // sublatticeIndex = index of the sub-lattice position (if matching a lattice site; -1 if not found)
  // return = true if the coordinates correspond to a lattice site
  virtual bool PositionToLatticeCoordinates(RealVector &position, int &tx, int &ty, int &sublatticeIndex);

  // get the eigenstates in real space, using CoreComputeBandStructureWithEmbedding
  // 
  // return value = tight binding eigenvectors
  ComplexMatrix GetRealSpaceTightBindingEigenstates();


  // get the tight binding hamiltonian in real space 
  // 
  // return value = tight binding hamiltonian
  HermitianMatrix GetRealSpaceTightBindingHamiltonian();


  HermitianMatrix  BuildTightBindingHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes);

  void ComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
				  int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
				  bool bosonFlag, double uPotential, double vPotential);

  // returns the single-particle wavefunction according to the definition phi_{n,k}=u_{n,alpha}(k)*exp{i k.r}
  // deducing the appropriate sublattice value alpha from the position
  //
  // FunctionValue[out] = return value phi_{nk}(r)
  // Position = overall position vector relative to the origin
  // indexK = linearised index of momentum sector
  // bandIndex = band index
  // return value = single-particle wavefunction in first argument
  void GetFunctionValue(Complex& FunctionValue, RealVector& Position, int indexK, int bandIndex);

  // returns the single-particle wavefunction according to the definition phi_{n,k}=u_{n,alpha}(k)*exp{i k.r}
  // defining position in multiples of lattice vectors and sublattice index
  // 
  // FunctionValue[out] = return value phi_{nk}(r)
  // Rx, Ry, alpha = defining position via r=Rx UnitCellX ex + Ry UnitCellY ey + rho_alpha
  // indexK = linearised index of momentum sector
  // bandIndex = band index
  // return value = single-particle wavefunction in first argument
  void GetFunctionValue(Complex& FunctionValue, int Rx, int Ry, int alpha, int indexK, int bandIndex);

  // decode single integer for sublattice index into set of quantum numbers/positions posx, posy
  // index = sublattice index
   // [out] posx = position along x-direction
  // [out] posy = position along y-direction
  //
  void DecodeSublatticeIndex(int index, int &posx, int &posy);  

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

  // version with real-space embedding of the wavefunctions
  void CoreComputeBandStructureWithEmbedding(long minStateIndex, long nbrStates);

  // initialize number of flux quanta
  // nbrFluxQuanta = number of quanta of flux piercing the unit cell
  void SetNbrFluxQuanta(int nbrFluxQuanta);

  // set natural embedding, i.e., at positions of a uniform lattice
  //
  void SetNaturalEmbedding();

  // code set of quantum numbers posx, posy into a single integer
  // posx = position along x-direction
  // posy = position along y-direction
  // KX = current momentum in x-direction
  // KY = current momentum in y-direction
  // translationPhase = phase factor associated with any crossings of unit cell boundary
  //
  int EncodeSublatticeIndex(int posx, int posy, double KX, double KY, Complex &translationPhase);


  // code set of quantum numbers posx, posy into a single integer
  // posx = position along x-direction
  // posy = position along y-direction
  // numXTranslations = number of translation in the x direction to get back to the unit cell 
  // numXTranslations = number of translation in the y direction to get back to the unit cell
  //
  int EncodeSublatticeIndex(int posx, int posy,int & numXTranslations,int &numYTranslations, Complex &translationPhase);

  int GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int orbitalIndex, int & numXTranslations, int &numYTranslations);

  // obtain dimensions of magnetic unit cell
  // numX = number of unit cells within MUC along x-direction
  // numY = number of unit cells within MUC along y-direction
  void GetMUCDimensions(int &numX, int &numY);

};


// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// numXTranslations = number of translation in the x direction to get back to the unit cell 
// numXTranslations = number of translation in the y direction to get back to the unit cell
//
inline int TightBindingModelHofstadterSquare::EncodeSublatticeIndex(int posx, int posy,int & numXTranslations,int &numYTranslations, Complex &translationPhase) 
{
  numXTranslations=0;
  numYTranslations=0;

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
  Complex tmpPhase(1.0,0.0);
  Complex tmpPhase2;
  translationPhase=tmpPhase;
  if (numXTranslations>0)
    tmpPhase2=LxTranslationPhase;
  else
    tmpPhase2=Conj(LxTranslationPhase);
  for (int i=0; i<abs(numXTranslations); ++i)
    tmpPhase*=tmpPhase2;
  tmpPhase=1.0;
  if (numYTranslations>0)
    tmpPhase2=LyTranslationPhase;
  else
    tmpPhase2=Conj(LyTranslationPhase);
  for (int i=0; i<abs(numYTranslations); ++i)
    tmpPhase*=tmpPhase2;
  return posx + this->UnitCellX*posy;
}


// decode single integer for sublattice index into set of quantum numbers/positions posx, posy
// index = sublattice index
// [out] posx = position along x-direction
// [out] posy = position along y-direction
//
inline void TightBindingModelHofstadterSquare::DecodeSublatticeIndex(int index, int &posx, int &posy)
{
  posx = index % this->UnitCellX;
  posy = index / this->UnitCellX;
}


// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  
//
inline int  TightBindingModelHofstadterSquare::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int orbitalIndex, int & numXTranslations, int &numYTranslations)
{
  numXTranslations=0;
  numYTranslations=0;
  orbitalIndex %= this->NbrBands;
  if (orbitalIndex < 0)
    orbitalIndex +=  this->NbrBands;

  if(x >= this->NbrSiteX)
  {
    x -=  this->NbrSiteX;
    numXTranslations--;
  }
  if (x < 0)
  {
    x +=  this->NbrSiteX;
    numXTranslations++;
  }
  if(y >= this->NbrSiteY)
  {
    y -=  this->NbrSiteY;
    numYTranslations--;
  }
  if (y < 0)
  {
    y +=  this->NbrSiteY;
    numYTranslations++;
  }
  return this->GetRealSpaceTightBindingLinearizedIndex(x, y, orbitalIndex);
}

// obtain dimensions of magnetic unit cell for case of Hofstadter model
// numX = number of unit cells within MUC along x-direction
// numY = number of unit cells within MUC along y-direction
inline void TightBindingModelHofstadterSquare::GetMUCDimensions(int &numX, int &numY)
{
    numX=this->UnitCellX;
    numY=this->UnitCellY;
}

#endif
