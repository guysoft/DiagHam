////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 2D tight binding model                  //
//                                                                            //
//                        last modification : 01/05/2012                      //
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


#ifndef ABSTRACT2DTIGHTBINDINGMODEL_H
#define ABSTRACT2DTIGHTBINDINGMODEL_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract1DTightBindingModel.h"
#include "Matrix/RealMatrix.h"


class Abstract2DTightBindingModel : public Abstract1DTightBindingModel
{

 protected:

   // number of sites in the y direction
  int NbrSiteY;

  // numerical factor for momentum along y
  double KyFactor;

  // boundary condition twisting angle along y
  double GammaY;

  // embedding of sublattices relative to the unit cell reference point along y
  RealVector EmbeddingY;

  // angle between the two primitive vectors
  double TwistAngle;

  // Chern number of each band
  int* Chern;

  // Berry curvature of each band over the BZ
  RealMatrix* Curvature;

  double *LLLGammaX;
  double *LLLGammaY;
    
  //first coordinate of the first spanning vector for a tilted lattice
  int Nx1;
  //second coordinate of the first spanning vector for a tilted lattice
  int Ny1;
  //first coordinate of the second spanning vector for a tilted lattice
  int Nx2;
  //second coordinate of the second spanning vector for a tilted lattice
  int Ny2;
  //array of projected momenta
  double** ProjectedMomenta;
  //second coordinate in momentum space of the second spanning vector of the reciprocal lattice for a tilted lattice
  int Offset;
  //second coordinate in real space of the second spanning vector of the direct lattice for a tilted lattice
  int OffsetReal;

  // unitary & involutory matrix where the (a, b) element Iab appears in: I|x,y,a> = sum_b Iab |-x-d_ax,-y-d_ay,b>
  ComplexMatrix Inversion;

 public:

  // default constructor
  //
  Abstract2DTightBindingModel();

  // destructor
  //
  ~Abstract2DTightBindingModel();

  // get the linearized momentum index
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // return value = linearized momentum index
  virtual int GetLinearizedMomentumIndex(int kx, int ky);

  // get momentum value from a linearized momentum index
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // return value = linearized momentum index
  virtual void GetLinearizedMomentumIndex(int index, int& kx, int& ky);

  // get the linearized momentum index, without assuming k to be in the first BZ
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // return value = linearized momentum index
  virtual int GetLinearizedMomentumIndexSafe(int kx, int ky);

  // get momentum value from a linearized momentum index, without assuming k to be in the first BZ
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // return value = inearized momentum index
  virtual void GetLinearizedMomentumIndexSafe(int index, int& kx, int& ky);

  // get the index of the real space tight binding model from the real space coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndex(int x, int y, int orbitalIndex);
  
  // get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int orbitalIndex);

  // get the real space coordinates from the index of the real space tight binding model
  //
  // index = linearized index of the real space tight binding model
  // x = reference on the x coordinate of the unit cell
  // y = reference on the y coordinate of the unit cell
  // orbitalIndex = reference on the index of the orbital / site within the unit cell
  virtual void GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y, int& orbitalIndex);

  // get the angle between the two primitive lattice vectors
  //
  // return value = angle between the two primitive lattice vectors
  virtual double GetTwistAngle();

  // get the number of sites in the y direction
  //
  // return value = number of sites in the y direction
  int GetNbrSiteY();

  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

  // write the full band structure information in an ASCII file
  //
  // fileName = name of the output file 
  // return value = true if no error occured  
  virtual bool WriteBandStructureASCII(char* fileName);

  // compute the exponentiated, unitary Abelian connection
  //
  // kx = momentum along x
  // ky = momentum along y
  // qx = momentum transfer along x
  // qy = momentum transfer along y
  // band = band index
  // return value = Phase of < u(k) | u(k+q) >
  virtual Complex GetAbelianConnection(int kx, int ky, int qx, int qy, int band);

  // compute the exponentiated, unitary Abelian connection times the quantum distance
  //
  // kx = momentum along x
  // ky = momentum along y
  // qx = momentum transfer along x
  // qy = momentum transfer along y
  // band = band index
  // return value = < u(k) | u(k+q) >
  virtual Complex GetAbelianConnectionQuantumDistance(int kx, int ky, int qx, int qy, int band);

  // compute the unitary Abelian Wilson loop
  //
  // ky = momentum along y
  // band = band index
  // return value = value of the Wilson loop
  virtual Complex GetAbelianWilsonLoopX(int ky, int band);

  // compute the unitary Abelian Wilson loop
  //
  // kx = momentum along x
  // band = band index
  // return value = value of the Wilson loop
  virtual Complex GetAbelianWilsonLoopY(int kx, int band);

  // get the total curvature over the plaquette specified by the lower-left corner. C > 0 corresponds to Bz < 0 magnetic field (holomorphic z = x + iy)
  //
  // kx = momentum along x
  // ky = momentum along y
  // band = band index
  // return value = Arg(<k|k+x><k+x|k+x+y><k+x+y|k+y><k+y|k>) / (2 Pi), range = (-0.5, 0.5]
  virtual double GetCurvature(int kx, int ky, int band);

  // get the Chern number of a specified band
  //
  // band = band index
  // return = Chern number
  virtual int GetChernNumber(int band);

  // compute the Chern number of several bands
  //
  // bands = band indices
  // nbrBands = number of bands that have to be taken into account
  // return value = Chern number
  virtual double ComputeChernNumber(int* bands, int nbrBands);

  // get the LLL twist angle along x for Bloch construction gauge fixing
  //
  // band = band index
  // return = gamma_x
  virtual double GetLLLGammaX(int band);

  // get the LLL twist angle along y for Bloch construction gauge fixing
  //
  // band = band index
  // return = gamma_y
  virtual double GetLLLGammaY(int band);

  // compute the total curvature over the plaquette specified by the lower-left corner. C > 0 corresponds to Bz < 0 magnetic field (holomorphic z = x + iy)
  //
  // kx = momentum along x
  // ky = momentum along y
  // band = band index
  // return value = Arg(<k|k+x><k+x|k+x+y><k+x+y|k+y><k+y|k>) / (2 Pi), range = (-0.5, 0.5]
  virtual double ComputeCurvatureSinglePlaquette(int kx, int ky, int band);

  // compute the stream function for the part of Berry connections that accounts for curvature fluctuations
  //
  // band = band index
  // phi = reference to the stream function over linearized BZ
  // vx = reference to the linear coefficent in phi along ky with a minus sign
  // vy = reference to the linear coefficent in phi along kx
  // return = 0 if succeed, otherwise fail
  virtual int ComputeStreamFunction(int band, RealVector& phi, double& vx, double& vy);

  // build the gauge transform such that gauge(k) * |k⟩_lat is in the "Γ"-shaped parallel-transport gauge
  // 
  // band = band index
  // gauge = reference to the gauge transform
  virtual void BuildParallelTransportGauge(int band, ComplexMatrix& gauge);

  // build the gauge transform such that gauge(k) * |k⟩_lat is in the generalized π/2-rotated Landau gauge
  //
  // band = band index
  // gauge = reference to the gauge transform
  // return = 0 if succeed, otherwise fail
  virtual int BuildGeneralizedLandauGauge(int band, ComplexMatrix& gauge);

  // compute the curvature over each plaquette in the BZ, and also Chern number
  //
  // band = band index
  // return = 0 if succeed, otherwise fail
  virtual void ComputeCurvature();

  // compute the Chern number of a given band
  //
  // band = band index
  // return value = Chern number
  virtual double ComputeChernNumber(int band);

  // compute the Berry curvature  of a given band
  //
  // band = band index
  // fileName = name of the output file 
  // return value = Chern number
  virtual double ComputeBerryCurvature(int band, char* fileName);
  
  
  //compute the complex eigenvalues of the D(ky) matrix (in order to compute the Z2 invariant)
  //
  //bandIndex = band index (corresponds to two bands that are related by time reversal symmetry)
  //nbrOccupiedBands = dimension of the D matrix
  //DMatrixEigenvalues = array of eigenvalues of the D Matrix, for all values of ky
  //kyMin = minimal value of ky for which the D matrix has to be diagonalized
  //kyMax = maximal value of ky for which the D matrix has to be diagonalized
  //nbrKy = number of ky values for which the D matrix has to be diagonalized
  //return value = array of eigenvalues of the D Matrix
  virtual Complex** ComputeDMatrixEigenvalues(int nbrOccupiedBands, int kyMin, int kyMax, int nbrKy);
  
  // write the eigenvalues of the D matrix in an ASCII file
  //
  // fileName = name of the ASCII file 
  //nbrOccupiedBands = nbr of occupied bands
  // return value = true if no error occured
  virtual bool WriteAsciiDMatrixEigenValues(char* fileName, int nbrOccupiedBands);
  
  //compute the Z2 topological invariant for a system with time reversal symmetry
  //
  //nbrOccupiedBands = number of occupied bands
  //return value = Z2 invariant
  virtual int ComputeZ2Invariant(int nbrOccupiedBands);
  
  //Computes value of projected momentum along the lattice directions
  //
  //kx = first coordinate of the given point in the Brillouin zone
  //ky = second coordinate of the given point in the Brillouin zone
  //latticeComponent = index of the lattice vector along which the projection is to be performed
  //return value = projected momentum
  virtual double GetProjectedMomentum(int kx, int ky, int latticeComponent);
  
  //get the value of the embedding vectors
  //
  virtual void GetEmbedding(RealVector& embeddingX, RealVector& embeddingY);
  
  //set the value of the embedding vectors
  //
  virtual void SetEmbedding(RealVector embeddingX, RealVector embeddingY);
  
   //set the value of the embedding vectors from an external ascii file
  //
  //embeddingFileName = name of the ascii file that defines the embedding
  //
  bool SetEmbeddingFromAsciiFile(char* embeddingFileName);


  virtual int  EncodeSublatticeIndex(int posx, int posy,int & numXTranslations,int &numYTranslations, Complex &translationPhase);
  
   // compute the index in real space lattice starting from the cartesian coordinates
  //
  // i = cartesian coordinate in the x direction of the Bravais lattice
  // j = cartesian coordinate in the y direction of the Bravais lattice
  // p = reference on the first lattice index
  // q = reference on the second lattice index
  virtual void GetRealSpaceIndex (int i, int j, int& p, int& q);
  
 protected:

  // write an header that describes the tight binding model
  // 
  // output = reference on the output stream
  // return value  = reference on the output stream
  virtual ofstream& WriteHeader(ofstream& output);
  
   
  //computes all the values of the momentum projected and stores them in a double array
  //
  virtual void ComputeAllProjectedMomenta();
  
  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

  // build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
  //
  // nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
  // orbitalIndices = array that gives the orbital indices of the connected orbitals
  // spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
  // hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
  // return value = tight binding hamiltonian in real space 
  virtual HermitianMatrix BuildTightBindingHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes);

  // build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions but without assuming its hermiticiy
  //
  // nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
  // orbitalIndices = array that gives the orbital indices of the connected orbitals
  // spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
  // hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
  // return value = tight binding hamiltonian in real space 
  virtual ComplexMatrix BuildTightBindingNonHermitianHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes);

  // build the tight binding hamiltonian in recirpocal space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
  //
  // kx = momentum along the x direction (in 2pi /N_x unit) for which the hamiltonian in recirpocal space has to be computed
  // ky = momentum along the y direction (in 2pi /N_y unit) for which the hamiltonian in recirpocal space has to be computed
  // nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
  // orbitalIndices = array that gives the orbital indices of the connected orbitals
  // spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
  // hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
  // return value = tight binding hamiltonian in real space 
  virtual HermitianMatrix BuildTightBindingHamiltonianReciprocalSpace(int kx, int ky, int* nbrConnectedOrbitals, int** orbitalIndices, 
								      int** spatialIndices, Complex** hoppingAmplitudes);
};

// get the linearized momentum index
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// return value = linearized momentum index

inline int Abstract2DTightBindingModel::GetLinearizedMomentumIndex(int kx, int ky)
{
  return ((kx * this->NbrSiteY) + ky);
}

// get momentum value from a linearized momentum index
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// return value = inearized momentum index

inline void Abstract2DTightBindingModel::GetLinearizedMomentumIndex(int index, int& kx, int& ky)
{
  kx = index / this->NbrSiteY;
  ky = index % this->NbrSiteY;
}

// get the linearized momentum index, without assuming k to be in the first BZ
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// return value = linearized momentum index

inline int Abstract2DTightBindingModel::GetLinearizedMomentumIndexSafe(int kx, int ky)
{
  while (kx < 0)
      kx += this->NbrSiteX;
  kx %= this->NbrSiteX;
  while (ky < 0)
      ky += this->NbrSiteY;
  ky %= this->NbrSiteY;
  return ((kx * this->NbrSiteY) + ky);
}

// get momentum value from a linearized momentum index, without assuming k to be in the first BZ
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// return value = inearized momentum index

inline void Abstract2DTightBindingModel::GetLinearizedMomentumIndexSafe(int index, int& kx, int& ky)
{
  int n = this->NbrSiteX * this->NbrSiteY;
  while (index < 0)
      index += n;
  index %= n;
  kx = index / this->NbrSiteY;
  ky = index % this->NbrSiteY;
}

// get the index of the real space tight binding model from the real space coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int Abstract2DTightBindingModel::GetRealSpaceTightBindingLinearizedIndex(int x, int y, int orbitalIndex)
{
  return (orbitalIndex + ((y  + x * this->NbrSiteY) * this->NbrBands)); 
}

// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int Abstract2DTightBindingModel::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int orbitalIndex)
{
  orbitalIndex %= this->NbrBands;
  if (orbitalIndex < 0)
    orbitalIndex +=  this->NbrBands;
  x %= this->NbrSiteX;
  if (x < 0)
    x +=  this->NbrSiteX;
  y %= this->NbrSiteY;
  if (y < 0)
    y +=  this->NbrSiteY;
  return this->GetRealSpaceTightBindingLinearizedIndex(x, y, orbitalIndex); 
}

// get the real space coordinates from the index of the real space tight binding model
//
// index = linearized index of the real space tight binding model
// x = reference on the x coordinate of the unit cell
// y = reference on the y coordinate of the unit cell
// orbitalIndex = reference on the index of the orbital / site within the unit cell

inline void Abstract2DTightBindingModel::GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y, int& orbitalIndex)
{
  orbitalIndex = index % this->NbrBands;
  index /= this->NbrBands;
  y = index % this->NbrSiteY;
  x = index / this->NbrSiteY;
}

// get the angle between the two primitive lattice vectors
//
// return value = angle between the two primitive lattice vectors

inline double Abstract2DTightBindingModel::GetTwistAngle()
{
  return this->TwistAngle;
}

// get the number of sites in the y direction
//
// return value = number of sites in the y direction

inline int Abstract2DTightBindingModel::GetNbrSiteY()
{
  return this->NbrSiteY;
}

// get the total curvature over the plaquette specified by the lower-left corner. C > 0 corresponds to Bz < 0 magnetic field (holomorphic z = x + iy)
//
// kx = momentum along x
// ky = momentum along y
// band = band index
// return value = Arg(<k|k+x><k+x|k+x+y><k+x+y|k+y><k+y|k>) / (2 Pi), range = (-0.5, 0.5]

inline double Abstract2DTightBindingModel::GetCurvature(int kx, int ky, int band)
{
    if (this->Curvature == NULL)
        this->ComputeCurvature();
    return this->Curvature[band][ky][kx];
}

// get the Chern number of a specified band
//
// band = band index
// return = Chern number

inline int Abstract2DTightBindingModel::GetChernNumber(int band)
{
    if (this->Chern == NULL)
        this->ComputeCurvature();
    return this->Chern[band];
}

// get the LLL twist angle along x for Bloch construction gauge fixing
//
// band = band index
// return = gamma_x

inline double Abstract2DTightBindingModel::GetLLLGammaX(int band)
{
    if (this->LLLGammaX == NULL)
        this->ComputeCurvature();
    return this->LLLGammaX[band];
}

// get the LLL twist angle along y for Bloch construction gauge fixing
//
// band = band index
// return = gamma_y
inline double Abstract2DTightBindingModel::GetLLLGammaY(int band)
{
    if (this->LLLGammaY == NULL)
        this->ComputeCurvature();
    return this->LLLGammaY[band];
}

// compute the curvature sum over the plaquette specified by the lower-left corner. C > 0 corresponds to Bz < 0 magnetic field (holomorphic z = x + iy)
//
// kx = momentum along x
// ky = momentum along y
// band = band index
// return value = Arg(<k|k+x><k+x|k+x+y><k+x+y|k+y><k+y|k>) / (2 Pi), range = (-0.5, 0.5]

inline double Abstract2DTightBindingModel::ComputeCurvatureSinglePlaquette(int kx, int ky, int band)
{
    Complex W = this->GetAbelianConnection(kx, ky, 1, 0, band) * this->GetAbelianConnection(kx + 1, ky, 0, 1, band);
    W *= this->GetAbelianConnection(kx + 1, ky + 1, -1, 0, band) * this->GetAbelianConnection(kx, ky + 1, 0, -1, band);
    return Arg(W) / (2 * M_PI);
}

// Computes value of projected momentum along the lattice directions
//
// kx = first coordinate of the given point in the Brillouin zone
// ky = second coordinate of the given point in the Brillouin zone
// latticeComponent = index of the lattice vector along which the projection is to be performed
// return value = projected momentum
   
  
inline double Abstract2DTightBindingModel::GetProjectedMomentum(int kx, int ky, int latticeComponent)
{
  return this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky)][latticeComponent];
}

 //get the value of the embedding vectors
 //
 inline void Abstract2DTightBindingModel::GetEmbedding(RealVector& embeddingX, RealVector& embeddingY)
 {
   embeddingX = this->EmbeddingX;
   embeddingY = this->EmbeddingY;
 }
 
 //get the value of the embedding vectors
 //
 inline void Abstract2DTightBindingModel::SetEmbedding(RealVector embeddingX, RealVector embeddingY)
 {
   this->EmbeddingX = embeddingX;
   this->EmbeddingY = embeddingY;
 }

 
  // compute the index in real space lattice starting from the cartesian coordinates
  //
  // i = cartesian coordinate in the x direction of the Bravais lattice
  // j = cartesian coordinate in the y direction of the Bravais lattice
  // p = reference on the first lattice index
  // q = reference on the second lattice index
  inline void Abstract2DTightBindingModel::GetRealSpaceIndex (int i, int j, int& p, int& q)
  {
    p = i - this->OffsetReal * j;
    q = j;
  }

// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// numXTranslations = number of translation in the x direction to get back to the unit cell 
// numXTranslations = number of translation in the y direction to get back to the unit cell
//
inline int  Abstract2DTightBindingModel::EncodeSublatticeIndex(int posx, int posy,int & numXTranslations,int &numYTranslations, Complex &translationPhase)
{
  std::cout <<"using dummy Abstract2DTightBindingModel::EncodeSublatticeIndex(int posx, int posy,int & numXTranslations,int &numYTranslations, Complex &translationPhase)"<<std::endl;
return -1;
}



#endif
