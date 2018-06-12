////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                         Class author : Cecile Repellin                     //
//                                                                            //
//         class of tight binding model for twisted bilayer graphene          //
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


#ifndef TIGHTBINDINGMODELBILAYERTWISTEDBILAYER_H
#define TIGHTBINDINGMODELBILAYERTWISTEDBILAYER_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelBilayerTwistBilayer : public Abstract2DTightBindingModel
{

 protected:

  // hopping amplitude between neareast neighbor sites in the same layer
  double NNHoppingIntra;
  // hopping amplitude between layers
  double NNHoppingInter;
  // hopping amplitude between next neareast neighbor sites
  double GammaOne;
  // spin orbit coupling to neareast neighbor sites
  double TrigonalWarping;
  // twisting angle in units of 
  double TwistingAngle;
  // voltage
  double UVoltage;
  
  // index of the band of where interaction has to be projected
  int BandIndex;
  
  //coordinates of the reciprocal lattice vectors
  double Gx1;
  double Gy1;
  double Gx2;
  double Gy2;
  
  double LatticeConstant;
  
  
  int NbrPointSuperlatticeX;
  int NbrPointSuperlatticeY;
  
  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;
  // nearest neighbor density-density potential strength

 public:

  // default constructor
  //
  TightBindingModelBilayerTwistBilayer();
     
     
  // constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t = tunneling amplitude within one layer
  // t3 = trigonal warping
  // g1 = tunneling amplitude within one layer
  // tM = tunneling amplitude between layers
  // uVoltage = amplitude of voltage
  // theta = twisting angle     
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelBilayerTwistBilayer(int nbrSiteX, int nbrSiteY, int nbrPointsX, int nbrPointsY, double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices);

  // destructor
  //
  ~TightBindingModelBilayerTwistBilayer();

  // compute the Bloch hamiltonian at a point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // ky = momentum along the x axis
  // return value = Bloch hamiltonian
  virtual HermitianMatrix ComputeBlochHamiltonian(double kx, double ky);

  // get the high symmetry points 
  //
  // pointNames = name of each high symmetry point
  // pointCoordinates = coordinates in the first Brillouin zone of the high symmetry points
  // return value = number of high symmetry points
  virtual int GetHighSymmetryPoints(char**& pointNames, double**& pointCoordinates);

  // compute the distance between two points in the first Brillouin zone, changing the coordinates the second one by a reciprocal lattice vector if needed
  //
  // kx1 = momentum of the first point along the x axis
  // ky1 = momentum of the first point along the y axis
  // kx2 = reference on the momentum of the second point along the x axis
  // ky2 = reference on the momentum of the second point along the y axis
  // return value = distance between the two points
  virtual double GetDistanceReciprocalSpace(double kx1, double ky1, double& kx2, double& ky2);
  
  // compute the distance between two points in the same unit cell, changing the coordinates the second one by a lattice vector if needed
  //
  // x1 = x coordinate of the first point
  // y1 = y coordinate of the first point
  // x2 = reference to the x coordinate of the second point
  // y2 = referenceto the y coordinate of the second point
  // return value = distance between the two points
  virtual double GetDistanceRealSpace(double x1, double y1, double& x2, double& y2);  
  
  // get the linearized momentum index
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // return value = linearized momentum index
  virtual int GetLinearizedMomentumIndexSuperlattice(int kx, int ky);
  
  // get momentum value from a linearized momentum index
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // return value = inearized momentum index
  virtual void GetLinearizedMomentumIndexSuperlattice(int index, int& kx, int& ky);
  
  // get the linearized momentum index, without assuming k to be in the first BZ
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // return value = linearized momentum index
  virtual int GetLinearizedMomentumIndexSuperlatticeSafe(int kx, int ky);
  
  // get momentum value from a linearized momentum index, without assuming k to be in the first BZ
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // return value = inearized momentum index
  virtual void GetLinearizedMomentumIndexSuperlatticeSafe(int index, int& kx, int& ky);
  
  
  // compute the form factor for the density operator 
  // 
  // kx = momentum along x of annihilation operator
  // ky = momentum along y of creation operator
  // qx = momentum transfer along x direction
  // qy = momentum transfer along y direction
  // valleyIndex = valley index of density operator
  virtual Complex ComputeDensityFormFactor(int kx, int ky, int qx, int qy, int valleyIndex);
  
  // evaluate the norm of a momentum space vector
  //
  // kx = component of momentum along first Bravais vector
  // ky = component of momentum along second Bravais vector
  // return value = norm of vector
  virtual double EvaluateNormQ(int kx, int ky);

 protected :

  // core part that computes the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


// get the linearized momentum index
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// return value = linearized momentum index

inline int TightBindingModelBilayerTwistBilayer::GetLinearizedMomentumIndexSuperlattice(int kx, int ky)
{
  return ((kx * (this->NbrPointSuperlatticeY)) + ky);
}

// get momentum value from a linearized momentum index
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// return value = inearized momentum index

inline void TightBindingModelBilayerTwistBilayer::GetLinearizedMomentumIndexSuperlattice(int index, int& kx, int& ky)
{
  kx = index / (this->NbrPointSuperlatticeY);
  ky = index % (this->NbrPointSuperlatticeY);
}

// get the linearized momentum index, without assuming k to be in the first BZ
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// return value = linearized momentum index

inline int TightBindingModelBilayerTwistBilayer::GetLinearizedMomentumIndexSuperlatticeSafe(int kx, int ky)
{
  while (kx < 0)
      kx += (this->NbrPointSuperlatticeX);
  kx %= (this->NbrPointSuperlatticeX);
  while (ky < 0)
      ky += (this->NbrPointSuperlatticeY);
  ky %= (this->NbrPointSuperlatticeY);
  return ((kx * (this->NbrPointSuperlatticeY)) + ky);
}

// get momentum value from a linearized momentum index, without assuming k to be in the first BZ
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// return value = inearized momentum index

inline void TightBindingModelBilayerTwistBilayer::GetLinearizedMomentumIndexSuperlatticeSafe(int index, int& kx, int& ky)
{
  int n = (this->NbrPointSuperlatticeX) * (this->NbrPointSuperlatticeY);
  while (index < 0)
      index += n;
  index %= n;
  kx = index / (this->NbrPointSuperlatticeY);
  ky = index % (this->NbrPointSuperlatticeY);
}

// evaluate the norm of a momentum space vector
//
// kx = component of momentum along first Bravais vector
// ky = component of momentum along second Bravais vector
// return value = norm of vector
inline double TightBindingModelBilayerTwistBilayer::EvaluateNormQ(int kx, int ky)
{
    double KX = (this->Gx1 * this->KxFactor * kx + this->Gx2 * this->KyFactor * ky) * this->LatticeConstant;
    double KY = (this->Gy1 * this->KxFactor * kx + this->Gy2 * this->KyFactor * ky) * this->LatticeConstant;
    return sqrt(KX*KX + KY*KY);
}

#endif
