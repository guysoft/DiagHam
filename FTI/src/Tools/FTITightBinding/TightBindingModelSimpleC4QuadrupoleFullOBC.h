////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the simple C4 quadrupole       //
//                  insulator with full open boundary conditions              //
//                                                                            //
//                        last modification : 10/03/2018                      //
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


#ifndef TIGHTBINDINGMODELSIMPLEC4QUADRUPOLEFULLOBC_H
#define TIGHTBINDINGMODELSIMPLEC4QUADRUPOLEFULLOBC_H


#include "config.h"
#include "Tools/FTITightBinding/AbstractTightBindingModel.h"


class TightBindingModelSimpleC4QuadrupoleFullOBC : public AbstractTightBindingModel
{

 protected:

   // number of sites in the x direction
  int NbrSiteX;
   // number of sites in the y direction
  int NbrSiteY;

  // hoping amplitude between neareast neighbor sites within the unit cell
  double NNHopping1;
  // phase (for the the hoping between neareast neighbor sites within the unit cell
  double NNHoppingPhase1;

  // hoping amplitude between neareast neighbor sites between unit cells
  double NNHopping2;
  // phase for the the hoping between neareast neighbor sites between unit cells
  double NNHoppingPhase2;
  
  // linearized coordiantes of the confining potential
  int* ConfiningPotentialCoordinates;
  // amplitudes of the confining potential on each sites
  double* ConfiningPotentialAmplitudes;
  // number of sites where there the confining potential has a non-zero amplitude
  int NbrConfiningPotentials;


 public:

  // default constructor
  //
  TightBindingModelSimpleC4QuadrupoleFullOBC();

  // basic constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction 
  // t1 = hoping amplitude between neareast neighbor sites within the unit cell
  // phi1 = phase (in pi units) for the the hoping between neareast neighbor sites within the unit cell
  // t2 = hoping amplitude between neareast neighbor sites between unit cells
  // phi2 = phase (in pi units) for the the hoping between neareast neighbor sites between unit cells
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelSimpleC4QuadrupoleFullOBC(int nbrSiteX, int nbrSiteY, double t1, double phi1, double t2, double phi2, 
					     AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // constructor with an additional confining potential
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction 
  // t1 = hoping amplitude between neareast neighbor sites within the unit cell
  // phi1 = phase (in pi units) for the the hoping between neareast neighbor sites within the unit cell
  // t2 = hoping amplitude between neareast neighbor sites between unit cells
  // phi2 = phase (in pi units) for the the hoping between neareast neighbor sites between unit cells
  // confiningPotentialXCoordinates = x coordiantes of the confining potential
  // confiningPotentialYCoordinates = y coordiantes of the confining potential
  // confiningPotentialAmplitudes = amplitudes of the confining potential on each sites
  // nbrConfiningPotentials = number of sites where there the confining potential has a non-zero amplitude
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelSimpleC4QuadrupoleFullOBC(int nbrSiteX, int nbrSiteY, double t1, double phi1, double t2, double phi2,
					      int* confiningPotentialXCoordinates, int* confiningPotentialYCoordinates, 
					      double* confiningPotentialAmplitudes, int nbrConfiningPotentials,
					      AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // destructor
  //
  ~TightBindingModelSimpleC4QuadrupoleFullOBC();

  // get the index of the real space tight binding model from the real space coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndex(int x, int y);
  
  // get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // return value = linearized index  (negative if non valid)
  virtual int GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y);

  // get the real space coordinates from the index of the real space tight binding model
  //
  // index = linearized index of the real space tight binding model
  // x = reference on the x coordinate of the unit cell
  // y = reference on the y coordinate of the unit cell
  virtual void GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y);

  // get the size (length / area / volume ... ) of the unit cell
  //
  // return value =  size
  virtual double GetUnitCellSize();

  // get the energy at a given momentum of the band structure
  //
  // bandIndex = index of the band to consider
  // momentumIndex = linearized momentum
  // return value = energy
  virtual double GetEnergy(int bandIndex, int momentumIndex);

  // ask if the one body transformation matrices are available
  //
  // return value = true if the one body transformation matrices are available
  virtual bool HaveOneBodyBasis();

  // get  the one body transformation matrix corresponding to a given momentum of the band structure
  //
  // momentumIndex = linearized momentum
  // return value = reference on the one body transformation matrix
  virtual ComplexMatrix& GetOneBodyMatrix(int momentumIndex);

  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

 protected :

  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

};

// get the index of the real space tight binding model from the real space coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int TightBindingModelSimpleC4QuadrupoleFullOBC::GetRealSpaceTightBindingLinearizedIndex(int x, int y)
{
  return (y  + (x * this->NbrSiteY)); 
}

// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// return value = linearized index (negative if non valid)

inline int TightBindingModelSimpleC4QuadrupoleFullOBC::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y)
{
  if ((x < 0) || (x >= this->NbrSiteX))
    {
      return -1;
    }
  if ((y < 0) || (y >= this->NbrSiteY))
    {
      return -1;
    }
  return this->GetRealSpaceTightBindingLinearizedIndex(x, y); 
}

// get the real space coordinates from the index of the real space tight binding model
//
// index = linearized index of the real space tight binding model
// x = reference on the x coordinate of the unit cell
// y = reference on the y coordinate of the unit cell

inline void TightBindingModelSimpleC4QuadrupoleFullOBC::GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y)
{
  y = index % this->NbrSiteY;
  x = index / this->NbrSiteY;
}

// get the size (length / area / volume ... ) of the unit cell
//
// return value =  size

inline double TightBindingModelSimpleC4QuadrupoleFullOBC::GetUnitCellSize()
{
  return 1.0;
}


// get the energy at a given momentum of the band structure
//
// bandIndex = index of the band to consider
// momentumIndex = linearized momentum
// return value = energy

inline double TightBindingModelSimpleC4QuadrupoleFullOBC::GetEnergy(int bandIndex, int momentumIndex)
{
  return this->EnergyBandStructure[bandIndex][0];
}

// ask if the one body transformation matrices are available
//
// return value = true if the one body transformation matrices are available

inline bool TightBindingModelSimpleC4QuadrupoleFullOBC::HaveOneBodyBasis()
{
  if (this->OneBodyBasis != 0)
    {
      return true;
    }
  else
    {
      return false;
    }
}

// get  the one body transformation matrix corresponding to a given momentum of the band structure
//
// momentumIndex = linearized momentum
// return value = reference on the one body transformation matrix

inline ComplexMatrix& TightBindingModelSimpleC4QuadrupoleFullOBC::GetOneBodyMatrix(int momentumIndex)
{
  return this->OneBodyBasis[0];
}


#endif
