////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2003 Duc-Phuong Nguyen                    //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 24/11/2003                        //
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


#ifndef PERIODICQUANTUMDOTS3DHAMILTONIAN_H
#define PERIODICQUANTUMDOTS3DHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"
#include "Tools/QuantumDot/Potential/ThreeDPotential.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class PeriodicQuantumDots3DHamiltonian : public AbstractHamiltonian
{

 protected:

  // Hilbert space associated to the system
  Periodic3DOneParticle* Space;

  // wave function basis dimension in the x direction
  int NbrStateX;

  int LowerImpulsionX;

  // wave function basis dimension in the y direction
  int NbrStateY;

  int LowerImpulsionY;

  // wave function basis dimension in the z direction
  int NbrStateZ;

  int LowerImpulsionZ;

  // number of cells in the x direction
  int NbrCellX;
  // number of cells in the y direction
  int NbrCellY;
  // number of cells in the z direction
  int NbrCellZ;

  // system dimension in the x direction (in Angstrom unit)
  double XSize;
  // system dimension in the y direction (in Angstrom unit)
  double YSize;
  // system dimension in the z direction (in Angstrom unit)
  double ZSize;

  // effective mass in the x direction (in electron mass unit)
  double Mux;
  // effective mass in the y direction (in electron mass unit)
  double Muy;
  // effective mass in the z direction (in electron mass unit)
  double Muz;

  // cache for hamiltonian diagonal elements
  double* KineticElements;

  // flag to indicate how many dimensions are precalculated for hamiltonian interaction elements
  int NbrPrecalculatedDimension;

  // tridimensionnal array containing all cell interaction factors
  double*** InteractionFactors;
  //ThreeDPotential* InteractionFactors;

  // wave function overlaps on a cell in a the x direction (with symmetric access type for the first two indices)
  double** RealWaveFunctionOverlapX;
  double** ImaginaryWaveFunctionOverlapX;
  // wave function overlaps on a cell in a the y direction (with symmetric access type for the first two indices)
  double** RealWaveFunctionOverlapY;
  double** ImaginaryWaveFunctionOverlapY;
  // wave function overlaps on a cell in a the z direction (with symmetric access type for the first two indices)
  double** RealWaveFunctionOverlapZ;
  double** ImaginaryWaveFunctionOverlapZ;  

  // tridimensionnal array (with symmetric access type i > j for the two first indices) to store partial calculation to construct hamiltonian
  // elements (integration over two dimensions, first and second indices are of the form n1 * dim + m1)
  double*** RealPrecalculatedHamiltonian;
  double*** ImaginaryPrecalculatedHamiltonian;

 public:
  
  // constructor from default data
  //
  // space = Hilbert space associated to the system
  // xSize = system dimension in the x direction (in Angstrom unit)
  // ySize = system dimension in the y direction (in Angstrom unit)
  // zSize = system dimension in the z direction (in Angstrom unit)
  // mux = effective mass in the x direction (in electron mass unit)
  // muy = effective mass in the y direction (in electron mass unit)
  // muz = effective mass in the z direction (in electron mass unit)
  // nbrCellX = number of cells in the x direction
  // nbrCellY = number of cells in the y direction
  // nbrCellZ = number of cells in the z direction
  // potentialInput = pointer to a 3D potential  
  PeriodicQuantumDots3DHamiltonian(Periodic3DOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDPotential* PotentialInput);

  // copy constructor (without duplicating datas)
  //
  // hamiltonian = reference on hamiltonian to copy  
  PeriodicQuantumDots3DHamiltonian(const PeriodicQuantumDots3DHamiltonian& hamiltonian);

  // destructor
  //  
  ~PeriodicQuantumDots3DHamiltonian();
  
  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian  
  AbstractHamiltonian* PeriodicQuantumDots3DHamiltonian::Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();

  // shift Hamiltonian from a given energy
  //
  // shift = shift value

  void ShiftHamiltonian (double shift);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  
  Complex MatrixElement (RealVector& V1, RealVector& V2);


  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element

  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);


  // multiply a vector by the current hamiltonian for a given range of idinces 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
				  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);
  
  // evaluate the all interaction factors
  // 
  void EvaluateInteractionFactors();

  // evaluate wave function overlap in a given direction
  //
  // nbrStep = number of subdivisions in the choosen direction
  // nbrState = number of states in the restrained Hilbert space in the choosen direction
  // realArray = reference to a 2D array to stock the real values of overlap
  // imaginaryArray = reference to a 2D array to stock the imaginary values of overlap
  bool EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray);
};

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

inline AbstractHilbertSpace* PeriodicQuantumDots3DHamiltonian::GetHilbertSpace ()
{
  return this->Space;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

inline int PeriodicQuantumDots3DHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Space->GetHilbertSpaceDimension ();
}

#endif
