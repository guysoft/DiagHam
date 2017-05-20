////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of MPS matrix for the PH-Pfaffian state               //
//                                                                            //
//                        last modification : 19/05/2017                      //
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


#ifndef FQHEMPSPHPFAFFIANMATRIX_H
#define FQHEMPSPHPFAFFIANMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Vector/LongRationalVector.h"


class LongRationalMatrix;
class RealMatrix;
class RealSymmetricMatrix;
class BosonOnDiskShort;
class AbstractArchitecture;

class FQHEMPSPHPfaffianMatrix : public FQHEMPSClustered2RMatrix
{

  friend class FQHEMPSEvaluateCFTOperation;

  // shift to apply to the p level to get a positive index
  int PLevelShift;

 public:
  
  // default constructor 
  //
  FQHEMPSPHPfaffianMatrix();

  // constructor 
  //
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cftDirectory = path to the directory where all the pure CFT matrices are stored
  // bosonicVersion = use a version of the code that is compatible with bosonic wave functions
  // useRational = use arbitrary precision numbers for all the CFT calculations
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  // architecture = architecture to use for precalculation
  FQHEMPSPHPfaffianMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, const char* cftDirectory, bool bosonicVersion = false, bool useRational = true,
			  bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
			  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0,
			  AbstractArchitecture* architecture = 0);


  // constructor from stored B matrices
  //
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  FQHEMPSPHPfaffianMatrix(int laughlinIndex, int pLevel, const char* fileName, bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
			  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0);

  // destructor
  //
  ~FQHEMPSPHPfaffianMatrix();
  
  // get the number of particles that fit the root configuration once the number of flux quanta is fixed
  // 
  // nbrFluxQuanta = number of flux quanta
  // padding = assume that the state has the extra padding
  // return value = number of partciles
  virtual int GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding = false);

  // create the B matrices for the laughlin state
  //
  // cftDirectory = an optional path to the directory where all the CFT matrices are stored
  // architecture = architecture to use for precalculation
  virtual void CreateBMatrices (const char* cftDirectory, AbstractArchitecture* architecture);

 protected:


};

  


#endif
