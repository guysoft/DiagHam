////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of manager for FQHE MPS matrices                  //
//                                                                            //
//                        last modification : 31/10/2012                      //
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


#ifndef FQHEMPSMATRIXMANAGER_H
#define FQHEMPSMATRIXMANAGER_H


#include "config.h"


class OptionManager;
class AbstractFQHEMPSMatrix;
class AbstractArchitecture;


class FQHEMPSMatrixManager
{

 protected:

  // pointer to the option manager
  OptionManager* Options;

  // indicates that the MPS matrix will be used to compute a transfer matrix
  bool EMatrixFlag;

 public:
  
  // default constructor 
  //
  // eMatrixFlag = indicates that the MPS matrix will be used to compute a transfer matrix
  FQHEMPSMatrixManager(bool eMatrixFlag = false);

  // destructor
  //
  ~FQHEMPSMatrixManager();
  
  // add an option group containing all options related to the MPS matrix construction 
  //
  // manager = pointer to the option manager
  // comment = additional comment that is displayed in the behind each option group
  void AddOptionGroup(OptionManager* manager, const char* comment = 0);

  // get the MPS matrice class defined by the running options
  //
  // nbrFluxQuanta = number of flux quanta
  // architecture = architecture to use for precalculation
  // return value = pointer to the MPS matrice class 
  AbstractFQHEMPSMatrix* GetMPSMatrices(int nbrFluxQuanta = 0, AbstractArchitecture* architecture = 0);

  // get the MPS matrice class defined by the running options for the right part of the transfer matrix
  //
  // nbrFluxQuanta = number of flux quanta
  // architecture = architecture to use for precalculation
  // return value = pointer to the MPS matrice class 
  AbstractFQHEMPSMatrix* GetRightMPSMatrices(int nbrFluxQuanta = 0, AbstractArchitecture* architecture = 0);

  // get the MPS matrice class defined by the running options for the right part of the transfer matrix
  //
  // nbrFluxQuanta = number of flux quanta
  // architecture = architecture to use for precalculation
  // return value = pointer to the MPS matrice class 
  AbstractFQHEMPSMatrix* GetLeftMPSMatrices(int nbrFluxQuanta = 0, AbstractArchitecture* architecture = 0);

  // get the cylinder perimeter (in magnetic length unit) if the cylinder geometry if used
  //
  // nbrFluxQuanta = number of flux quanta
  // return value = cylinder perimeter (negative if another geometry is used)
  double GetCylinderPerimeter(int nbrFluxQuanta = 0);

 protected:

  // get the MPS matrice class defined by the running options and additional external parameters
  //
  // quasiholeSectorFlag = use the quasihole sector
  // topologicalSectorValue = select a specific topological sector when grouping B matrices
  // importBMatrices = import the B matrices from a file
  // nbrFluxQuanta = number of flux quanta
  // architecture = architecture to use for precalculation
  // return value = pointer to the MPS matrice class  
  AbstractFQHEMPSMatrix* GetMPSMatrices(bool quasiholeSectorFlag, int topologicalSectorValue, char* importBMatrices, int nbrFluxQuanta, AbstractArchitecture* architecture);

};

#endif
