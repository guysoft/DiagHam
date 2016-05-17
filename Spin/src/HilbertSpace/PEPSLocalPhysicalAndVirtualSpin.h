////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-20016 Nicolas Regnault                 //
//                         class author: Antoine Sterdyniak                   //
//                                                                            //
//                                                                            //
//                 class of local spins around a single PEPS tensor           //
//                                                                            //
//                        last modification : 09/05/2016                      //
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


#ifndef PEPSLOCALPHYSICALANDVIRTUALSPIN_H
#define PEPSLOCALPHYSICALANDVIRTUALSPIN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"


class PEPSLocalPhysicalAndVirtualSpin : public AbstractSpinChain
{


 protected:

  int PhysicalSpinValue;
  int NbrSpinStatesForPhysicalSpin;

  int NbrVirtualSpinRepresentation;
  int NbrSpinStatesForVirtualSpins;

  int Sz;

  int * TableOfVirtualSpinRepresentations;
  int * TableOfSzValuesForVirtualSpin;
  int * TableOfSValuesForVirtualSpin;
  int * TableOfSzValuesForPhysicalSpin;
  int * PowerOfNbrSpinStatesForVirtualSpins;

  unsigned long * StateDescription;
  
  
 public:
  
  
  // virtual destructor
  //
  PEPSLocalPhysicalAndVirtualSpin ();
  
  // virtual destructor
  //
  PEPSLocalPhysicalAndVirtualSpin (int physicalSpinValue, int nbrVirtualSpinRepresentation, int * tableOfVirtualSpinRepresentations,int sz);
  
  
  //  virtual destructor
  //
  ~PEPSLocalPhysicalAndVirtualSpin ();
  
  
  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);
  
  virtual AbstractHilbertSpace* Clone();
  void PrintConversionTable();

  // get the value of the spin (i.e. S) at a given site
  // 
  // site = site index
  // return value = twice the spin
  virtual int GetLocalSpin(int site,int state);

  // return eigenvalue of Sz_i Sz_j associated to a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  virtual double SziSzj (int i, int j, int state);


  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state 
  virtual int SmiSpj (int i, int j, int state, double& coefficient);

  void ApplyDMinusXReflexion (ComplexVector & initialState, ComplexVector &  FinalState);    
  void ApplyDXReflexion (ComplexVector & initialState, ComplexVector &  FinalState);
  void ApplyVReflexion (ComplexVector & initialState, ComplexVector &  FinalState);
  void ApplyHReflexion (ComplexVector & initialState, ComplexVector &  FinalState);
  void ApplyC4Rotation (ComplexVector & initialState, ComplexVector &  FinalState);

 protected:
  
  // find state index
  //
  // state = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long state);
    
  int GeneratesStates(int nbrVirtualSpinToBeAttributed, int szToBeRealized, int beginningOfStateRepresentation, int pos);
  
  int EvaluateHilbertSpaceDimension(int nbrVirtualSpinToBeAttributed, int szToBeRealized);
  

};

#endif


