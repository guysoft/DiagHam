////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
// class of hilbert space of one particle in 3d with planar rotation symmetry //
//                                                                            //
//                        last modification : 11/08/2004                      //
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


#ifndef ABSTRACTPLANARROTATIONSYMMETRYTHREEDONEPARTICLE_H
#define ABSTRACTPLANARROTATIONSYMMETRYTHREEDONEPARTICLE_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace/QuantumDotHilbertSpace/AbstractOneDOneParticle.h"


class AbstractPlanarRotationSymmetryThreeDOneParticle : public AbstractHilbertSpace
{

 protected:
  
  // quantum number of the Z kinetic momentum component
  int Lz;
  // wave function basis dimension of the radial coordinate
  AbstractOneDOneParticle* StateR;
  // wave function basis dimension in the Z direction
  AbstractOneDOneParticle* StateZ;

 public:

  // default constructor
  //
  AbstractPlanarRotationSymmetryThreeDOneParticle ();

  // constructor
  //
  // lz = quantum number of the Z kinetic momentum component
  // nbrStateR = wave function basis dimension of the radial coordinate
  // nbrStateZ = wave function basis dimension in the z direction
  AbstractPlanarRotationSymmetryThreeDOneParticle (int lz, int nbrStateR, int nbrStateZ);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  AbstractPlanarRotationSymmetryThreeDOneParticle (const AbstractPlanarRotationSymmetryThreeDOneParticle& space);

  // destructor
  //
  ~AbstractPlanarRotationSymmetryThreeDOneParticle ();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  AbstractPlanarRotationSymmetryThreeDOneParticle& operator = (const AbstractPlanarRotationSymmetryThreeDOneParticle& space);

  // get the quantum number of the Z kinetic momentum component
  //
  // return = value of quantum number of the Z kinetic momentum component
  int GetLz ();

  // get wave function basis dimension in the plane
  //
  // return value = wave function basis dimension in the plane
  int GetNbrStateR ();

  // get wave function basis dimension in the z direction
  //
  // return value = wave function basis dimension in the z direction
  int GetNbrStateZ ();

  // get the wave function basis in the plane
  //
  // return = pointer to 1D one particle basis
  virtual AbstractOneDOneParticle* GetStateR ();

  // get the wave function basis in Z direction
  //
  // return = pointer to 1D one particle basis
  virtual AbstractOneDOneParticle* GetStateZ ();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
					 SubspaceSpaceConverter& converter);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

};

// get the quantum number of the Z kinetic momentum component
//
// return = value of quantum number of the Z kinetic momentum component

inline int AbstractPlanarRotationSymmetryThreeDOneParticle::GetLz ()
{
  return this->Lz;
}

// get wave function basis dimension in the plane
//
// return value = wave function basis dimension in the plane

inline int AbstractPlanarRotationSymmetryThreeDOneParticle::GetNbrStateR ()
{
  return this->StateR->GetNbrState ();
}

// get wave function basis dimension in the z direction
//
// return value = wave function basis dimension in the z direction

inline int AbstractPlanarRotationSymmetryThreeDOneParticle::GetNbrStateZ ()
{
  return this->StateZ->GetNbrState ();
}

// get the wave function basis in the plane
//
// return = pointer to 1D one particle basis

inline AbstractOneDOneParticle* AbstractPlanarRotationSymmetryThreeDOneParticle::GetStateR ()
{
  AbstractOneDOneParticle* stateR = (AbstractOneDOneParticle*) this->StateR->Clone ();
  return stateR;
}

// get the wave function basis in Z direction
//
// return = pointer to 1D one particle basis

inline AbstractOneDOneParticle* AbstractPlanarRotationSymmetryThreeDOneParticle::GetStateZ ()
{
  AbstractOneDOneParticle* stateZ = (AbstractOneDOneParticle*) this->StateZ->Clone ();
  return stateZ;
}

#endif
