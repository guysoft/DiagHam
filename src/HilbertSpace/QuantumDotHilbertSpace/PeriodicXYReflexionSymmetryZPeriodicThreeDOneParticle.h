////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//class of hilbert space of one 3d periodic box particle with reflexion symmetry in the plane//
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


#ifndef PERIODICXYREFLEXIONSYMMETRYZPERIODICTHREEDONEPARTICLE_H
#define PERIODICXYREFLEXIONSYMMETRYZPERIODICTHREEDONEPARTICLE_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PeriodicOneDOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PeriodicReflexionSymmetryOneDOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/AbstractThreeDOneParticle.h"


class PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle : public AbstractThreeDOneParticle
{

 protected:

 public:

  // default constructor
  //
  PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle ();

  // constructor
  //
  // nbrStateX = wave function basis dimension in the x direction without symmetry redundancy
  // evenX = true if the wave functions in X is even, else odd
  // nbrStateY = wave function basis dimension in the y direction without symmetry redundancy
  // evenY = true if the wave functions in Y is even, else odd
  // nbrStateZ = wave function basis dimension in the z direction
  // lowZ = lower impulsion in Z direction (in unit of 2 * Pi / Lz)
  PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle (int nbrStateX, bool evenX, int nbrStateY, bool evenY, int nbrStateZ, int lowZ);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle (const PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle& space);

  // destructor
  //
  virtual ~PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle& operator = (const PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle& space);

  // get lower impulsion in X
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsionX ();

  // get lower impulsion in Y
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsionY ();

  // get lower impulsion in Z
  //
  // return value = lower value in unit of 2 * Pi / L
  virtual int GetLowerImpulsionZ ();

  // get the wave function basis in X direction
  //
  // return = pointer to 1D one particle basis
  virtual AbstractOneDOneParticle* GetStateX ();

  // get the wave function basis in Y direction
  //
  // return = pointer to 1D one particle basis
  virtual AbstractOneDOneParticle* GetStateY ();

  // get the wave function basis in Z direction
  //
  // return = pointer to 1D one particle basis
  virtual AbstractOneDOneParticle* GetStateZ ();
};

// get lower impulsion in X
//
// return value = lower value in unit of 2 * Pi / L

inline int PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle::GetLowerImpulsionX ()
{
  return ((PeriodicReflexionSymmetryOneDOneParticle*) this->StateX)->GetLowerImpulsion ();
}

// get lower impulsion in Y
//
// return value = lower value in unit of 2 * Pi / L

inline int PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle::GetLowerImpulsionY ()
{
  return ((PeriodicReflexionSymmetryOneDOneParticle*) this->StateY)->GetLowerImpulsion ();
}

// get lower impulsion in Z
//
// return value = lower value in unit of 2 * Pi / L

inline int PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle::GetLowerImpulsionZ ()
{
  return ((PeriodicOneDOneParticle*) this->StateZ)->GetLowerImpulsion ();
}

// get the wave function basis in X direction
//
// return = pointer to 1D one particle basis

inline AbstractOneDOneParticle* PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle::GetStateX ()
{
  AbstractOneDOneParticle* stateX = (PeriodicReflexionSymmetryOneDOneParticle*) this->StateX->Clone ();
  return stateX;
}

// get the wave function basis in X direction
//
// return = pointer to 1D one particle basis

inline AbstractOneDOneParticle* PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle::GetStateY ()
{
  AbstractOneDOneParticle* stateY = (PeriodicReflexionSymmetryOneDOneParticle*) this->StateY->Clone ();
  return stateY;
}

// get the wave function basis in Z direction
//
// return = pointer to 1D one particle basis
inline AbstractOneDOneParticle* PeriodicXYReflexionSymmetryZPeriodicThreeDOneParticle::GetStateZ ()
{
  AbstractOneDOneParticle* stateZ = (PeriodicOneDOneParticle*) this->StateZ->Clone ();
  return stateZ;
}

#endif
