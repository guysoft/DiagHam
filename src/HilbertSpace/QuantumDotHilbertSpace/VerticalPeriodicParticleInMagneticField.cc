////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2002-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//            class of hilbert space of one particle in magnetic field        //
//                                                                            //
//                        last modification : 04/22/2004                      //
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


#include "config.h"
#include "HilbertSpace/QuantumDotHilbertSpace/VerticalPeriodicParticleInMagneticField.h"


// constructor
//
// nbrStateR = number of Landau states in plane 
// nbrStateZ = wave function basis dimension in the z direction
// lowerImpulsionZ = LowerImpulsionZ

VerticalPeriodicParticleInMagneticField::VerticalPeriodicParticleInMagneticField(int nbrStateR, int nbrStateZ, int lowerImpulsionZ)
{
  this->NbrStateR = nbrStateR;
  this->NbrStateZ = nbrStateZ;
  this->LowerImpulsionZ = lowerImpulsionZ;
  this->HilbertSpaceDimension = this->NbrStateR * this->NbrStateZ;
}

// copy constructor
//
// space = reference on Hilbert space to copy

VerticalPeriodicParticleInMagneticField::VerticalPeriodicParticleInMagneticField(const VerticalPeriodicParticleInMagneticField& space)
{
  this->NbrStateR = space.NbrStateR;
  this->NbrStateZ = space.NbrStateZ;
  this->LowerImpulsionZ = space.LowerImpulsionZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
}

// destructor
//

VerticalPeriodicParticleInMagneticField::~VerticalPeriodicParticleInMagneticField()
{
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

VerticalPeriodicParticleInMagneticField& VerticalPeriodicParticleInMagneticField::operator = (const VerticalPeriodicParticleInMagneticField& space)
{
  this->NbrStateR = space.NbrStateR;
  this->NbrStateZ = space.NbrStateZ;
  this->LowerImpulsionZ = space.LowerImpulsionZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* VerticalPeriodicParticleInMagneticField::Clone()
{
  return new VerticalPeriodicParticleInMagneticField(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> VerticalPeriodicParticleInMagneticField::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* VerticalPeriodicParticleInMagneticField::GetQuantumNumber (int index)
{
  return 0;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* VerticalPeriodicParticleInMagneticField::ExtractSubspace (AbstractQuantumNumber& q, 
							      SubspaceSpaceConverter& converter)
{
  return 0;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& VerticalPeriodicParticleInMagneticField::PrintState (ostream& Str, int state)
{
  return Str;
}
