////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//              class of hilbert space of one 1d periodic box particle        //
//                                                                            //
//                        last modification : 05/07/2004                      //

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
#include "HilbertSpace/QuantumDotHilbertSpace/Periodic1DOneParticle.h"


// default constructor
//

Periodic1DOneParticle::Periodic1DOneParticle()
{
}

// constructor
//
// nbrState = wave function basis dimension
// low = lower impulsion

Periodic1DOneParticle::Periodic1DOneParticle(int nbrState, int low)
{
  this->NbrState = nbrState;
  this->LowerImpulsion = low;
  this->HilbertSpaceDimension = this->NbrState;
}

// copy constructor
//
// space = reference on Hilbert space to copy

Periodic1DOneParticle::Periodic1DOneParticle(const Periodic1DOneParticle& space)
{
  this->NbrState = space.NbrState;
  this->LowerImpulsion = space.LowerImpulsion;
}

// destructor
//

Periodic1DOneParticle::~Periodic1DOneParticle()
{
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

Periodic1DOneParticle& Periodic1DOneParticle::operator = (const Periodic1DOneParticle& space)
{
  this->NbrState = space.NbrState;
  this->LowerImpulsion = space.LowerImpulsion;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Periodic1DOneParticle::Clone()
{
  return new Periodic1DOneParticle(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Periodic1DOneParticle::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* Periodic1DOneParticle::GetQuantumNumber (int index)
{
  return 0;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Periodic1DOneParticle::ExtractSubspace (AbstractQuantumNumber& q, 
							      SubspaceSpaceConverter& converter)
{
  return 0;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Periodic1DOneParticle::PrintState (ostream& Str, int state)
{
  
  Str << "(" << state << ")";
  return Str;
}


