////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//                     class of hilbert space of one 1d particle              //
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


#include "config.h"
#include "HilbertSpace/QuantumDotHilbertSpace/AbstractOneDOneParticle.h"


// default constructor
AbstractOneDOneParticle::AbstractOneDOneParticle()
{
}

// constructor
//
// nbrState = number of states

AbstractOneDOneParticle::AbstractOneDOneParticle(int nbrState)
{
  this->NbrState = nbrState;
  this->HilbertSpaceDimension = this->NbrState;
}

// copy constructor
//
// space = reference on Hilbert space to copy

AbstractOneDOneParticle::AbstractOneDOneParticle(const AbstractOneDOneParticle& space)
{
  this->NbrState = space.NbrState;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
}

// destructor
//

AbstractOneDOneParticle::~AbstractOneDOneParticle()
{
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

AbstractOneDOneParticle& AbstractOneDOneParticle::operator = (const AbstractOneDOneParticle& space)
{
  this->NbrState = space.NbrState;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* AbstractOneDOneParticle::Clone()
{
  return new AbstractOneDOneParticle(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> AbstractOneDOneParticle::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* AbstractOneDOneParticle::GetQuantumNumber (int index)
{
  return 0;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* AbstractOneDOneParticle::ExtractSubspace (AbstractQuantumNumber& q, 
							      SubspaceSpaceConverter& converter)
{
  return 0;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& AbstractOneDOneParticle::PrintState (ostream& Str, int state)
{
  Str << "(" << state << ")";
  return Str;
}


