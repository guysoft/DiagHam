////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//                     class of hilbert space of one 3d particle              //
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
#include "HilbertSpace/QuantumDotHilbertSpace/AbstractTwoDOneParticle.h"


// default constructor
//

AbstractTwoDOneParticle::AbstractTwoDOneParticle ()
{
}

// constructor
//
// nbrStateX = wave function basis dimension in the x direction
// nbrStateY = wave function basis dimension in the y direction

AbstractTwoDOneParticle::AbstractTwoDOneParticle (int nbrStateX, int nbrStateY)
{
  this->StateX = new AbstractOneDOneParticle (nbrStateX);
  this->StateY = new AbstractOneDOneParticle (nbrStateY);
  this->HilbertSpaceDimension = this->StateX->GetNbrState () * this->StateY->GetNbrState ();
}

// copy constructor
//
// space = reference on Hilbert space to copy

AbstractTwoDOneParticle::AbstractTwoDOneParticle (const AbstractTwoDOneParticle& space)
{
  this->StateX = space.StateX;
  this->StateY = space.StateY;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
}

// destructor
//

AbstractTwoDOneParticle::~AbstractTwoDOneParticle ()
{
  delete this->StateX;
  delete this->StateY;
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

AbstractTwoDOneParticle& AbstractTwoDOneParticle::operator = (const AbstractTwoDOneParticle& space)
{
  this->StateX = space.StateX;
  this->StateY = space.StateY;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* AbstractTwoDOneParticle::Clone ()
{
  return new AbstractTwoDOneParticle (*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> AbstractTwoDOneParticle::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* AbstractTwoDOneParticle::GetQuantumNumber (int index)
{
  return 0;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* AbstractTwoDOneParticle::ExtractSubspace (AbstractQuantumNumber& q, 
							      SubspaceSpaceConverter& converter)
{
  return 0;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& AbstractTwoDOneParticle::PrintState (ostream& Str, int state)
{ 
  int NbrStateY = this->StateY->GetNbrState ();
  int i1 = state / NbrStateY;
  Str << "(" << i1 << ", " << (state - (i1 * NbrStateY)) << ")";
  return Str;
}


