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
#include "HilbertSpace/QuantumDotHilbertSpace/AbstractThreeDOneParticle.h"

using std::cout;
using std::endl;


// default constructor
//

AbstractThreeDOneParticle::AbstractThreeDOneParticle ()
{
}

// constructor
//
// nbrStateX = wave function basis dimension in the x direction
// nbrStateY = wave function basis dimension in the y direction
// nbrStateZ = wave function basis dimension in the z direction

AbstractThreeDOneParticle::AbstractThreeDOneParticle (int nbrStateX, int nbrStateY, int nbrStateZ)
{
  this->StateX = new AbstractOneDOneParticle (nbrStateX);
  this->StateY = new AbstractOneDOneParticle (nbrStateY);
  this->StateZ = new AbstractOneDOneParticle (nbrStateZ);
  this->HilbertSpaceDimension = this->StateX->GetNbrState () * this->StateY->GetNbrState () * this->StateZ->GetNbrState ();
}

// copy constructor
//
// space = reference on Hilbert space to copy

AbstractThreeDOneParticle::AbstractThreeDOneParticle (const AbstractThreeDOneParticle& space)
{
  this->StateX = space.StateX;
  this->StateY = space.StateY;
  this->StateZ = space.StateZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
}

// destructor
//

AbstractThreeDOneParticle::~AbstractThreeDOneParticle ()
{
  //cout << "AbstractThreeDOneParticle destructor is being called ..." << endl;
  delete this->StateX;
  delete this->StateY;
  delete this->StateZ;
  //cout << "AbstractThreeDOneParticle destructor is finished ..." << endl;
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

AbstractThreeDOneParticle& AbstractThreeDOneParticle::operator = (const AbstractThreeDOneParticle& space)
{
  this->StateX = space.StateX;
  this->StateY = space.StateY;
  this->StateZ = space.StateZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* AbstractThreeDOneParticle::Clone ()
{
  return new AbstractThreeDOneParticle (*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> AbstractThreeDOneParticle::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* AbstractThreeDOneParticle::GetQuantumNumber (int index)
{
  return 0;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* AbstractThreeDOneParticle::ExtractSubspace (AbstractQuantumNumber& q, 
							      SubspaceSpaceConverter& converter)
{
  return 0;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& AbstractThreeDOneParticle::PrintState (ostream& Str, int state)
{ 
  int NbrStateY = this->StateY->GetNbrState (); int NbrStateZ = this->StateZ->GetNbrState (); 
  int i1 = state / (NbrStateY * NbrStateZ);
  int j1 = (state - (i1 * NbrStateY * NbrStateZ)) / NbrStateZ;
  Str << "(" << i1 << ", " << j1 << ", " << (state - (i1 * NbrStateY * NbrStateZ) - (j1 * NbrStateZ)) << ")";
  return Str;
}


