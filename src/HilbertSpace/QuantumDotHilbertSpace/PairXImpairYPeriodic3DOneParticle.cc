////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//class of hilbert space of 3d periodic box particle with XY reflexion symmetry  //
//                                                                            //
//                        last modification : 29/03/2004                      //
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
#include "HilbertSpace/QuantumDotHilbertSpace/PairXImpairYPeriodic3DOneParticle.h"


// constructor
//
// space = a more general space without symmetry reduction

PairXImpairYPeriodic3DOneParticle::PairXImpairYPeriodic3DOneParticle(XYReflexionSymmetricPeriodic3DOneParticle& space)
{
  this->LowerImpulsionX = 0;
  this->NbrStateX = space.GetNbrCosinusStateX();
  this->LowerImpulsionY = 1;
  this->NbrStateY = space.GetNbrSinusStateY();;
  this->NbrStateZ = space.GetNbrStateZ();
  this->LowerImpulsionZ = space.GetLowerImpulsionZ();
  this->HilbertSpaceDimension = this->NbrStateX * this->NbrStateY * this->NbrStateZ;
}

// copy constructor
//
// space = reference on Hilbert space to copy

PairXImpairYPeriodic3DOneParticle::PairXImpairYPeriodic3DOneParticle(const PairXImpairYPeriodic3DOneParticle& space)
{
  this->NbrStateX = space.NbrStateX;
  this->LowerImpulsionX = space.LowerImpulsionX;
  this->NbrStateY = space.NbrStateY;
  this->LowerImpulsionY = space.LowerImpulsionY;
  this->NbrStateZ = space.NbrStateZ;
  this->LowerImpulsionZ = space.LowerImpulsionZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
}

// destructor
//

PairXImpairYPeriodic3DOneParticle::~PairXImpairYPeriodic3DOneParticle()
{
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

PairXImpairYPeriodic3DOneParticle& PairXImpairYPeriodic3DOneParticle::operator = (const PairXImpairYPeriodic3DOneParticle& space)
{
  this->NbrStateX = space.NbrStateX;
  this->LowerImpulsionX = space.LowerImpulsionX;
  this->NbrStateY = space.NbrStateY;
  this->LowerImpulsionY = space.LowerImpulsionY;
  this->NbrStateZ = space.NbrStateZ;
  this->LowerImpulsionZ = space.LowerImpulsionZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  return *this;
}
/*
// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PairXImpairYPeriodic3DOneParticle::Clone()
{
  return new PairXImpairYPeriodic3DOneParticle(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> PairXImpairYPeriodic3DOneParticle::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* PairXImpairYPeriodic3DOneParticle::GetQuantumNumber (int index)
{
  return 0;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* PairXImpairYPeriodic3DOneParticle::ExtractSubspace (AbstractQuantumNumber& q, 
							      SubspaceSpaceConverter& converter)
{
  return 0;
}
*/
// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& PairXImpairYPeriodic3DOneParticle::PrintState (ostream& Str, int state)
{
  int k1 = state / (this->NbrStateY * this->NbrStateX);
  int j1 = (state - (k1 * this->NbrStateY * this->NbrStateX)) / this->NbrStateX;
  Str << "(" << (state - (k1 * this->NbrStateY * this->NbrStateX) - (j1 * this->NbrStateX)) << ", " << j1 << ", " << k1 << ")";
  return Str;
}
