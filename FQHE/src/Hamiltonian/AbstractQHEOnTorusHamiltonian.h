////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                         to particles on a torus with                       //
//                                                                            //
//                        last modification : 27/06/2003                      //
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


#ifndef ABSTRACTQHEONTORUSHAMILTONIAN_H
#define ABSTRACTQHEONTORUSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorus.h"
#include "Hamiltonian/AbstractQHEOnSphereFullHamiltonian.h"

#include <iostream>


using std::ostream;


class AbstractArchitecture;


class AbstractQHEOnTorusHamiltonian : public AbstractQHEOnSphereFullHamiltonian
{

  friend class QHEParticlePrecalculationOperation;
  
 protected:
  
  // ratio between the width in the x direction and the width in the y direction
  double Ratio;
  // ratio between the width in the y direction and the width in the x direction
  double InvRatio;

  // number of Lz values in a state
  int NbrLzValue;
  
 public:

  // default constructor
  //
  AbstractQHEOnTorusHamiltonian();

  // destructor
  //
  virtual ~AbstractQHEOnTorusHamiltonian() = 0;

 protected:

};


#endif
