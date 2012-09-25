////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of abstract spin chain                       //
//                                                                            //
//                        last modification : 18/04/2001                      //
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


#include "HilbertSpace/AbstractSpinChain.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"


// virtual destructor
//

AbstractSpinChain::~AbstractSpinChain () 
{
}


// return matrix representation of Sx
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& AbstractSpinChain::Sxi (int i, Matrix& M)
{
  return M;
}

// return matrix representation of i * Sy
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& AbstractSpinChain::Syi (int i, Matrix& M)
{
  return M;
}

// return matrix representation of Sz
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& AbstractSpinChain::Szi (int i, Matrix& M)
{
  return M;
}

