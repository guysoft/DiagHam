////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin chain hamiltonian with translations             //
//                                                                            //
//                        last modification : 04/03/2002                      //
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


#include "HilbertSpace/AbstractSpinChainWithTranslations.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"

#include <iostream>

using std::cout;
using std::endl;


// virtual destructor
//

AbstractSpinChainWithTranslations::~AbstractSpinChainWithTranslations () 
{
}


// return index of resulting state from application of S-_i1 S+_j1 S-_i2 S+_j2 operator on a given state
//
// i1 = position of leftmost S- operator
// j1 = position of leftmost S+ operator
// i2 = position of rightmost S- operator
// j2 = position of rightmost S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state (orbit index)

int AbstractSpinChainWithTranslations::SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::SmiSpjSmiSpj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i1 Sz_j1 S-_i2 S+_j2 operator on a given state
//
// i1 = position of first Sz operator
// j1 = position of second Sz operator
// i2 = position of S- operator
// j2 = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state (orbit index)

int AbstractSpinChainWithTranslations::SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::SmiSpjSmiSpj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChainWithTranslations::Smi (int i, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::Smi" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

int AbstractSpinChainWithTranslations::Szi (int i, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::Szi" << endl;
  return this->HilbertSpaceDimension;
}
  
// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChainWithTranslations::SpiSpj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::SpiSpj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChainWithTranslations::SmiSmj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::SmiSmj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChainWithTranslations::SpiSzj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::SpiSzj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChainWithTranslations::SmiSzj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::SmiSzj" << endl;
  return this->HilbertSpaceDimension;
}
  
// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

int AbstractSpinChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::SmiSpj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChainWithTranslations::Spi (int i, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::Spi" << endl;
  return this->HilbertSpaceDimension;
}
