////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of random number generator based on the stdlib rand function      //
//                                                                            //
//                        last modification : 15/09/2004                      //
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


#ifndef RANLUXRANDOMNUMBERGENERATOR_H
#define RANLUXRANDOMNUMBERGENERATOR_H


#include "config.h"
#include "MathTools/RandomNumber/AbstractRandomNumberGenerator.h"

#include <stdlib.h>


#define INV_RAND_MAX 1.0 / (1.0 + RAND_MAX)


class RanluxRandomNumberGenerator: public AbstractRandomNumberGenerator
{

 public:

  // constructor
  //
  // seed = optional seed definition (0 if none)
  RanluxRandomNumberGenerator(unsigned int& seed);

  // copy constructor
  //
  // generator = generator to copy
  RanluxRandomNumberGenerator(const RanluxRandomNumberGenerator& generator);

  // destructor
  //
  ~RanluxRandomNumberGenerator();

  // clone random number generator 
  //
  // return value = clone of the random number generator
  AbstractRandomNumberGenerator* Clone ();

  // set seed of the random number generator
  //
  // seed = new seed
  void SetSeed(const unsigned long& seed);
  
  // get real random number between 0 and 1
  //
  // return value = random number
  double GetRealRandomNumber();

  // get integer random number between 0 and GetMaxInteger
  //
  // return value = random number
  unsigned long GetIntegerRandomNumber();
  
  // get maximum integer value that can be returned by GetIntegerRandomNumber
  //
  // return value = maximum integer
  unsigned long GetMaxInteger();

};

// set seed of the random number generator
//
// seed = new seed

inline void RanluxRandomNumberGenerator::SetSeed(const unsigned long& seed)
{
  srand((unsigned int) seed);
}
  
// get real random number between 0 and 1
//
// return value = random number

inline double RanluxRandomNumberGenerator::GetRealRandomNumber()
{
  return (rand() * INV_RAND_MAX);
}

// get integer random number between 0 and GetMaxInteger
//
// return value = random number

inline unsigned long RanluxRandomNumberGenerator::GetIntegerRandomNumber()
{
  return rand();
}

// get maximum integer value that can be returned by GetIntegerRandomNumber
//
// return value = maximum integer

inline unsigned long RanluxRandomNumberGenerator::GetMaxInteger()
{
  return RAND_MAX;
}

#endif
