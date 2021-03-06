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


#include "config.h"
#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include <sys/time.h>


// constructor
//
// seed = optional seed definition (0 if none)

StdlibRandomNumberGenerator::StdlibRandomNumberGenerator(const unsigned int& seed)
{
  this->NbrGeneratedNumbers = 0ul;
  if (seed != 0)
    srand (seed);
  this->iset = 0;
}

// copy constructor
//
// generator = generator to copy

StdlibRandomNumberGenerator::StdlibRandomNumberGenerator(const StdlibRandomNumberGenerator& generator)
{
  this->NbrGeneratedNumbers = generator.NbrGeneratedNumbers;
  this->iset=generator.iset;
}

// destructor
//

StdlibRandomNumberGenerator::~StdlibRandomNumberGenerator()
{
}

// clone random number generator 
//
// return value = clone of the random number generator

AbstractRandomNumberGenerator* StdlibRandomNumberGenerator::Clone ()
{
  return new StdlibRandomNumberGenerator(*this);
}

// set seed of the random number generator to system time
//

void StdlibRandomNumberGenerator::UseTimeSeed()
{
  struct timeval TmpTime;
  gettimeofday(&TmpTime,NULL);
  srand((unsigned int) (TmpTime.tv_sec + TmpTime.tv_usec));
}

