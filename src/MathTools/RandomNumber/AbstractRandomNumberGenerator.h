////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of abstract random number generator                //
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


#ifndef ABSTRACTRANDOMNUMBERGENERATOR_H
#define ABSTRACTRANDOMNUMBERGENERATOR_H


#include "config.h"


class AbstractRandomNumberGenerator
{

 public:

  // virtual destructor
  //
  virtual ~AbstractRandomNumberGenerator();

  // clone random number generator 
  //
  // return value = clone of the random number generator
  virtual AbstractRandomNumberGenerator* Clone () = 0;

  // set seed of the random number generator
  //
  // seed = new seed
  virtual void SetSeed(const unsigned long& seed) = 0;
  
  // get real random number between 0 and 1
  //
  // return value = random number
  virtual double GetRealRandomNumber() = 0;

  // get integer random number between 0 and GetMaxInteger
  //
  // return value = random number
  virtual unsigned long GetIntegerRandomNumber() = 0;
  
  // get maximum integer value that can be returned by GetIntegerRandomNumber
  //
  // return value = maximum integer
  virtual unsigned long GetMaxInteger() = 0;

};

#endif
