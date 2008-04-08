////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of random number generator from a external generated file       //
//                                                                            //
//                        last modification : 08/04/2008                      //
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


#ifndef FILERANDOMNUMBERGENERATOR_H
#define FILERANDOMNUMBERGENERATOR_H


#include "config.h"
#include "MathTools/RandomNumber/AbstractRandomNumberGenerator.h"
#include "GeneralTools/GarbageFlag.h"


class FileRandomNumberGenerator: public AbstractRandomNumberGenerator
{

 protected:
  
  // array that contains all the real random numbers
  double* RealRandomNumbers;
  // array that contains the integer random number couterpart of RealRandomNumbers
  unsigned long* IntegerRandomNumbers;

  // numvber of random numbers available
  int NbrRandomNumbers;
  // index of the current random number
  int Position;

  // largest integer that can be generated
  unsigned long MaxInteger;

  // garbage flag to avoid memory duplication
  GarbageFlag Flag;

 public:

  // constructor
  //
  // inputFile = name of the file where random numbers are stored
  // nbrRandomNumbers = number of random numbers to read from inputFile
  // seekPosition = position of random numbers to read
  FileRandomNumberGenerator(char* inputFile, int nbrRandomNumbers, int seekPosition);

  // constructor from another random generator, used to generate a file of random numbers
  //
  // generator = generator to use to generate the random numbers
  // nbrRandomNumbers = number of random numbers to generate and to write into outputFile
  // outputFile = name of the file where random numbers will be stored
  FileRandomNumberGenerator(AbstractRandomNumberGenerator* generator, int nbrRandomNumbers, char* outputFile);

  // copy constructor
  //
  // generator = generator to copy
  FileRandomNumberGenerator(const FileRandomNumberGenerator& generator);

  // destructor
  //
  ~FileRandomNumberGenerator();

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

inline void FileRandomNumberGenerator::SetSeed(const unsigned long& seed)
{
}
  
// get real random number between 0 and 1
//
// return value = random number

inline double FileRandomNumberGenerator::GetRealRandomNumber()
{
  return this->RealRandomNumbers[this->Position++];
}

// get integer random number between 0 and GetMaxInteger
//
// return value = random number

inline unsigned long FileRandomNumberGenerator::GetIntegerRandomNumber()
{
  return this->IntegerRandomNumbers[this->Position++];
}

// get maximum integer value that can be returned by GetIntegerRandomNumber
//
// return value = maximum integer

inline unsigned long FileRandomNumberGenerator::GetMaxInteger()
{
  return this->MaxInteger;
}

#endif
