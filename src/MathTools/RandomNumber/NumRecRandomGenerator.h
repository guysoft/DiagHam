
#ifndef NUMRECRANDOMGENERATOR_H
#define NUMRECRANDOMGENERATOR_H


#include "config.h"
#include "MathTools/RandomNumber/AbstractRandomNumberGenerator.h"

#include <stdlib.h>


// size of internal storage

#define NTAB 32


class NumRecRandomGenerator:public AbstractRandomNumberGenerator

{

 private:
      
  long idum;    //global seed for the random number generator
  long iy;   
  long iv[NTAB];
  long iset;
  double gset;

 public:

  // constructor
  NumRecRandomGenerator(long seed=-1);

  // copy constructor
  NumRecRandomGenerator(const NumRecRandomGenerator&);
  
  // virtual destructor
  //
  virtual ~NumRecRandomGenerator();

  // clone random number generator 
  //
  // return value = clone of the random number generator
  virtual AbstractRandomNumberGenerator* Clone ();

  // set seed of the random number generator
  //
  // seed = new seed
  virtual void SetSeed(const unsigned long& seed);
  
  // get real random number between 0 and 1
  //
  // return value = random number
  virtual double GetRealRandomNumber();

  // get standard gaussian distributed real random number
  //
  // return value = random number
  virtual double GetGaussianRandomNumber();

  // get integer random number between 0 and GetMaxInteger
  //
  // return value = random number
  unsigned long GetIntegerRandomNumber();
  
  // get maximum integer value that can be returned by GetIntegerRandomNumber
  //
  // return value = maximum integer
  unsigned long GetMaxInteger();

};

// get integer random number between 0 and GetMaxInteger
//
// return value = random number

inline unsigned long NumRecRandomGenerator::GetIntegerRandomNumber()
{
  return rand();
}

// get maximum integer value that can be returned by GetIntegerRandomNumber
//
// return value = maximum integer

inline unsigned long NumRecRandomGenerator::GetMaxInteger()
{
  return RAND_MAX;
}

#endif // NUMRECRANDOMGENERATOR_H
