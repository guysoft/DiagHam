////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system defined throught          //  
//                                the product rules                           //
//                                                                            //
//                        last modification : 16/03/2011                      //
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

#ifndef BOSONONSPHEREFUSEDBASIS
#define BOSONONSPHEREFUSEDBASIS

#include "config.h"
#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "MathTools/LongRational.h"
#include "Vector/LongRationalVector.h"
#include "Vector/RealVector.h"

#include <iostream>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class BosonOnSphereFusedBasis : public AbstractQHEParticle
{

 protected:

  // number of particles
  int NbrParticles;
  // twice the momentum total value
  int TotalLz;
  // momentum total value shifted by LzMax / 2 * NbrBosons
  int ShiftedTotalLz;
  // twice the maximum Lz value reached by a boson
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;

  // one body wavefunction normalization
  int Normalization;

  // temporary state used when applying operators
  unsigned long* TemporaryState;
  int TemporaryStateLzMax;

 public:

  enum 
    {
      Unnormalized = 0x01,
      ConformalLimit = 0x02,
      Sphere = 0x03,
      Disk = 0x04
    };

  // basic constructor
  //
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // normalization = normalization type 
  BosonOnSphereFusedBasis (int nbrBosons, int totalLz, int lzMax, int normalization);

  //  constructor from a configuration file
  //
  // fileName = configuration file name
  // normalization = normalization type 
  BosonOnSphereFusedBasis (char* fileName, int normalization);

  // destructor
  //
  ~BosonOnSphereFusedBasis ();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // convert the state associated to the fused basis into a vector defined in a full basis
  //
  // targetBasis = pointer to the full basis
  // outputVector = reference on the vector where the state will be stored
  // return value = true if no error occured
  virtual bool ConvertToFullBasis(BosonOnSphereShort* targetBasis, RealVector& outputVector);

  // convert the state associated to the fused basis into a vector defined in a full basis
  //
  // targetBasis = pointer to the full basis
  // outputVector = reference on the vector where the state will be stored
  // return value = true if no error occured
  virtual bool ConvertToFullBasis(BosonOnSphereShort* targetBasis, LongRationalVector& outputVector);

 protected:

  // get the target Hilbert space data from a configuration file
  //
  // inputFileName = name of the file that describes the states to fuse
  // nbrParticles = reference on the number of particles
  // lzMax = reference on twice the angular momentum 
  // totalLz = reference on twice the total Lz value
  // statistics = reference on the particle statistics
  // return value = true if no error occured
  virtual bool GetTargetHilbertSpaceData(char* inputFileName, int& nbrParticles, int& lzMax, int& totalLz, bool& statistics);

};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnSphereFusedBasis::GetParticleStatistic()
{
  return AbstractQHEParticle::BosonicStatistic;
}

#endif

