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

#include "config.h"
#include "HilbertSpace/BosonOnSphereFusedBasis.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>

using std::cout;
using std::endl;
using std::dec;
using std::hex;


// basic constructor
//
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// normalization = normalization type 

BosonOnSphereFusedBasis::BosonOnSphereFusedBasis (int nbrBosons, int totalLz, int lzMax, int normalization)
{
  this->NbrParticles = nbrBosons;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + (this->NbrParticles * this->LzMax)) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->Normalization = normalization;
}

//  constructor from a configuration file
//
// fileName = configuration file name
// normalization = normalization type 

BosonOnSphereFusedBasis::BosonOnSphereFusedBasis (char* fileName, int normalization)
{
  bool Statistics = true;
  this->NbrParticles = 0;
  this->TotalLz = 0;
  this->LzMax = 0;
  this->GetTargetHilbertSpaceData(fileName, this->NbrParticles, this->LzMax, this->TotalLz, Statistics);
  this->ShiftedTotalLz = (this->TotalLz + (this->NbrParticles * this->LzMax)) >> 1;
  this->NbrLzValue = this->LzMax + 1;  
  this->Normalization = normalization;
}

// destructor
//

BosonOnSphereFusedBasis::~BosonOnSphereFusedBasis ()
{
}

// convert the state associated to the fused basis into a vector defined in a full basis
//
// targetBasis = pointer to the full basis
// outputVector = reference on the vector where the state will be stored
// return value = true if no error occured

bool BosonOnSphereFusedBasis::ConvertToFullBasis(BosonOnSphereShort* targetBasis, RealVector& outputVector)
{
  return true;
}

// convert the state associated to the fused basis into a vector defined in a full basis
//
// targetBasis = pointer to the full basis
// outputVector = reference on the vector where the state will be stored
// return value = true if no error occured

bool BosonOnSphereFusedBasis::ConvertToFullBasis(BosonOnSphereShort* targetBasis, LongRationalVector& outputVector)
{
  return true;
}

// get the target Hilbert space data from a configuration file
//
// inputFileName = name of the file that describes the states to fuse
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the angular momentum 
// totalLz = reference on twice the total Lz value
// statistics = reference on the particle statistics
// return value = true if no error occured

bool BosonOnSphereFusedBasis::GetTargetHilbertSpaceData(char* inputFileName, int& nbrParticles, int& lzMax, int& totalLz, bool& statistics)
{
  MultiColumnASCIIFile InputVectors;
  if (InputVectors.Parse(inputFileName) == false)
    {
      InputVectors.DumpErrors(cout) << endl;
      return false;
    }
  int NbrFusedStates = InputVectors.GetNbrColumns();
  if (((NbrFusedStates % 4) != 0) || (NbrFusedStates < 8))
    {
      cout << "wrong number of columns in " << inputFileName << endl;
      return false;
    }
  NbrFusedStates /= 4;
  nbrParticles = 0;
  lzMax = 0;
  totalLz = 0;
  for (int j = 0; j < NbrFusedStates; ++j)
    {
      int TmpNbrParticles = 0;
      int TmpLzMax = 0;
      int TmpTotalLz = 0;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(InputVectors(1 + (4 * j), 0), TmpNbrParticles, TmpLzMax, TmpTotalLz, statistics) == false)
	{
	  cout << "error while retrieving system parameters from state name " << InputVectors(1 + (4 * j), 0) << endl;
	  return false;
	}
      int LocalPadding = 0;
      if (j != (NbrFusedStates - 1))
	{
	  LocalPadding = atoi(InputVectors(4 + (4 * j), 0));
	}
      else
	{
	  LocalPadding = -1;
	}
      totalLz += ((TmpTotalLz +  (TmpLzMax * TmpNbrParticles)) >> 1) + (lzMax * TmpNbrParticles);
      nbrParticles += TmpNbrParticles;
      lzMax += TmpLzMax + LocalPadding + 1;
    }
  totalLz <<= 1;
  totalLz -= lzMax * nbrParticles;
  return true;
}

