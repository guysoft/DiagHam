////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             set of functions used to squeezed (aka Haldane) basis          //
//                                                                            //
//                        last modification : 24/02/2009                      //
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
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include <iostream>


using std::cout;
using std::endl;


// get the root parition from a file
// 
// rootFileName = name of the file that contains the root description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum Lz value
// referenceState = array where the root partition description will be stored
// return value = true if no error occured

bool FQHEGetRootPartition (char* rootFileName, int& nbrParticles, int& lzMax, int*& referenceState)
{
  ConfigurationParser ReferenceStateDefinition;
  if (ReferenceStateDefinition.Parse(rootFileName) == false)
    {
      ReferenceStateDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", nbrParticles) == false) || (nbrParticles <= 0))
    {
      cout << "NbrParticles is not defined or as a wrong value" << endl;
      return false;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", lzMax) == false) || (lzMax < 0))
    {
      cout << "LzMax is not defined or as a wrong value" << endl;
      return false;
    }
  int MaxNbrLz;
  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', referenceState, MaxNbrLz) == false)
    {
      cout << "error while parsing ReferenceState in " << rootFileName << endl;
      return false;     
    }
  if (MaxNbrLz != (lzMax + 1))
    {
      cout << "wrong LzMax value in ReferenceState" << endl;
      return false;     
    }
  return true;
}
