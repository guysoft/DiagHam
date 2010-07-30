////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     set of functions used to managed files related to QHE on torus         //
//                                                                            //
//                        last modification : 12/04/2010                      //
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


#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& kyMax, bool& statistics)
{
  char* StrNbrParticles;

  StrNbrParticles = strstr(filename, "_n_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  nbrParticles = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess number of particles from file name " << filename << endl;
      return false;            
    }

  StrNbrParticles = strstr(filename, "_2s_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  kyMax = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
	  else
	    StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess maximum momentum from file name " << filename << endl;
      return false;            
    }

  if (strstr(filename, "fermion") == 0)
    {
      if (strstr(filename, "boson") == 0)
	{
	  cout << "can't guess particle statistics from file name " << filename << endl;
	  return false;	  
	}
      else
	{
	  statistics = false;
	}
    }
  else
    {
      statistics = true;
    }
  return true;
}


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// kyMax = reference to the maximum momentum for a single particle
// ky = reference to the y projection of the angular momentum
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, bool& statistics)
{
  if (FQHEOnTorusFindSystemInfoFromFileName(filename, nbrParticles, kyMax, statistics) == false)
    return false;
  char* StrNbrParticles;

  StrNbrParticles = strstr(filename, "_ky_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      if (StrNbrParticles[SizeString] == '-')
	++SizeString;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '.') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  ky = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '.';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess ky momentum from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name for system suth an SU(2) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles
// kyMax = reference to the maximum momentum for a single particle
// ky = reference to the y projection of the angular momentum
// sz = reference to twice the z projection of the total spin
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured

bool FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& kyMax, int& ky, int& sz, bool& statistics)
{
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(filename, nbrParticles, kyMax, ky, statistics) == false)
    return false;
  char* StrNbrParticles;

  StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      if (StrNbrParticles[SizeString] == '-')
	++SizeString;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  ky = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess sz from file name " << filename << endl;
      return false;            
    }
  return true;
}

