////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     set of functions used to managed files related to Hubbard models       //
//                                                                            //
//                        last modification : 19/06/2014                      //
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


#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference to the number sites
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// return value = true if no error occured

bool FTIHubbardModelFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSites, bool& statistics)
{
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_n_");
  if (StrNbrParticles == 0)
    StrNbrParticles = strstr(filename, "_p_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrParticles = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
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
  StrNbrParticles = strstr(filename, "_x_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  nbrSites = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSites = 0;
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
// nbrSites = reference to the number sites
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference to flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured

bool FTIHubbardModelFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, bool& statistics, bool& gutzwiller)
{
  if (FTIHubbardModelFindSystemInfoFromFileName(filename, nbrParticles, nbrSites, statistics) == false)
    {
      return false;
    }
  char* GutzwillerFlag = strstr(filename, "_gutzwiller_");
  if (GutzwillerFlag != 0)
    gutzwiller = true;
  else
    gutzwiller = false;
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSites = reference on the number sites
// xMomentum = reference on the momentum sector in the x direction
// xPeriodicity = reference on the periodicity in the x direction with respect to site numbering 
// statistics = reference on flag for fermionic statistics (true for fermion, false for bosons)
// gutzwiller = reference on flag  that indicated if the Gutzwiller projection was implemented within the Hilbert space
// return value = true if no error occured
bool FTIHubbardModelWith1DTranslationFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSites, int& xMomentum, int& xPeriodicity, bool& statistics, bool& gutzwiller)
{
  if (FTIHubbardModelFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSites, statistics, gutzwiller) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_kx_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '_') || (StrNbrParticles[SizeString] == '.')) && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  xMomentum = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess momentum sector from file name " << filename << endl;
      return false;            
    }
  StrNbrParticles = strstr(filename, "_xmomentum_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 11;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	{
          char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  xPeriodicity = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      cout << "can't guess periodicity from file name " << filename << endl;
      return false;            
    }
  return true;
}

