////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             set of functions used to managed files related to FQHE         //
//                          on square lattice disk                            //
//                                                                            //
//                        last modification : 01/03/2011                      //
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


#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSquareLatticeFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, bool& statistics)
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
	  StrNbrParticles[SizeString] = '\0';
	  nbrSiteX = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSiteX = 0;
    }
  StrNbrParticles = strstr(filename, "_y_");
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
	  nbrSiteY = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSiteY = 0;
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
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_kx_");
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
	  momentumX = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      momentumX = 0;
    }
  StrNbrParticles = strstr(filename, "_ky_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  momentumY = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '.';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      momentumY = 0;
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// totalSz = reference to the Sz value
// statistics = reference to flag for fermionic statistics
// return value = true if no error occured
bool FQHEOnSquareLatticeWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, int& totalSz, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, momentumX, momentumY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles = strstr(filename, "_sz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  totalSz = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '.';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
     totalSz  = 0;
    }
  return true;
}

// try to guess system information from file name 
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumX = reference to the momentum along the y direction
// mass = reference to the mass term
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& momentumX, int& momentumY, double& mass, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, momentumX, momentumY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_m_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (SizeString != 0)
	{
	  char TmpChar = StrNbrParticles[SizeString];
	  StrNbrParticles[SizeString] = '\0';
	  mass = atof(StrNbrParticles);
	  StrNbrParticles[SizeString] = TmpChar;
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      mass = 0.0;
    }
  return true;
}

// try to guess system information from file name for a cubic lattice
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// nbrSiteZ = reference to the number sites along the z direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnCubicLatticeFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& nbrSiteZ, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles = strstr(filename, "_z_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 3;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if (((StrNbrParticles[SizeString] == '.') || (StrNbrParticles[SizeString] == '_')) && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  nbrSiteZ = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSiteZ = 0;
    }
  return true;
}

// try to guess system information from file name for a cubic lattice
//
// filename = vector file name
// nbrParticles = reference to the number of particles 
// nbrSiteX = reference to the number sites along the x direction
// nbrSiteY = reference to the number sites along the y direction
// nbrSiteZ = reference to the number sites along the y direction
// momentumX = reference to the momentum along the x direction
// momentumY = reference to the momentum along the y direction
// momentumZ = reference to the momentum along the z direction
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnCubicLatticeFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrSiteX, int& nbrSiteY, int& nbrSiteZ, int& momentumX, int& momentumY, int& momentumZ, bool& statistics)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(filename, nbrParticles, nbrSiteX, nbrSiteY, momentumX, momentumY, statistics) == false)
    {
      return false;
    }
  char* StrNbrParticles;
  StrNbrParticles = strstr(filename, "_ky_");
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
	  momentumY = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      momentumY = 0;
    }
  StrNbrParticles = strstr(filename, "_kz_");
  if (StrNbrParticles != 0)
    {
      StrNbrParticles += 4;
      int SizeString = 0;
      while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
	     && (StrNbrParticles[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrParticles[SizeString] == '.') && (SizeString != 0))
	{
	  StrNbrParticles[SizeString] = '\0';
	  momentumZ = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '.';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      momentumZ = 0;
    }
  StrNbrParticles = strstr(filename, "_z_");
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
	  nbrSiteZ = atoi(StrNbrParticles);
	  StrNbrParticles[SizeString] = '_';
	  StrNbrParticles += SizeString;
	}
      else
	StrNbrParticles = 0;
    }
  if (StrNbrParticles == 0)
    {
      nbrSiteZ = 0;
    }
  return true;
}

