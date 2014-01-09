////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     set of functions used to managed files related to spin system          //
//                                                                            //
//                        last modification : 22/06/2010                      //
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


#include "Tools/SpinFiles/SpinFileTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins
// return value = true if no error occured

bool SpinFindSystemInfoFromFileName(char* filename, int& nbrSpins)
{
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_n_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 3;
      int SizeString = 0;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if (SizeString != 0)
	{
	  char Tmp = StrNbrSpins[SizeString];
	  StrNbrSpins[SizeString] = '\0';
	  nbrSpins = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = Tmp;
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess number of spins from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins
// spin = reference to twice the spin value per site
// return value = true if no error occured

bool SpinFindSystemInfoFromFileName(char* filename, int& nbrSpins, int& spin)
{
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_n_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 3;
      int SizeString = 0;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] != '_') && (StrNbrSpins[SizeString] >= '0') 
	     && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrSpins[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrSpins[SizeString] = '\0';
	  nbrSpins = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '_';
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess number of spins from file name " << filename << endl;
      return false;            
    }

  StrNbrSpins = strstr(filename, "spin_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 5;
      int SizeString = 0;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] != '_') && (StrNbrSpins[SizeString] >= '0') 
	     && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrSpins[SizeString] == '_') && (SizeString != 0))
	{
	  StrNbrSpins[SizeString] = '\0';
	  spin = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '_';
	  StrNbrSpins += SizeString;
	  if (strstr((StrNbrSpins + SizeString), "_2_") != (StrNbrSpins + SizeString))
	    {
	      spin *= 2;
	    }
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess spin value from file name " << filename << endl;
      return false;            
    }

  return true;
}

// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// return value = true if no error occured

bool SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin)
{
  if (SpinFindSystemInfoFromFileName(filename, nbrSpins, spin) == false)
    return false;
  char* StrNbrSpins;

  StrNbrSpins = strstr(filename, "_sz_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 4;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && (StrNbrSpins[SizeString] != '.') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((StrNbrSpins[SizeString] == '.') && (SizeString != 0))
	{
	  StrNbrSpins[SizeString] = '\0';
	  sz = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '.';
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess total Sz projection from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = file name
// nbrSpins = reference to the number of spins
// sz = reference to twice the Sz value
// spin = reference to twice the spin value per site
// momentum = reference on the momentum
// return value = true if no error occured

bool SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrSpins, int& sz, int& spin, int& momentum)
{
  if (SpinFindSystemInfoFromVectorFileName(filename, nbrSpins, sz, spin) == false)
    return false;
  char* StrMomentum;

  StrMomentum = strstr(filename, "_k_");
  if (StrMomentum != 0)
    {
      StrMomentum += 4;
      int SizeString = 0;
      if (StrMomentum[SizeString] == '-')
	++SizeString;
      while ((StrMomentum[SizeString] != '\0') && (StrMomentum[SizeString] != '_') && (StrMomentum[SizeString] >= '0') 
	     && (StrMomentum[SizeString] <= '9'))
	++SizeString;
      if ((StrMomentum[SizeString] == '_') && (SizeString != 0))
	{
	  StrMomentum[SizeString] = '\0';
	  momentum = atoi(StrMomentum);
	  StrMomentum[SizeString] = '_';
	  StrMomentum += SizeString;
	}
      else
	StrMomentum = 0;
    }
  if (StrMomentum == 0)
    {
      cout << "can't guess momentum from file name " << filename << endl;
      return false;            
    }
  return true;
}

// try to guess system information from file name
//
// filename = file name
// nbrSites = reference to the number of sites
// qValue = reference to Zn charge
// return value = true if no error occured

bool PottsFindSystemInfoFromVectorFileName(char* filename, int& nbrSites, int& qValue)
{
  if (SpinFindSystemInfoFromFileName(filename, nbrSites) == false)
    return false;
  char* StrNbrSpins;
  StrNbrSpins = strstr(filename, "_q_");
  if (StrNbrSpins != 0)
    {
      StrNbrSpins += 3;
      int SizeString = 0;
      if (StrNbrSpins[SizeString] == '-')
	++SizeString;
      while ((StrNbrSpins[SizeString] != '\0') && 
	     (StrNbrSpins[SizeString] >= '0') && (StrNbrSpins[SizeString] <= '9'))
	++SizeString;
      if ((SizeString > 1) || ((SizeString == 1) && (StrNbrSpins[0] != '-')))
	{
	  StrNbrSpins[SizeString] = '\0';
	  qValue = atoi(StrNbrSpins);
	  StrNbrSpins[SizeString] = '.';
	  StrNbrSpins += SizeString;
	}
      else
	StrNbrSpins = 0;
    }
  if (StrNbrSpins == 0)
    {
      cout << "can't guess the total Zn charge from file name " << filename << endl;
      return false;            
    }
  return true; 
}

