////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    set of functions used to managed files related to QHE on sphere         //
//                                                                            //
//                        last modification : 06/06/2006                      //
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


#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// try to guess internal symmetry group from file name
//
// su2Flag = reference on the flag that indicates if the internal symmetry group is SU(2)
// su3Flag = reference on the flag that indicates if the internal symmetry group is SU(3)
// su4Flag = reference on the flag that indicates if the internal symmetry group is SU(4)
// return value = true if no error occured

bool FQHEOnSphereFindInternalSymmetryGroupFromFileName(char* filename, bool& su2Flag, bool& su3Flag, bool& su4Flag)
{
  su2Flag = false;
  su3Flag = false;
  su4Flag = false;
  if (strstr(filename, "_su2_"))
    su2Flag = true;
  else
    if (strstr(filename, "_su3_"))
      su3Flag = true;
    else
      if (strstr(filename, "_su4_"))
	su4Flag = true;
  return true;
}

// try to guess system information from file name
//
// filename = file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereFindSystemInfoFromFileName(char* filename, int& nbrParticles, int& lzMax, bool& statistics)
{
  char* StrNbrParticles;
  if (nbrParticles == 0)
    {
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
    }
  if (lzMax == 0)
    {
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
	      lzMax = atoi(StrNbrParticles);
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
    }
  if (statistics == true)
    {
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
    }
  return true;
}

// try to guess system information from PES file name
//
// filename = file name
// nbrParticles = reference to the total number of particles (grab it only if initial value is 0)
// nbrParticlesA = reference to the number of particles na in subsystem A (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereFindSystemInfoFromPESFileName(char* filename, int& nbrParticles, int& nbrParticlesA, int& lzMax, bool& statistics)
{
  char* StrNbrParticles;
  if (nbrParticles == 0)
    {
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
    }
  if (nbrParticlesA == 0)
    {
      StrNbrParticles = strstr(filename, "_na_");
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
	      nbrParticlesA = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess number of particles in subsystem A from file name " << filename << endl;
	  return false;            
	}
    }
  if (lzMax == 0)
    {
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
	      lzMax = atoi(StrNbrParticles);
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
    }
  if (statistics == true)
    {
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
    }
  return true;
}

// try to guess system information from file name
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, bool& statistics)
{
  if (FQHEOnSphereFindSystemInfoFromFileName(filename, nbrParticles, lzMax, statistics) == false)
    return false;
  char* StrNbrParticles;
  if (lz == 0)
    {
      StrNbrParticles = strstr(filename, "_lz_");
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
	      lz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '.';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the angular momentum from file name " << filename << endl;
	  return false;            
	}
    }
  return true;
}

// try to guess system information from file name for system suth an SU(2) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// sz = reference to twice the z projection of the total spin (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, int& sz, bool& statistics)
{
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(filename, nbrParticles, lzMax, lz, statistics) == false)
    return false;
  char* StrNbrParticles;
  if (sz == 0)
    {
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
	      sz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the total spin from file name " << filename << endl;
	  return false;            
	}
    }
  return true;
}


// try to guess system information from file name for system with an SU(2) degree of freedom and discrete symmetries
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// sz = reference to twice the z projection of the total spin (grab it only if initial value is 0)
// szSymmetry = reference on the flag for the Sz<->-Sz symmetry
// szSymmetryMinusParity = reference on the flag for the minus parity sector of the Sz<->-Sz symmetry
// lzSymmetry = reference on the flag for the Lz<->-Lz symmetry
// lzSymmetryMinusParity = reference on the flag for the minus parity sector of the Lz<->-Lz symmetry
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, int& sz, 
							  bool& szSymmetry, bool& szSymmetryMinusParity, bool& lzSymmetry, bool& lzSymmetryMinusParity, bool& statistics)
{
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(filename, nbrParticles, lzMax, lz, sz, statistics) == false)
    return false;
  szSymmetry = false;
  szSymmetryMinusParity = false;
  lzSymmetry = false;
  lzSymmetryMinusParity = false;
  if (strstr(filename, "_szpsym_") != 0)
    {
      szSymmetry = true;
      return true;
    }
  if (strstr(filename, "_szmsym_") != 0)
    {
      szSymmetry = true;
      szSymmetryMinusParity = true;
      return true;
    }
  if (strstr(filename, "_lzpsym_") != 0)
    {
      lzSymmetry = true;
      return true;
    }
  if (strstr(filename, "_lzmsym_") != 0)
    {
      lzSymmetry = true;
      lzSymmetryMinusParity = true;
      return true;
    }
  if (strstr(filename, "_lzpszpsym_") != 0)
    {
      szSymmetry = true;
      lzSymmetry = true;
      return true;
    }
  if (strstr(filename, "_lzmszpsym_") != 0)
    {
      szSymmetry = true;
      lzSymmetry = true;
      lzSymmetryMinusParity = true;
      return true;
    }
  if (strstr(filename, "_lzpszmsym_") != 0)
    {
      szSymmetry = true;
      szSymmetryMinusParity = true;
      lzSymmetry = true;
      return true;
    }
  if (strstr(filename, "_lzmszmsym_") != 0)
    {
      szSymmetry = true;
      szSymmetryMinusParity = true;
      lzSymmetry = true;
      lzSymmetryMinusParity = true;
      return true;
    }
  return true;
  
}


// try to guess system information from file name for system with an SU(3) degree of freedom and discrete symmetries
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// tz = reference to twice the Tz quantum number (grab it only if initial value is 0)
// y = reference to three time the Y quantum number (grab it only if initial value is 0)
// tzSymmetry = reference on the flag for the Tz<->-Tz symmetry
// tzSymmetryMinusParity = reference on the flag for the minus parity sector of the Y<->-Y symmetry
// ySymmetry = reference on the flag for the Z3 symmetry
// lzSymmetry = reference on the flag for the Lz<->-Lz symmetry
// lzSymmetryMinusParity = reference on the flag for the minus parity sector of the Lz<->-Lz symmetry
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured
bool FQHEOnSphereWithSU3SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, int& tz, int& y,
							     bool& tzSymmetry, bool& tzSymmetryMinusParity,  bool& ySymmetry, 
							     bool& lzSymmetry, 
							     bool& lzSymmetryMinusParity, bool& statistics)
{
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(filename, nbrParticles, lzMax, lz, statistics) == false)
    return false;
  char* StrNbrParticles;
  if (tz == 0)
    {
      StrNbrParticles = strstr(filename, "_tz_");
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
	      tz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
    }
  if (y == 0)
    {
      StrNbrParticles = strstr(filename, "_y_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 3;
	  int SizeString = 0;
	  if (StrNbrParticles[SizeString] == '-')
	    ++SizeString;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      y = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the total spin from file name " << filename << endl;
	  return false;            
	}
    }
  tzSymmetry = false;
  tzSymmetryMinusParity = false;
  ySymmetry = false;
  lzSymmetry = false;
  lzSymmetryMinusParity = false;
  if (strstr(filename, "tzp") != 0)
    {
      tzSymmetry = true;
    }
  if (strstr(filename, "tzm") != 0)
    {
      tzSymmetry = true;
      tzSymmetryMinusParity = true;
    }
  if (strstr(filename, "z3sym") != 0)
    {
      ySymmetry = true;
    }
  if (strstr(filename, "lzp") != 0)
    {
      lzSymmetry = true;
    }
  if (strstr(filename, "lzm") != 0)
    {
      lzSymmetry = true;
      lzSymmetryMinusParity = true;
    }
  return true;
}

// try to guess system information from file name for system suth an SU(4) degree of freedom
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// lz = reference to twice the z projection of the angular momentum (grab it only if initial value is 0)
// sz = reference to twice the z projection of the total spin (grab it only if initial value is 0)
// iz = reference to twice the z projection of the total isospin (grab it only if initial value is 0)
// pz = reference to twice the entanglement (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false fro bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnSphereWithSU4SpinFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& lzMax, int& lz, int& sz, int& iz, int& pz, bool& statistics)
{
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(filename, nbrParticles, lzMax, lz, sz, statistics) == false)
    return false;
  char* StrNbrParticles;
  if (iz == 0)
    {
      StrNbrParticles = strstr(filename, "_iz_");
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
	      sz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the total isospin from file name " << filename << endl;
	  return false;            
	}
    }
  if (pz == 0)
    {
      StrNbrParticles = strstr(filename, "_pz_");
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
	      sz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess the entropy from file name " << filename << endl;
	  return false;            
	}
    }
  return true;
}

// try to guess system information from file name for a system of bosons on the 4D sphere
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// jz = reference to twice the z projection of the  jz angular momentum (grab it only if initial value is 0)
// kz = reference to twice the z projection of the  kz angular momentum (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOn4DSphereFindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrFluxQuanta, int& jz, int& kz, bool& statistics)
{
  if (FQHEOnSphereFindSystemInfoFromFileName(filename, nbrParticles, nbrFluxQuanta, statistics) == false)
    return false;
  
  char* StrNbrParticles;
  if (jz == 0)
    {
      StrNbrParticles = strstr(filename, "_jz_");
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
	      jz = atoi(StrNbrParticles);
	      cout << "jz = " << jz << endl;
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the jz angular momentum from file name " << filename << endl;
	  return false;            
	}
    }
  if (kz == 0)
    {
      StrNbrParticles = strstr(filename, "_kz_");
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
	      kz = atoi(StrNbrParticles);
	      cout << "kz = " << kz << endl;
	      StrNbrParticles[SizeString] = '.';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the kz angular momentum from file name " << filename << endl;
	  return false;            
	}
    }
  return true;
}

// try to guess system information from file name for a system of bosons on the CP2
//
// filename = vector file name
// nbrParticles = reference to the number of particles (grab it only if initial value is 0)
// nbrFluxQuanta = reference to the number of flux of quanta
// tz = reference to twice the z projection of the  tz angular momentum (grab it only if initial value is 0)
// y = reference to three times the z projection of the  y angular momentum (grab it only if initial value is 0)
// tzSymmetry = reference on the flag for the Tz<->-Tz symmetry
// tzSymmetryMinusParity = reference on the flag for the minus parity sector of the Tz<->-Tz symmetry
//tzZ3Symmetry = reference on the flag of the permutation symmetry 
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOnCP2FindSystemInfoFromVectorFileName(char* filename, int& nbrParticles, int& nbrFluxQuanta, int& tz, int& y, bool& tzSymmetry, bool& tzSymmetryMinusParity, bool& tzZ3Symmetry, bool& statistics)
{
  if (FQHEOnSphereFindSystemInfoFromFileName(filename, nbrParticles, nbrFluxQuanta, statistics) == false)
    return false;
  
  char* StrNbrParticles;
  if (tz == 0)
    {
      StrNbrParticles = strstr(filename, "_tz_");
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
	      tz = atoi(StrNbrParticles);
	      cout << "tz = " << tz << endl;
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the tz angular momentum from file name " << filename << endl;
	  return false;            
	}
    }
  if (y == 0)
    {
      StrNbrParticles = strstr(filename, "_y_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 3;
	  int SizeString = 0;
	  if (StrNbrParticles[SizeString] == '-')
	    ++SizeString;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '.') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      y = atoi(StrNbrParticles);
	      cout << "y = " << y << endl;
	      StrNbrParticles[SizeString] = '.';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the y angular momentum from file name " << filename << endl;
	  return false;            
	}
    }
  tzSymmetry = false;
  tzSymmetryMinusParity = false;
  tzZ3Symmetry = false;
  if (strstr(filename, "tzp") != 0)
    {
      tzSymmetry = true;
    }
  if (strstr(filename, "tzm") != 0)
    {
      tzSymmetry = true;
      tzSymmetryMinusParity = true;
    }
  if (strstr(filename, "tzZ3") != 0)
    {
      tzSymmetry = true;
      tzZ3Symmetry = true;
    }
  return true;
}

// try to guess system information from PES vector file name for a system of bosons on the 4D sphere
//
// filename = vector file name
// nbrParticles = reference to the total number of particles (grab it only if initial value is 0)
// lzMax = reference to twice the maximum momentum for a single particle (grab it only if initial value is 0)
// jz = reference to twice the z projection of the  jz angular momentum (grab it only if initial value is 0)
// kz = reference to twice the z projection of the  kz angular momentum (grab it only if initial value is 0)
// statistics = reference to flag for fermionic statistics (true for fermion, false for bosons, grab it only if initial value is true)
// return value = true if no error occured

bool FQHEOn4DSphereFindSystemInfoFromPESVectorFileName(char* filename, int& nbrParticles, int& nbrParticlesA, int& nbrFluxQuanta, int& jz, int& kz, bool& statistics)
{
  if (FQHEOnSphereFindSystemInfoFromPESFileName(filename, nbrParticles, nbrParticlesA, nbrFluxQuanta, statistics) == false)
    return false;
  
  char* StrNbrParticles;
  if (jz == 0)
    {
      StrNbrParticles = strstr(filename, "_jza_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 5;
	  int SizeString = 0;
	  if (StrNbrParticles[SizeString] == '-')
	    ++SizeString;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      jz = atoi(StrNbrParticles);
	      cout << "jza = " << jz << endl;
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the jza angular momentum from file name " << filename << endl;
	  return false;            
	}
    }
  if (kz == 0)
    {
      StrNbrParticles = strstr(filename, "_kza_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 5;
	  int SizeString = 0;
	  if (StrNbrParticles[SizeString] == '-')
	    ++SizeString;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '.') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      kz = atoi(StrNbrParticles);
	      cout << "kza = " << kz << endl;
	      StrNbrParticles[SizeString] = '.';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess z projection of the kza angular momentum from file name " << filename << endl;
	  return false;            
	}
    }
  return true;
}
