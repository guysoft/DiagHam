////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.08                           //
//                                                                            //
//                  Copyright (C) 1998-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         a set of functions used for little<->big endian convertion         //
//                                                                            //
//                          first release : 05/06/2003                        //
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
#include <fstream>


using std::ofstream;
using std::ifstream;
using std::fstream;


// function to read Little Endian encoded variable from a file
//
// file = reference on the input file stream
// var = reference on the variable to store the result

template<class ClassName>
void ReadLittleEndian (fstream& file, ClassName& var)
{
  file.read ((char*) &var, sizeof(ClassName));
#ifdef __BIGENDIAN__
  ClassName TmpVar = var;
  unsigned char* TmpBin1 = (unsigned char*) &var;
  unsigned char* TmpBin2 = (unsigned char*) &TmpVar;
  int max = sizeof(ClassName) >> 1;
  for (int i = 0; i < max; i++)
    {
      TmpBin1[i] = TmpBin2[sizeof(ClassName) - i -1];
      TmpBin1[sizeof(ClassName) - i -1] = TmpBin2[i];
    }
#endif
}

// function to write Little Endian encoded variable from a file using 
//
// file = reference on the output file stream
// var = reference on the variable to store the result

template<class ClassName>
void WriteLittleEndian (fstream& file, ClassName& var)
{
#ifdef __BIGENDIAN__
  ClassName TmpVar = var;
  unsigned char* TmpBin2 = (unsigned char*) &var;
  unsigned char* TmpBin1 = (unsigned char*) &TmpVar;
  int max = sizeof(ClassName) >> 1;
  for (int i = 0; i < max; i++)
    {
      TmpBin1[i] = TmpBin2[sizeof(ClassName) - i -1];
      TmpBin1[sizeof(ClassName) - i -1] = TmpBin2[i];
    }
  file.write ((char*) &TmpVar, sizeof(ClassName));
#else
  file.write ((char*) &var, sizeof(ClassName));
#endif
}

// function to read Little Endian encoded variable from a file
//
// file = reference on the input file stream
// var = reference on the variable to store the result

template<class ClassName>
void ReadLittleEndian (ifstream& file, ClassName& var)
{
  file.read ((char*) &var, sizeof(ClassName));
#ifdef __BIGENDIAN__
  ClassName TmpVar = var;
  unsigned char* TmpBin1 = (unsigned char*) &var;
  unsigned char* TmpBin2 = (unsigned char*) &TmpVar;
  int max = sizeof(ClassName) >> 1;
  for (int i = 0; i < max; i++)
    {
      TmpBin1[i] = TmpBin2[sizeof(ClassName) - i -1];
      TmpBin1[sizeof(ClassName) - i -1] = TmpBin2[i];
    }
#endif
}

// function to write Little Endian encoded variable from a file using 
//
// file = reference on the output file stream
// var = reference on the variable to store the result

template<class ClassName>
void WriteLittleEndian (ofstream& file, ClassName& var)
{
#ifdef __BIGENDIAN__
  ClassName TmpVar = var;
  unsigned char* TmpBin2 = (unsigned char*) &var;
  unsigned char* TmpBin1 = (unsigned char*) &TmpVar;
  int max = sizeof(ClassName) >> 1;
  for (int i = 0; i < max; i++)
    {
      TmpBin1[i] = TmpBin2[sizeof(ClassName) - i -1];
      TmpBin1[sizeof(ClassName) - i -1] = TmpBin2[i];
    }
  file.write ((char*) &TmpVar, sizeof(ClassName));
#else
  file.write ((char*) &var, sizeof(ClassName));
#endif
}


#ifndef ENDIAN_H
#define ENDIAN_H

//#define NAIVE_ENDIAN

// function to read Little Endian encoded variable from a file
//
// file = reference on the input file stream
// var = reference on the variable to store the result

inline void ReadBlockLittleEndian (ifstream& file, double *var, long size)
{
  file.read ((char*) var, size*sizeof(double));
#ifdef __BIGENDIAN__

#ifdef NAIVE_ENDIAN
  int max = sizeof(double) >> 1;
  for (int s=0; s<size; ++s)
    {
      double TmpVar = var[s];
      unsigned char* TmpBin1 = (unsigned char*) &var;
      unsigned char* TmpBin2 = (unsigned char*) &TmpVar;
      for (int i = 0; i < max; i++)
	{
	  TmpBin1[i] = TmpBin2[sizeof(double) - i -1];
	  TmpBin1[sizeof(double) - i -1] = TmpBin2[i];
	}
    }
#else

#ifdef __64_BITS__
  unsigned long Mask3=0x00ff00ff00ff00fful;
  unsigned long Mask4=0x0000ffff0000fffful;
  unsigned long TmpVar;
  for (int s=0; s<size; ++s)
    {
      TmpVar = (unsigned long) var[s];
      TmpVar = ((TmpVar & Mask3)<< 8) | ((TmpVar >> 8) & Mask3);  // swap bytes
      TmpVar = ((TmpVar & Mask4)<< 16) | ((TmpVar >> 16) & Mask4);  // swap wydes
      var[s] = (double) ( (TmpVar << 32) | (TmpVar >> 32) );
    }
#else
  unsigned Mask3=0x00ff00fful;
  unsigned *TmpVar;
  unsigned Swap;
  for (int s=0; s<size; ++s)
    {
      TmpVar =  &(var[s]);
      TmpVar[0] = ((TmpVar[0] & Mask3)<< 8) | ((TmpVar[0] >> 8) & Mask3);  // swap bytes on word 1
      TmpVar[0] = ( (TmpVar[0] << 16) | (TmpVar[0] >> 16) ); // swap wydes on word 1
      TmpVar[1] = ((TmpVar[1] & Mask3)<< 8) | ((TmpVar[1] >> 8) & Mask3);  // swap bytes on word 2
      TmpVar[1] = ( (TmpVar[1] << 16) | (TmpVar[1] >> 16) ); // swap wydes on word 2
      Swap = TmpVar[1];
      TmpVar[1] = TmpVar[0];
      TmpVar[0] = Swap;
    }
#endif // __64_BITS__
  
#endif // NAIVE_ENDIAN

#endif // __BIGENDIAN__
}

// function to write Little Endian encoded variable from a file using 
//
// file = reference on the output file stream
// var = reference on the variable to store the result

inline void WriteBlockLittleEndian (ofstream& file, double *var, long size)
{
#ifdef __BIGENDIAN__
  for (int s=0; s<size; ++s)
    {
      double TmpVar = var[s];
      unsigned char* TmpBin2 = (unsigned char*) &var;
      unsigned char* TmpBin1 = (unsigned char*) &TmpVar;
      int max = sizeof(double) >> 1;
      for (int i = 0; i < max; i++)
	{
	  TmpBin1[i] = TmpBin2[sizeof(double) - i -1];
	  TmpBin1[sizeof(double) - i -1] = TmpBin2[i];
	}
      file.write ((char*) &TmpVar, sizeof(double));
    }
#else
  file.write ((char*) var, size*sizeof(double));
#endif
}


#endif // ENDIAN_H
