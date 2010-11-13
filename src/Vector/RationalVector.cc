////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class for n dimensional rational vector                   //
//                                                                            //
//                        last modification : 04/11/2010                      //
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


#include "Vector/RationalVector.h"
#include "GeneralTools/Endian.h"

#include <fstream>
#include <iostream>
#include <cstdlib>

using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;


// default constructor  
// 

RationalVector::RationalVector()
{
  this->VectorType = Vector::RationalData;
  this->Dimension = 0;
  this->TrueDimension = 0;
  this->LargeDimension = 0l;
  this->LargeTrueDimension = 0l;
  this->VectorId = 0;
  this->Numerators = 0;
  this->Denominators = 0;
}

// default constructor for an empty vector
// 
// size = Vector Dimension 
// zeroFlag = true if all coordinates have to be set to zero

RationalVector::RationalVector(int size, bool zeroFlag)
{
  this->VectorType = Vector::RationalData;
  this->Dimension = size;
  this->TrueDimension = size;
  this->LargeDimension = (long) size;
  this->LargeTrueDimension = (long) size;
  this->VectorId = 0;
  this->Numerators = new long [this->LargeDimension];
  this->Denominators = new long [this->LargeDimension];
  this->Flag.Initialize();
  if (zeroFlag == true)
    {
      for (long i = 0l; i < this->LargeTrueDimension; ++i)
	{
	  this->Numerators[i] = 0l;
	  this->Denominators[i] = 1l;
	}
    }
  else
    {
      for (long i = 0l; i < this->LargeTrueDimension; ++i)
	{
	  this->Denominators[i] = 1l;
	}
    }
}

// constructor for an empty rational vector bigger than 2^31
//
// size = Vector Dimension 
// zeroFlag = true if all coordinates have to be set to zero

RationalVector::RationalVector(long size, bool zeroFlag)
{
  this->VectorType = Vector::RationalData;
  this->LargeDimension = size;
  this->LargeTrueDimension = this->LargeDimension;
#ifdef  __64_BITS__
  if (this->LargeDimension < (1l << 31))
    this->Dimension = (int) size;
  else
    {
      this->Dimension = -1;
      this->VectorType |= Vector::LargeData;
    }
#else
  this->Dimension = (int) size;
#endif
  this->TrueDimension = this->Dimension;
  this->Numerators = new long [this->LargeDimension];
  this->Denominators = new long [this->LargeDimension];
  this->Flag.Initialize();
  if (zeroFlag == true)
    {
      for (long i = 0l; i < this->LargeTrueDimension; ++i)
	{
	  this->Numerators[i] = 0l;
	  this->Denominators[i] = 1l;
	}
    }
  else
    {
      for (long i = 0l; i < this->LargeTrueDimension; ++i)
	{
	  this->Denominators[i] = 1l;
	}
    }
}

// copy constructor
//
// vector = vector to copy
// duplicateFlag = true if datas have to be duplicated

RationalVector::RationalVector(const RationalVector& vector, bool duplicateFlag)
{
  this->VectorType = vector.VectorType;
  this->VectorId = vector.VectorId;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->LargeDimension = vector.LargeDimension;
  this->LargeTrueDimension = vector.LargeTrueDimension;
  if (duplicateFlag == false)
    {
      this->Numerators = vector.Numerators;
      this->Denominators = vector.Denominators;
      this->Flag = vector.Flag;
    }
  else
    {
      if (vector.Dimension > 0)
	{
	  this->Flag.Initialize();
	  this->Numerators = new long [this->TrueDimension + 1]; 
	  this->Denominators = new long [this->TrueDimension + 1]; 
	  for (int i = 0; i < this->Dimension; i++)
	    {
	      this->Numerators[i] = vector.Numerators[i];
	      this->Denominators[i] = vector.Denominators[i];
	    }
	}
      else
	{
	  if (this->LargeDimension > 0l)
	    {
	      this->Flag.Initialize();
	      this->Numerators = new long [this->LargeTrueDimension + 1]; 
	      this->Denominators = new long [this->LargeTrueDimension + 1]; 
	      for (long i = 0; i < this->LargeDimension; i++)
		{
		  this->Numerators[i] = vector.Numerators[i];
		  this->Denominators[i] = vector.Denominators[i];
		}
	    }
	  else
	    {
	      this->Numerators = 0;
	      this->Denominators = 0;
	    }
	}
    }
}

// destructor
//

RationalVector::~RationalVector ()
{
  if (((this->Dimension != 0) || (this->LargeDimension != 0l)) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Numerators;
      delete[] this->Denominators;
    }
}

// Resize vector
//
// dimension = new dimension

void RationalVector::Resize (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      this->LargeDimension = dimension;
      return;
    }
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->LargeDimension = dimension;
  this->LargeTrueDimension = dimension;
  long* TmpNumerators = new long [this->TrueDimension + 1];
  long* TmpDenominators = new long [this->TrueDimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    {
      TmpNumerators[i] = this->Numerators[i];
      TmpDenominators[i] = this->Denominators[i];
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Numerators;
      delete[] this->Denominators;
    }
  this->Numerators = TmpNumerators;
  this->Denominators = TmpDenominators;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
}

// Resize vector and set to zero all components that have been added
//
// dimension = new dimension

void RationalVector::ResizeAndClean (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      this->LargeDimension = dimension;
      return;
    }
  long* TmpNumerators = new long [dimension + 1];
  long* TmpDenominators = new long [dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    {
      TmpNumerators[i] = this->Numerators[i];
      TmpDenominators[i] = this->Denominators[i];
    }
  for (int i = this->Dimension; i < dimension; i++)
    {
      TmpNumerators[i] = 0l;
      TmpDenominators[i] = 1l;
    }
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->LargeDimension = dimension;
  this->LargeTrueDimension = dimension;
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Numerators;
      delete[] this->Denominators;
    }
  this->Numerators = TmpNumerators;
  this->Denominators = TmpDenominators;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
}

// put all vector components to zero
//
// return value = reference on current vector

Vector& RationalVector::ClearVector ()
{
  for (long i = 0l; i < this->LargeTrueDimension; ++i)
    {
      this->Numerators[i] = 0l;
      this->Denominators[i] = 1l;
    }
  return *this;
}

// put select vector components to zero
// start = start index
// nbrComponent = number of components to set to zero
// return value = reference on current vector

Vector& RationalVector::ClearVectorSegment (long start, long nbrComponent)
{
  nbrComponent += start;
  for (; start < nbrComponent; ++start)
    {
      this->Numerators[start] = 0l;
      this->Denominators[start] = 1l;
    }
  return *this;
}

// create a new vector with same size and same type but non-initialized components
//
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* RationalVector::EmptyClone(bool zeroFlag)
{
  return new RationalVector(this->LargeDimension, zeroFlag);
}

// create an array of new vectors with same size and same type but non-initialized components
//
// nbrVectors = number of vectors to sreate
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to the array of new vectors
  
Vector* RationalVector::EmptyCloneArray(int nbrVectors, bool zeroFlag)
{
  RationalVector* TmpVectors = new RationalVector [nbrVectors];
  for (int i = 0; i < nbrVectors; ++i)
    TmpVectors[i] = RationalVector(this->LargeDimension, zeroFlag);
  return TmpVectors;
}

// write vector in a file 
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool RationalVector::WriteVector (const char* fileName)
{
  this->Localize();
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  WriteLittleEndian(File, this->Dimension);
  if (this->Dimension == -1)
    {
      WriteLittleEndian(File, this->LargeDimension);
      for (long i = 0; i < this->LargeDimension; ++i)
	{
	  WriteLittleEndian(File, this->Numerators[i]);
	  WriteLittleEndian(File, this->Denominators[i]);
	}
    }
  else
    for (int i = 0; i < this->Dimension; ++i)
      {
	WriteLittleEndian(File, this->Numerators[i]);
	WriteLittleEndian(File, this->Denominators[i]);
      }
  File.close();
  this->Delocalize();
  return true;
}

// write vector in a file in ascii mode
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool RationalVector::WriteAsciiVector (const char* fileName)
{
  this->Localize();
  ofstream File;
  File.precision(14);
  File.open(fileName, ios::binary | ios::out);
  if (this->Dimension == -1)
    {
      long ReducedDimension = this->LargeDimension - 1;
      for (long i = 0; i < ReducedDimension; ++i)
	File << this->Numerators[i] << " " << this->Denominators[i] << "  ";
      File << this->Numerators[ReducedDimension] << " " << this->Denominators[ReducedDimension] << endl;  
    }
  else
    {
      int ReducedDimension = this->LargeDimension - 1;
      for (int i = 0; i < ReducedDimension; ++i)
	File << this->Numerators[i] << " " << this->Denominators[i] << "  ";
      File << this->Numerators[ReducedDimension] << " " << this->Denominators[ReducedDimension] << endl;  
    }
  File.close();
  this->Delocalize();
  return true;
}

// read vector from a file 
//
// fileName = name of the file where the vector has to be read
// return value = true if no error occurs

bool RationalVector::ReadVector (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  
  std::streampos ZeroPos, MaxPos;
  File.seekg (0, ios::beg);
  ZeroPos = File.tellg();
  File.seekg (0, ios::end);
  MaxPos = File.tellg ();

  std::streampos Length = MaxPos-ZeroPos-sizeof(int);  
  File.seekg (0, ios::beg);
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);

  if (TmpDimension > 0)
    {
      if (Length / (2.0 * sizeof(long)) > ((unsigned) TmpDimension))
	{      
	  cout << "Error reading real vector "<<fileName<<": estimated length "<< (Length / (2l * sizeof(long))) <<" vs dimension "<<TmpDimension<<endl;
	  if ((unsigned)TmpDimension == Length/sizeof(double))
	    cout << "This could be a real vector!"<<endl;
	  exit(1);
	}
      this->Resize(TmpDimension);
      for (int i = 0; i < this->Dimension; ++i)
	{
	  ReadLittleEndian(File, this->Numerators[i]);
	  ReadLittleEndian(File, this->Denominators[i]);
	}
    }
  else
    {
      long TmpLargeDimension;
      ReadLittleEndian(File, TmpLargeDimension);
      this->Resize(TmpLargeDimension);
      for (long i = 0; i < this->LargeDimension; ++i)
	{
	  ReadLittleEndian(File, this->Numerators[i]);
	  ReadLittleEndian(File, this->Denominators[i]);
	}
    }
  File.close();
  return true;
}

// read vector from a file, only within a given range of indices
//
// fileName = name of the file where the vector has to be read
// minIndex = index of the first component to read (if negative, start from the end of vector)
// maxIndex = index of the last component to read (negative or zero is considered as the last component)
// return value = true if no error occurs

bool RationalVector::ReadVector (const char* fileName, long minIndex, long maxIndex)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  
  unsigned ZeroPos, MaxPos;
  File.seekg (0, ios::beg);
  ZeroPos = File.tellg();
  File.seekg (0, ios::end);
  MaxPos = File.tellg ();

  unsigned Length = MaxPos-ZeroPos-sizeof(int);  
  File.seekg (0, ios::beg);
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);

  if (TmpDimension > 0)
    {
      if (Length / (2l * sizeof(long)) > (unsigned)TmpDimension)
	{      
	  cout << "Error reading real vector "<<fileName<<": estimated length "<< (Length / (2l * sizeof(long))) <<" vs dimension "<<TmpDimension<<endl;
	  if ((unsigned)TmpDimension == Length/sizeof(double))
	    cout << "This could be a real vector!"<<endl;
	  exit(1);
	}
      if ((maxIndex >= TmpDimension) || (maxIndex <= 0l))
	maxIndex = TmpDimension - 1;
      if (minIndex < 0l)      
	minIndex += TmpDimension;
      this->Resize(maxIndex - minIndex + 1l);
      File.seekg (minIndex * sizeof(double), ios::cur);    
      for (int i = 0; i < this->Dimension; ++i)
	{
	  ReadLittleEndian(File, this->Numerators[i]);
	  ReadLittleEndian(File, this->Denominators[i]);
	}
    }
  else
    {
      long TmpLargeDimension;
      ReadLittleEndian(File, TmpLargeDimension);
      if ((maxIndex >= TmpLargeDimension) || (maxIndex <= 0l))
	maxIndex = TmpLargeDimension - 1l;
      if (minIndex < 0l)      
	minIndex += TmpLargeDimension;
      this->Resize(maxIndex - minIndex + 1l);
      File.seekg (minIndex * sizeof(double), ios::cur);    
      for (long i = 0; i < this->LargeDimension; ++i)
	{
	  ReadLittleEndian(File, this->Numerators[i]);
	  ReadLittleEndian(File, this->Denominators[i]);
	}
    }
  File.close();
  return true;
}

// read vector dimension from a file, without loading the full vector 
//
// fileName = name of the file where the vector has to be read
// return value = vector dimension

long RationalVector::ReadVectorDimension (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);
  if (TmpDimension > 0)
    {
      File.close();
      return ((long) TmpDimension);
    }
  long TmpLargeDimension = 0l;
  ReadLittleEndian(File, TmpLargeDimension);
  File.close();
  return TmpLargeDimension;
}

// test if a vector can be read from a file (matching the right type), without loading the full vector 
//
// fileName = name of the file where the vector has to be read
// return value = true if the vector can be read

bool RationalVector::ReadVectorTest (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  
  std::streampos ZeroPos, MaxPos;
  File.seekg (0, ios::beg);
  ZeroPos = File.tellg();
  File.seekg (0, ios::end);
  MaxPos = File.tellg ();

  long Length = (long) (MaxPos - ZeroPos) - sizeof(int);  
  File.seekg (0, ios::beg);
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);

  if (TmpDimension > 0)
    {
      File.close();
      Length /= 2l * sizeof(long);
      if (Length == ((long) TmpDimension))
	return true;
      else
	return false;
    }
  else
    {
      long TmpLargeDimension;
      ReadLittleEndian(File, TmpLargeDimension);
      File.close();
      Length -= sizeof (long);
      Length /= 2l * sizeof(long);
      if (Length == TmpLargeDimension)
	return true;
      else
	return false;
    }
  File.close();
  return false;
}

// Output Stream overload
//
// str = reference on output stream
// v = vector to print
// return value = referenceint GetNumSites(){return this->NSites;} on output stream

ostream& operator << (ostream& str, RationalVector& v)
{
  for (long i = 0; i < v.LargeDimension; ++i)
    {
      if (v.Denominators[i] == 1l)
	str << v.Numerators[i] << endl;
      else
	str << v.Numerators[i] << "/" << v.Denominators[i] << endl;
    }
  return str;
}
