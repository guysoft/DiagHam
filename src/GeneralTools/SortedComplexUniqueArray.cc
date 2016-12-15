////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Nicolas Regnault                  //
//                                                                            //
//                         class author: Gunnar Moeller                       //
//                                                                            //
//                 class for an array which has unique entries                //
//                                                                            //
//                        last modification : 13/02/2008                      //
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


#include "SortedComplexUniqueArray.h"

#include "GeneralTools/Endian.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"

#include <iostream>
#include <limits>
#include <cstdlib>
#include <cmath>

using std::cout;
using std::endl;


// standard constructor
SortedComplexUniqueArray::SortedComplexUniqueArray(unsigned internalSize, double tolerance)
{
  this->InternalSize=internalSize;
  this->Tolerance=tolerance;
  this->ToleranceSqr=tolerance*tolerance;
  this->NbrElements=0;
  if (internalSize>0)
    {
      this->Elements=new Complex[internalSize];
      this->Flag.Initialize();
    }
#ifdef __SMP__
  this->BufferMutex = new pthread_mutex_t;
  pthread_mutex_init(this->BufferMutex, NULL);
#endif
  this->Sorted = true;
}

SortedComplexUniqueArray::SortedComplexUniqueArray(SortedComplexUniqueArray &array, bool duplicateFlag)
{
  this->InternalSize=array.InternalSize;
  this->Tolerance=array.Tolerance;
  this->ToleranceSqr=array.ToleranceSqr;
  this->NbrElements=array.NbrElements;
  if (duplicateFlag == false)
    {
      this->Elements = array.Elements;
      this->Flag = array.Flag;
    }
  else
    {
      if (this->InternalSize>0)
	{
	  this->Elements=new Complex[InternalSize];
	  for (unsigned i=0; i<NbrElements; ++i)
	    this->Elements[i]=array.Elements[i];
	  this->Flag.Initialize();
	}
    }
  this->Sorted = array.Sorted;
#ifdef __SMP__
  this->BufferMutex = new pthread_mutex_t;
  pthread_mutex_init(this->BufferMutex, NULL);
#endif
}

// destructor
SortedComplexUniqueArray::~SortedComplexUniqueArray()
{
  if ( (this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] Elements;
    }
#ifdef __SMP__
  pthread_mutex_destroy(this->BufferMutex);
  delete this->BufferMutex;
#endif
}

// Insert element
// element = new element to be inserted
// returns : index of this element  
unsigned SortedComplexUniqueArray::InsertElement(const Complex& element)
{
  this->Sorted=false;
#ifdef __SMP__
  pthread_mutex_lock(this->BufferMutex);
#endif
  unsigned TmpNbrElements = this->NbrElements; // safely get the current number of elements
#ifdef __SMP__
  pthread_mutex_unlock(this->BufferMutex);
#endif
  // start searching without locking memory access
  unsigned i=0;
  for (; i<TmpNbrElements; ++i)
    {
      if (SqrNorm(Elements[i]-element)<this->ToleranceSqr)
	return i;
    }
  // element not found among previously existing entries, but another thread may have added further elements in the meantime.
#ifdef __SMP__
  pthread_mutex_lock(this->BufferMutex);
#endif
  for (; i<this->NbrElements; ++i)
    {
      if (SqrNorm(Elements[i]-element)<this->ToleranceSqr)
	return i;
    }
  // element not found
  if (NbrElements < InternalSize)
    {
      this->Elements[NbrElements]=element;
      ++NbrElements;      
    }
  else
    {
      if (this->InternalSize < (std::numeric_limits<unsigned>::max() / 2))
	this->InternalSize*=2;
      else
	{
	  if (this->InternalSize == std::numeric_limits<unsigned>::max())
	    {
	      cout << "Array overflow in SortedComplexUniqueArray: cannot store more entries"<<endl;
	      exit(1);
	    }
	  else
	    this->InternalSize = std::numeric_limits<unsigned>::max();
	}
      Complex *newElements= new Complex[InternalSize];
      for (unsigned i=0; i<NbrElements; ++i)
	newElements[i]=Elements[i];
      newElements[NbrElements]=element;
      ++NbrElements;
      Complex *tmpElements=this->Elements;
      this->Elements=newElements;
      if ((this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
	delete [] tmpElements;
    }
  unsigned Result=NbrElements-1;
#ifdef __SMP__
  pthread_mutex_unlock(this->BufferMutex);
#endif
  return Result;
}

// search entry
// value = value to be searched for
// @param[out] index : index of the element, if found.
// return : true if element was found, false otherwise.
bool SortedComplexUniqueArray::SearchElement(const Complex &value, unsigned &index)
{
  if (this->Sorted)
    {
      unsigned PosMax = this->NbrElements-1;
      unsigned PosMin = 0;
      unsigned PosMid = (PosMin + PosMax) >> 1;
      Complex CurrentState = this->Elements[PosMid];
      cout << "Searching "<<value<<"...";
      while ((PosMin != PosMid) && (SqrNorm(CurrentState - value) >= this->ToleranceSqr))
	{
	  if (CurrentState > value)
	    {
	      PosMax = PosMid;
	    }
	  else
	    {
	      PosMin = PosMid;
	    } 
	  PosMid = (PosMin + PosMax) >> 1;
	  CurrentState = this->Elements[PosMid];
	  //cout << "PosMid="<<PosMid<<", CurrentState="<<CurrentState<<", PosMin="<<PosMin<<", PosMax="<<PosMax<< endl;
	}
      if (SqrNorm(CurrentState - value) < this->ToleranceSqr)
	{
	  index = PosMid;
	  cout << "Found "<<this->Elements[index]<<endl;
	  return true;
	}
      else
	{
	  index = PosMax; // ?? 
	  if (this->Elements[PosMax] == value)
	    {
	      cout << "Found "<<this->Elements[PosMax]<<endl;
	      return true;
	    }
	  else 
	    {
	      cout << "Not found."<<endl;
	      return false;
	    }
	}
    }
  else for (unsigned i=0; i<NbrElements; ++i)
    {
      if (SqrNorm(Elements[i]-value)<this->ToleranceSqr)
	{
	  index = i;
	  return true;
	}
    }
  // element not found
  index = 0;
  return false;
}

// empty all elements
// disallocate = flag indicating whether all memory should be unallocated
// internalSize = minimum table size to allocate (only used if disallocating)
void SortedComplexUniqueArray::Empty(bool disallocate, unsigned internalSize)
{
#ifdef __SMP__
  pthread_mutex_lock(this->BufferMutex);
#endif
  if (disallocate)
    {
      if ((this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
	{
	  delete [] Elements;
	}
      this->InternalSize=internalSize;
      this->Elements = new Complex[InternalSize];
    }
  this->NbrElements=0;
#ifdef __SMP__
  pthread_mutex_unlock(this->BufferMutex);
#endif
}

// apply a simple sort algorithm to the existing entries of the table
void SortedComplexUniqueArray::SortEntries()
{
  if (this->Sorted==true) return;
  unsigned inc = std::round(NbrElements/2.0);
  Complex tmpC;
  while (inc > 0)
    {
      for (unsigned i = inc; i< NbrElements; ++i)
	{
	  tmpC = this->Elements[i];
	  unsigned j = i;
	  while ((j>=inc) && (this->Elements[j-inc] > tmpC) )
	    {
	      this->Elements[j] = this->Elements[j - inc];
	      j = j - inc;
	    }
	  this->Elements[j] = tmpC;
	}
      inc = std::round(inc / 2.2);
    }
  this->Sorted=true;
}

// Merge data with another UniqueArray
void SortedComplexUniqueArray::MergeArray(SortedComplexUniqueArray &a)
{
  a.SortEntries();
  this->SortEntries();

  cout << "this ="<<*this<<endl;
  cout << "a ="<<a << endl;
  cout << "Merging arrays with "<<this->NbrElements<<" and "<<a.NbrElements<<" entries "<<endl;
  long newPos = a.NbrElements+this->NbrElements;
  if (newPos > std::numeric_limits<unsigned>::max())
    {
      cout << "Error merged array size exceeds maximum"<< endl;
      exit(1);
    }
  Complex *newElements = new Complex[newPos];
  newPos=0;
  unsigned myPos = 0;
  for (unsigned theirPos=0; theirPos<a.NbrElements; ++theirPos)
    {
      // definite insert elements in this-> that are smaller than the next element of a and not within tolerance
      while (myPos < this->NbrElements && this->Elements[myPos] < a.Elements[theirPos] && SqrNorm(this->Elements[myPos] - a.Elements[theirPos]) >= this->ToleranceSqr)
	{
	  cout << "Insert this->Elements["<<myPos<<"] at 1 with this->Elements["<<myPos<<"] < a.Elements["<<theirPos<<"]"<<endl;
	  newElements[newPos++] = this->Elements[myPos++];
	}
      // also insert any elements that are equal or approximately equal
      while (SqrNorm(this->Elements[myPos] - a.Elements[theirPos]) < this->ToleranceSqr)
	{
	  cout << "Insert this->Elements["<<myPos<<"] at 2"<<endl;
	  newElements[newPos++] = this->Elements[myPos++];
	}
      cout << "myPos = "<<myPos << " this->NbrElements = "<<this->NbrElements << " this->Elements[myPos]="<<this->Elements[myPos]<<", a.Elements[theirPos] ="<<a.Elements[theirPos]<<endl;
      if (newPos > 0)
	{
	  if (SqrNorm(a.Elements[theirPos] - newElements[newPos-1]) < this->ToleranceSqr)
	    {
	      // next element on a is identical within accuracy to the last entry on the new array
	      cout << "Omitting identical element a.Elements["<<theirPos<<"] at 1"<<endl;
	      // do nothing - counter is increased in loop.
	    }
	  else
	    {
	      cout << "Insert a.Elements["<<theirPos<<"] at 1"<<endl;
	      newElements[newPos++] = a.Elements[theirPos];
	    }
	}
      else
	{
	  cout << "Insert a.Elements["<<theirPos<<"] at 2"<<endl;
	  newElements[newPos++] = a.Elements[theirPos];
	}
    }
  
  while (myPos < this->NbrElements) //  && this->Elements[myPos] < a.Elements[theirPos] && SqrNorm(this->Elements[myPos] - a.Elements[theirPos]) >= this->ToleranceSqr)
    {
      cout << "Insert this->Elements["<<myPos<<"] at 3 with this->Elements["<<myPos<<"]"<<endl;
      newElements[newPos++] = this->Elements[myPos++];
    }

  if ( (this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] Elements;
    }
  this->Elements = newElements;
  this->NbrElements = newPos;
  cout << "Unique entries retained: "<<this->NbrElements<<endl;
  this->Sorted=true;
}


// write to file
// file = open stream to write to
void SortedComplexUniqueArray::WriteArray(ofstream &file)
{
  WriteLittleEndian(file, this->NbrElements);
  for (unsigned i = 0; i < this->NbrElements; ++i)
    WriteLittleEndian(file, this->Elements[i]);  
}

// Read from file
// file = open stream to read from
void SortedComplexUniqueArray::ReadArray(ifstream &file)
{
#ifdef __SMP__
  pthread_mutex_lock(this->BufferMutex);
#endif
  if ( (this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] Elements;
    }
  unsigned TmpDimension;
  ReadLittleEndian(file, TmpDimension);
  this->InternalSize=TmpDimension;
  this->NbrElements=TmpDimension;
  this->Elements=new Complex[TmpDimension];
  for (unsigned i = 0; i < this->NbrElements; ++i)
    ReadLittleEndian(file, this->Elements[i]);
#ifdef __SMP__
  pthread_mutex_unlock(this->BufferMutex);
#endif

}

// Test object
void SortedComplexUniqueArray::TestClass(unsigned samples)
{
  SortedComplexUniqueArray a1(samples>>1, 1e-13);
  SortedComplexUniqueArray a2(samples>>1, 1e-13);


  NumRecRandomGenerator gen;
  
  // insert half the elements as independent numbers
  for (int i=0; i<samples; ++i)
    {
      a1.InsertElement( Complex(gen.GetRealRandomNumber(),gen.GetRealRandomNumber()) );
      a2.InsertElement( Complex(gen.GetRealRandomNumber(),gen.GetRealRandomNumber()) );
    }
  // insert other half as same number
  Complex TmpC;
  for (int i=0; i<samples; ++i)
    {
      TmpC.Re=gen.GetRealRandomNumber();
      TmpC.Im=gen.GetRealRandomNumber();
      a1.InsertElement(TmpC);
      a2.InsertElement(TmpC);
    }
  
  SortedComplexUniqueArray a3(a1, true);
  SortedComplexUniqueArray a4(a2, true);

  unsigned common=0, index;
  for (int i=0; i<samples; ++i)
    {
      if (a2.SearchElement(a1[i],index))
	{
	  ++common;
	  cout << "Random common entry at index "<<index<<endl;
	}
    }
  
  a3.MergeArray(a2);

  cout << "Merged array has "<< a3.GetNbrElements() <<" entries"<<endl;
  
  if (a3.GetNbrElements() < 3*samples - common)
    {
      cout << "Unexpected number of entries in a3: "<< a3.GetNbrElements() << " vs " << 3*samples - common << " expected "<<endl;
    }

  for (int i=0; i<2*samples; ++i)
    {
      a4.InsertElement(a1[i]);
    }

  cout << "Manually merged array has "<< a4.GetNbrElements() <<" entries"<<endl;
  
  if (a4.GetNbrElements() < 3*samples - common)
    {
      cout << "Unexpected number of entries in a4: "<< a4.GetNbrElements() << " vs " << 3*samples - common << " expected "<<endl;
    }

  // check all entries are present:
  for (unsigned i=0; i<a4.GetNbrElements(); ++i)
    {
      if (! a3.SearchElement(a4[i], index))
	{
	  cout << "Element " << i << " of array 4, "<<a4[i] << ", missing from array 3"<<endl;
	}
    }

  for (unsigned i=0; i<a3.GetNbrElements(); ++i)
    {
      if (! a4.SearchElement(a3[i], index))
	{
	  cout << "Element "<<i<<" of array 3, "<<a3[i] << ", missing from array 4"<<endl;
	}
    }

  // check all entries are found:
  for (unsigned i=0; i<a3.GetNbrElements(); ++i)
    {
      if (! a3.SearchElement(a3[i], index))
	{
	  cout << "Element " << i << " not found by search on array 3"<<endl;
	}
      if (index != i)
	cout << "Discrepancy in search for index "<< index <<" on element " << i << " (array 3)."<<endl;
    }


    for (unsigned i=0; i<a4.GetNbrElements(); ++i)
    {
      if (! a4.SearchElement(a4[i], index))
	{
	  cout << "Element " << i << " not found by search on array 4"<<endl;
	}
      if (index != i)
	cout << "Discrepancy in search for index "<< index <<" on element " << i << " (array 4)."<<endl;
    }

}


// Output Stream overload
//

ostream& operator << (ostream& Str, const SortedComplexUniqueArray& A)
{
  for (unsigned i = 0; i < A.NbrElements; ++i)
    Str << A.Elements[i]<<endl;
  return Str;
}

