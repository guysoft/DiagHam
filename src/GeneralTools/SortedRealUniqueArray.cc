////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Nicolas Regnault                  //
//                                                                            //
//                         class author: Gunnar Moeller                       //
//                                                                            //
// class for an array which has unique entries for single-threaded insertion  //
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


#include "SortedRealUniqueArray.h"

#include "GeneralTools/Endian.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"

#include <iostream>
#include <limits>
#include <cstdlib>
#include <cmath>

using std::cout;
using std::endl;
using std::max;

// flag for testing
#define TESTING_SRUA

// standard constructor
SortedRealUniqueArray::SortedRealUniqueArray(unsigned internalSize, double tolerance, bool keepSorted)
{
  this->InternalSize=internalSize;
  this->Tolerance=tolerance;
  this->NbrElements=0;
  if (internalSize>0)
    {
      this->Elements=new double[internalSize];
      this->Flag.Initialize();
    }
  this->Sorted = 0;
  this->KeepSorted = keepSorted;
}

SortedRealUniqueArray::SortedRealUniqueArray(SortedRealUniqueArray &array, bool duplicateFlag)
{
  this->InternalSize=array.InternalSize;
  this->Tolerance=array.Tolerance;
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
	  this->Elements=new double[InternalSize];
	  for (unsigned i=0; i<NbrElements; ++i)
	    this->Elements[i]=array.Elements[i];
	  this->Flag.Initialize();
	}
    }
  this->Sorted = array.Sorted;
  this->KeepSorted = array.KeepSorted;
}

// destructor
SortedRealUniqueArray::~SortedRealUniqueArray()
{
  if ( (this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] Elements;
    }
}

// Insert element
// element = new element to be inserted
// returns : index of this element  
unsigned SortedRealUniqueArray::InsertElement(const double& element)
{
  unsigned index;
  if (this->SearchElement(element, index))
    return index;
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
	      cout << "Array overflow in SortedRealUniqueArray: cannot store more entries"<<endl;
	      exit(1);
	    }
	  else
	    this->InternalSize = std::numeric_limits<unsigned>::max();
	}
      double *newElements= new double[InternalSize];
      for (unsigned i=0; i<NbrElements; ++i)
	newElements[i]=Elements[i];
      newElements[NbrElements]=element;
      ++NbrElements;
      double *tmpElements=this->Elements;
      this->Elements=newElements;
      if ((this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
	delete [] tmpElements;
    }
  unsigned Result=NbrElements-1;
  if (this->KeepSorted && this->NbrElements > this->Sorted + 12)
    {
      this->SortEntries();
      if (!this->SearchElement(element, Result))
	{
	  cout << "Error: did not find the element that was just inserted"<<endl;
	  exit(1);
	}
    }
  return Result;
}

// search entry
// value = value to be searched for
// @param[out] index : index of the element, if found.
// return : true if element was found, false otherwise.
bool SortedRealUniqueArray::SearchElement(const double &value, unsigned &index)
{
  unsigned start=0;
  if (this->Sorted>3)
    {
      unsigned PosMax = this->Sorted - 1;
      unsigned PosMin = 0;
      unsigned PosMid = (PosMin + PosMax) >> 1;
      double CurrentState = this->Elements[PosMid];
      // cout << "Searching "<<value<<"...";
      while ((PosMin != PosMid) && (fabs(CurrentState - value) >= this->Tolerance))
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
      if (fabs(CurrentState - value) < this->Tolerance)
	{
	  index = PosMid;
	  // cout << "Found "<<this->Elements[index]<<endl;
	  return true;
	}
      else
	{
	  index = PosMax; // ?? 
	  if (fabs(this->Elements[PosMax] - value) < this->Tolerance)
	    {
	      // cout << "Found "<<this->Elements[PosMax]<<endl;
	      return true;
	    }
	  else 
	    {
	      // cout << "Not found."<<endl;
	      start = this->Sorted-1;
	    }
	}
    }
  for (unsigned i=start; i<this->NbrElements; ++i)
    {
      if (fabs(Elements[i]-value)<this->Tolerance)
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
void SortedRealUniqueArray::Empty(bool disallocate, unsigned internalSize)
{
  if (disallocate)
    {
      if ((this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
	{
	  delete [] Elements;
	}
      this->InternalSize=internalSize;
      this->Elements = new double[InternalSize];
    }
  this->NbrElements=0;
}

// apply a simple sort algorithm to the existing entries of the table
void SortedRealUniqueArray::SortEntries()
{
  if (this->Sorted==this->NbrElements) return;
  unsigned inc = std::floor(NbrElements/2.0+0.5);
  // if (this->Sorted>inc) inc=this->Sorted-1;
  double tmpC;
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
      inc = std::floor(inc / 2.2 + 0.5);
    }
  this->Sorted=this->NbrElements;

#ifdef TESTING_SRUA
  if (!this->IsSorted())
    {
      cout << "Error: sort algorithm left unsorted array."<<endl;
      exit(1);
    }
#endif //TESTING_SRUA
}


// Test if the array is sorted
bool SortedRealUniqueArray::IsSorted()
{
  if (this->NbrElements<2) return true;
  for (unsigned i=0; i<this->NbrElements-1; ++i)
    if (this->Elements[i]>=this->Elements[i+1])
      {
	cout << "Unexpected order at position "<<i<<": (this->Elements["<<i<<"]<this->Elements["<<i+1<<"]"<<endl; 
	return false;
      }
  return true;
}


// Merge data with another UniqueArray
void SortedRealUniqueArray::MergeArray(SortedRealUniqueArray &a)
{
  a.SortEntries();
  this->SortEntries();

  // cout << "this ="<<*this<<endl;
  // cout << "a ="<<a << endl;
  // cout << "Merging arrays with "<<this->NbrElements<<" and "<<a.NbrElements<<" entries "<<endl;
  long newPos = a.NbrElements+this->NbrElements;
  if (newPos > std::numeric_limits<unsigned>::max())
    {
      cout << "Error merged array size exceeds maximum"<< endl;
      exit(1);
    }
  double *newElements = new double[newPos];
  newPos=0;
  unsigned myPos = 0;
  for (unsigned theirPos=0; theirPos<a.NbrElements; ++theirPos)
    {
      // definite insert elements in this-> that are smaller than the next element of a and not within tolerance
      while (myPos < this->NbrElements && this->Elements[myPos] < a.Elements[theirPos] && fabs(this->Elements[myPos] - a.Elements[theirPos]) >= this->Tolerance)
	{
	  // cout << "Insert this->Elements["<<myPos<<"] at 1 with this->Elements["<<myPos<<"] < a.Elements["<<theirPos<<"]"<<endl;
	  newElements[newPos++] = this->Elements[myPos++];
	}
      // also insert any elements that are equal or approximately equal
      while (fabs(this->Elements[myPos] - a.Elements[theirPos]) < this->Tolerance)
	{
	  //  cout << "Insert this->Elements["<<myPos<<"] at 2"<<endl;
	  newElements[newPos++] = this->Elements[myPos++];
	}
      // cout << "myPos = "<<myPos << " this->NbrElements = "<<this->NbrElements << " this->Elements[myPos]="<<this->Elements[myPos]<<", a.Elements[theirPos] ="<<a.Elements[theirPos]<<endl;
      if (newPos > 0)
	{
	  if (fabs(a.Elements[theirPos] - newElements[newPos-1]) < this->Tolerance)
	    {
	      // next element on a is identical within accuracy to the last entry on the new array
	      // cout << "Omitting identical element a.Elements["<<theirPos<<"] at 1"<<endl;
	      // do nothing - counter is increased in loop.
	    }
	  else
	    {
	      //cout << "Insert a.Elements["<<theirPos<<"] at 1"<<endl;
	      newElements[newPos++] = a.Elements[theirPos];
	    }
	}
      else
	{
	  //cout << "Insert a.Elements["<<theirPos<<"] at 2"<<endl;
	  newElements[newPos++] = a.Elements[theirPos];
	}
    }
  
  while (myPos < this->NbrElements) //  && this->Elements[myPos] < a.Elements[theirPos] && fabs(this->Elements[myPos] - a.Elements[theirPos]) >= this->Tolerance)
    {
      // cout << "Insert this->Elements["<<myPos<<"] at 3 with this->Elements["<<myPos<<"]"<<endl;
      newElements[newPos++] = this->Elements[myPos++];
    }

  if ( (this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] Elements;
    }
  this->Elements = newElements;
  this->NbrElements = newPos;
  // cout << "Unique entries retained: "<<this->NbrElements<<endl;
  this->Sorted=true;
}

// Test all entries
// search for entries and make sure their indices match the search result
// result: true if all entries are found, false otherwise
bool SortedRealUniqueArray::TestAllEntries()
{
  unsigned index;
  bool success=true;
  for (unsigned i=0; i<this->GetNbrElements(); ++i)
    {
      if (! this->SearchElement(this->Elements[i], index))
	{
	  cout << "Element " << i << " not found during self-check"<<endl;
	  success=false;
	}
      if (index != i)
	{
	  cout << "Discrepancy in search for index "<< index <<" on element " << i << " during self-check."<<endl;
	  success=false;
	}
    }
  return success;
}


// write to file
// file = open stream to write to
void SortedRealUniqueArray::WriteArray(ofstream &file)
{
  WriteLittleEndian(file, this->NbrElements);
  for (unsigned i = 0; i < this->NbrElements; ++i)
    WriteLittleEndian(file, this->Elements[i]);  
}

// Read from file
// file = open stream to read from
void SortedRealUniqueArray::ReadArray(ifstream &file)
{
  if ( (this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] Elements;
    }
  unsigned TmpDimension;
  ReadLittleEndian(file, TmpDimension);
  this->InternalSize=TmpDimension;
  this->NbrElements=TmpDimension;
  this->Elements=new double[TmpDimension];
  for (unsigned i = 0; i < this->NbrElements; ++i)
    ReadLittleEndian(file, this->Elements[i]);
}

// Test object
void SortedRealUniqueArray::TestClass(unsigned samples, bool keepSorted)
{
  double precision = 1e-13;
  SortedRealUniqueArray a1(samples>>1, precision, keepSorted);
  SortedRealUniqueArray a2(samples>>1, precision, keepSorted);


  NumRecRandomGenerator gen;
  
  // insert half the elements as independent numbers
  for (int i=0; i<samples; ++i)
    {
      a1.InsertElement( gen.GetRealRandomNumber() );
      a2.InsertElement( gen.GetRealRandomNumber() );
    }
  // count identical entries
  unsigned common=0, index;
  for (int i=0; i<samples; ++i)
    {
      if (a2.SearchElement(a1[i],index))
	{
	  ++common;
	  cout << "Random common entry at index "<<index<<endl;
	}
    }

  // insert other half as same number
  double Tmp;
  for (int i=0; i<samples; ++i)
    {
      Tmp = gen.GetRealRandomNumber();
      a1.InsertElement(Tmp);
      a2.InsertElement(Tmp);
    }
  
  SortedRealUniqueArray a3(a1, true);
  SortedRealUniqueArray a4(a2, true);
  
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

    cout << "Tests completed."<<endl;

}


// Output Stream overload
//

ostream& operator << (ostream& Str, const SortedRealUniqueArray& A)
{
  for (unsigned i = 0; i < A.NbrElements; ++i)
    Str << A.Elements[i]<<endl;
  return Str;
}

