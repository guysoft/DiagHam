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

#ifndef SORTEDCOMPLEXUNIQUEARRAY_H
#define SORTEDCOMPLEXUNIQUEARRAY_H

#include "config.h"

#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"

#include <iostream>

using std::ostream;
using std::ofstream;
using std::ifstream;
class SortedComplexUniqueArray;

class SortedComplexUniqueArray
{
 protected:
  // array with elements
  Complex *Elements;
  // tolerance for taking elements to be the same (squared)
  double ToleranceSqr;
  // size of array
  unsigned InternalSize;
  // number of elements stored
  unsigned NbrElements;

  // garbage flag
  GarbageFlag Flag;

  // flag indicating how many entries have been sorted
  unsigned Sorted;

  // flag indicating whether to keep elements sorted
  bool KeepSorted;

 public:
  // standard constructor
  SortedComplexUniqueArray(double tolerance = MACHINE_PRECISION, unsigned internalSize=128, bool keepSorted=true);

  // copy constructor
  SortedComplexUniqueArray(SortedComplexUniqueArray &array, bool duplicateFlag = false);

  // destructor
  ~SortedComplexUniqueArray();

  // Insert element
  // element = new element to be inserted
  // returns : index of this element  
  unsigned InsertElement(const Complex &element);

  // search entry using a binary search
  // value = value to be searched for
  // @param[out] index : index of the element, or -1 if not found
  // return : true if element was found, false otherwise.
  bool SearchElement(const Complex &value, unsigned &index);

  // search entry, performing a linear search
  // value = value to be searched for
  // @param[out] index : index of the element, or -1 if not found
  // @param enhanceTolerance : increase tolerance by this factor
  // return : true if element was found, false otherwise.
  bool CarefulSearchElement(const Complex &value, unsigned &index, double enhanceTolerance=1.0);

  // search entry closest to a given value
  // value = value to be searched for
  // @param[out] index : index of a nearby element, or the element itself, if found.
  // return : true if the exact element was found, false otherwise.
  bool NearbyEntry(const Complex &value, unsigned &index);

  // get number of elements
  unsigned GetNbrElements(){ return NbrElements;}

  // empty all elements
  // disallocate = flag indicating whether all memory should be unallocated
  // internalSize = minimum table size to allocate (only used if disallocating)
  void Empty(bool disallocate = false, unsigned internalSize = 100);

  // Access an element
  Complex& operator [] (unsigned i);

  // Sort the entries
  void SortEntries();
  
  // Test if the array is sorted
  bool IsSorted();
  
  // Merge data with another UniqueArray
  void MergeArray(SortedComplexUniqueArray &a);

  // Test all entries
  // search for entries and make sure their indices match the search result
  // result: true if all entries are found, false otherwise
  bool TestAllEntries();

  // Write to file
  // file = open stream to write to
  void WriteArray(ofstream &file);

  // Read from file
  // file = open stream to read from
  void ReadArray(ifstream &file);
   
  // Test object
  static void TestClass(unsigned samples=2048, bool keepSorted=true);

  // output stream overload
  friend ostream& operator << (ostream& Str, const SortedComplexUniqueArray &A);

};


// return array's i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline Complex& SortedComplexUniqueArray::operator [] (unsigned i)
{
  return this->Elements[i];
}

#endif
