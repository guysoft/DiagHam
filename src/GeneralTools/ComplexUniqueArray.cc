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


#include "ComplexUniqueArray.h"

#include "GeneralTools/Endian.h"

#include <iostream>

using std::cout;
using std::endl;


// standard constructor
ComplexUniqueArray::ComplexUniqueArray(unsigned internalSize)
{
  this->InternalSize=internalSize;
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
}

ComplexUniqueArray::ComplexUniqueArray(ComplexUniqueArray &array, bool duplicateFlag)
{
  this->InternalSize=array.InternalSize;
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
#ifdef __SMP__
  this->BufferMutex = new pthread_mutex_t;
  pthread_mutex_init(this->BufferMutex, NULL);
#endif
}

// destructor
ComplexUniqueArray::~ComplexUniqueArray()
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
unsigned ComplexUniqueArray::InsertElement(const Complex& element)
{
  for (unsigned i=0; i<NbrElements; ++i)
    {
      if (SqrNorm(Elements[i]-element)<1e-30)
	return i;
    }
#ifdef __SMP__
  pthread_mutex_lock(this->BufferMutex);
#endif
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
	      cout << "Array overflow in ComplexUniqueArray: cannot store more entries"<<endl;
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
      if ((this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
	delete [] Elements;
      this->Elements=newElements;      
    }
  unsigned Result=NbrElements-1;
#ifdef __SMP__
  pthread_mutex_unlock(this->BufferMutex);
#endif
  return Result;
}

// search entry
// value = value to be searched for
// returns : index of the element, or -1 if not found
unsigned ComplexUniqueArray::SearchElement(const Complex &value)
{
  for (unsigned i=0; i<NbrElements; ++i)
    {
      if (SqrNorm(Elements[i]-value)<1e-30)
	return i;
    }
  // element not found
  return -1;
}

// empty all elements
// disallocate = flag indicating whether all memory should be unallocated
// internalSize = minimum table size to allocate (only used if disallocating)
void ComplexUniqueArray::Empty(bool disallocate, unsigned internalSize)
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


// write to file
// file = open stream to write to
void ComplexUniqueArray::WriteArray(ofstream &file)
{
  WriteLittleEndian(file, this->NbrElements);
  for (unsigned i = 0; i < this->NbrElements; ++i)
    WriteLittleEndian(file, this->Elements[i]);  
}

// Read from file
// file = open stream to read from
void ComplexUniqueArray::ReadArray(ifstream &file)
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
