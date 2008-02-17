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

#include <iostream>

using std::cout;
using std::endl;


// standard constructor
ComplexUniqueArray::ComplexUniqueArray(int internalSize)
{
  this->InternalSize=internalSize;
  this->NbrElements=0;
  if (internalSize>0)
    {
      this->Elements=new Complex[internalSize];
      this->Flag.Initialize();
    }
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
	  for (int i=0; i<NbrElements; ++i)
	    this->Elements[i]=array.Elements[i];
	  this->Flag.Initialize();
	}
    }
}

// destructor
ComplexUniqueArray::~ComplexUniqueArray()
{
  if ( (this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] Elements;
    }
}

// Insert element
// element = new element to be inserted
// returns : index of this element  
int ComplexUniqueArray::InsertElement(const Complex& element)
{
  for (int i=0; i<NbrElements; ++i)
    {
      if (Elements[i]==element)
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
      this->InternalSize*=2;      
      Complex *newElements= new Complex[InternalSize];
      for (int i=0; i<NbrElements; ++i)
	newElements[i]=Elements[i];
      newElements[NbrElements]=element;
      ++NbrElements;
      if ((this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
	delete [] Elements;
      this->Elements=newElements;      
    }
  return (NbrElements-1);
}

// search entry
// value = value to be searched for
// returns : index of the element, or -1 if not found
int ComplexUniqueArray::SearchElement(const Complex &value)
{
  for (int i=0; i<NbrElements; ++i)
    {
      if (Elements[i]==value)
	return i;
    }
  // element not found
  return -1;
}

