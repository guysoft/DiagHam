////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      various generic tools about arrays                    //
//                                                                            //
//                        last modification : 15/09/2003                      //
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

#ifndef ARRAYTOOLS_H
#define ARRAYTOOLS_H

#include "config.h"

// up ordering array sort using quick sort
//
// array = pointer to the array
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayUpOrdering(ClassName* array, int nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] > array[1])
	  {
	    ClassName TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	  }
	if (array[1] > array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	  }	
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ClassName TmpElement;
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	  }
	if (array[i] >  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	  }
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	  }
	--j;
	ClassName Pivot = array[i];
	array[i] = array[j];
	array[j] = Pivot;
	i = 0;
	while (true)
	  {
	    while (array[++i] < Pivot);
	    while (array[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	SortArrayUpOrdering(array, i);
	SortArrayUpOrdering(&(array[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}


#endif
