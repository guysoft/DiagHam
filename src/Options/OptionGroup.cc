////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of group of options                         //
//                                                                            //
//                        last modification : 19/05/2004                      //
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
#include "GeneralTools/ListIterator.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"

#include <iostream>
#include <string.h>


using std::ostream;
using std::endl;


// constructor
//
// groupName = group full name

OptionGroup::OptionGroup(char* groupName)
{
  this->GroupName = new char [strlen(groupName) + 1];
  strcpy (this->GroupName, groupName);
}

// destructor
//

OptionGroup::~OptionGroup()
{
  delete[] this->GroupName;
  ListIterator<AbstractOption*> IterOption(this->Options);
  AbstractOption** TmpOption;
  while ((TmpOption = IterOption()))
    {
      delete *TmpOption;
    }
}

// add an option to the group
// 
// option = pointer to the option to add 
// return value = reference on the current option group

OptionGroup& OptionGroup::operator += (AbstractOption* option)
{
  this->Options += option;
  return *this;
}

// get option from its name
//
// optionName = string containing option name
// return value = poitner to the option if it has been found, 0 either

AbstractOption* OptionGroup::operator[] (char* optionName)
{
  ListIterator<AbstractOption*> IterOption(this->Options);
  AbstractOption** TmpOption;
  while ((TmpOption = IterOption()))
    {
      if ((*TmpOption)->IsOptionName(optionName) == true)
	{
	  return *TmpOption;
	}
    }
  return 0;
}

// Test if an argument corresponds to the current option and read its content
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// argumentPosition = position of the first argument to read
// return value = number of arguments that have been read (-1 if an error occured)

int OptionGroup::ReadOption(char** argumentValues, int nbrArgument, int argumentPosition)
{
  int Pos = 0;
  ListIterator<AbstractOption*> IterOption(this->Options);
  AbstractOption** TmpOption;
  while ((Pos == 0) && ((TmpOption = IterOption())))
    {
      Pos = (*TmpOption)->ReadOption(argumentValues, nbrArgument, argumentPosition);
    }
  return Pos;
}

// print error message on output stream
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& OptionGroup::PrintError (ostream& output)
{
  ListIterator<AbstractOption*> IterOption(this->Options);
  AbstractOption** TmpOption;
  while ((TmpOption = IterOption()))
    {
      (*TmpOption)->PrintError(output) << endl;
    }
  return output;
}

// print help concerning current option group
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& OptionGroup::DisplayHelp (ostream& output)
{
  output << this->GroupName << ":" << endl;
  ListIterator<AbstractOption*> IterOption(this->Options);
  AbstractOption** TmpOption;
  while ((TmpOption = IterOption()))
    {
      output << "    ";
      (*TmpOption)->DisplayHelp(output) << endl;
    }
  return output;
}
