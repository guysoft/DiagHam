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


#ifndef OPTIONGROUP_H
#define OPTIONGROUP_H


#include "config.h"
#include "GeneralTools/List.h"

#include <iostream>


using std::ostream;


class AbstractOption;


class OptionGroup
{

 protected:
  
  // ordered list containing all options belonging to the group
  List<AbstractOption*> Options;

  // group full name
  char* GroupName;

 public:

  // constructor
  //
  // groupName = group full name
  OptionGroup(char* groupName);

  // destructor
  //
  ~OptionGroup();

  // add an option to the group
  // 
  // option = pointer to the option to add 
  // return value = reference on the current option group
  OptionGroup& operator += (AbstractOption* option);

  // get option from its name
  //
  // optionName = string containing option name
  // return value = poitner to the option if it has been found, 0 either
  AbstractOption* operator[] (char* optionName);

  // Test if an argument corresponds to the current option and read its content
  //
  // argumentValues = string array of arguments
  // nbrArgument = number of arguments in argumentValues array
  // argumentPosition = position of the first argument to read
  // return value = number of arguments that have been read (-1 if an error occured)
  int ReadOption(char** argumentValues, int nbrArgument, int argumentPosition);

  // print error message on output stream
  //
  // output = reference on output stream;
  // return value = reference on current output stream
  ostream& PrintError (ostream& output);

  // print help concerning current option group
  //
  // output = reference on output stream;
  // return value = reference on current output stream
  ostream& DisplayHelp (ostream& output);

};

#endif
