////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of option manager                         //
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


#ifndef OPTIONMANAGER_H
#define OPTIONMANAGER_H


#include "config.h"
#include "GeneralTools/List.h"

#include <iostream>


using std::ostream;


class OptionGroup;
class AbstractOption;

class OptionManager
{

 protected:
  
  // ordered list containing all option groups which have to be handled by the option manager
  List<OptionGroup*> Groups;

  // program name 
  char* ProgramName;
 
  // string containing program version
  char* ProgramVersion;

  // an optional string that can be displayed with help
  char* ProgramAdditionalInformations;

 public:

  // constructor
  //
  // groupName = group full name
  OptionManager(char* programName, char* programVersion = 0, char* programAdditionalInformations = 0);

  // destructor
  //
  ~OptionManager();

  // add an option group to the manager
  // 
  // group = pointer to the option group to add 
  // return value = reference on the current option manager
  OptionManager& operator += (OptionGroup* group);

  // get option from its name
  //
  // optionName = string containing option name
  // return value = poitner to the option if it has been found, 0 either
  AbstractOption* operator[] (char* optionName);

  // Proceed running options from command line arguments
  //
  // argumentValues = string array of arguments
  // nbrArgument = number of arguments in argumentValues array
  // output = reference on output stream used to display errors
  // return value = true if proceeding succeded, false if an error occurs
  bool ProceedOptions (char** argumentValues, int nbrArgument, ostream& output);

  // print help concerning current option group
  //
  // output = reference on output stream;
  // return value = reference on current output stream
  ostream& DisplayHelp (ostream& output);

};

#endif
