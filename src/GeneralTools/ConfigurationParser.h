////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class that provide parser for configuration files             //
//                                                                            //
//                        last modification : 06/01/2005                      //
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


#ifndef CONFIGURATIONPARSER_H
#define CONFIGURATIONPARSER_H

#include "config.h"
#include "GeneralTools/List.h"

#include <iostream>


using std::ostream;


class ConfigurationParser
{

 protected:

  // list of all parameter names
  List<char*> ParameterNames;
  // list of all parameter values (same order than the ParameterNames list)
  List<char*> ParameterValues;

  // list of all error messages 
  List<char*> ErrorLog;

 public:

  // default constructor
  //
  ConfigurationParser();

  // destructor
  //
  ~ConfigurationParser();

  // parse configuration from a file
  //
  // filename = name of the file to parse
  // return value = true if no error occurs
  bool Parse(char* filename);

  // get the string corresponding to a configuration parameter 
  //
  // parameterName = string corresponding to a parameter name
  // return value = string (must not be altered) corresponding to a configuration parameter, null if the parameter is not defined
  char* operator [] (char* parameterName);

  // dump configuration file
  //
  // str = reference on the output stream
  // return value = reference on the output stream
  ostream& DumpConfiguration (ostream& str);

  // print last error encountered during parsing operation
  //
  // str = reference on the output stream
  // return value = reference on the output stream
  ostream& PrintLastError (ostream& str);

  // dump all errors encountered during parsing operation
  //
  // str = reference on the output stream
  // return value = reference on the output stream
  ostream& DumpErrors (ostream& str);

 private:

  // clean a line from useless comments and spaces
  // 
  // line = pointer to the line to clean
  // return value = true if the line still contains usefull information
  bool CleanLine (char* line);

};

#endif
