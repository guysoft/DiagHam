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


#include "config.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/ListIterator.h"

#include <string.h>
#include <fstream>


using std::ifstream;
using std::ios;
using std::endl;


// default constructor
//

ConfigurationParser::ConfigurationParser()
{
}

// destructor
//

ConfigurationParser::~ConfigurationParser()
{
  if (this->ParameterNames.GetNbrElement() > 0)
    {
      ListIterator<char*> IterNames(this->ParameterNames);
      ListIterator<char*> IterValues(this->ParameterValues);
      char** TmpName;
      char** TmpValue;
      while ((TmpName = IterNames()))
	{
	  TmpValue = IterValues();
	  delete[] *TmpName;
	  delete[] *TmpValue;
	}
    }
  if (this->ErrorLog.GetNbrElement() > 0)
    {
      ListIterator<char*> IterError(this->ErrorLog);
      char** TmpError;
      while ((TmpError = IterError()))
	{
	  delete[] (*TmpError);
	}
    }
}

// parse configuration from a file
//
// filename = name of the file to parse
// return value = true if no error occurs

bool ConfigurationParser::Parse(char* filename)
{
  ifstream File;
  File.open(filename, ios::binary | ios::in);
  if (!File.is_open())
    {
      char* TmpString = new char [strlen (filename) + 32]; 
      sprintf (TmpString, "cannot open file: %s\n", filename);
      this->ErrorLog += TmpString;
      return false;
    }
  File.seekg(0, ios::end);
  unsigned int Size = File.tellg();
  File.seekg(0, ios::beg);
  char* TmpBuffer = new char [Size + 1];
  if (TmpBuffer == 0)
    {
      char* TmpString = new char [strlen (filename) + 32]; 
      sprintf (TmpString, "%s is too big\n", filename);
      this->ErrorLog += TmpString;
      return false;
    }
  File.read(TmpBuffer, Size);
  TmpBuffer[Size] = '\0';
  File.close();
  unsigned int Pos = 0;
  char* Start = TmpBuffer;
  int LineNumber = 1;
  bool Flag = true;
  while (Pos < Size)
    {
      Start = TmpBuffer + Pos;
      while ((Pos < Size) && (TmpBuffer[Pos] != '\n'))
	{
	  Pos++;
	}
      TmpBuffer[Pos] = '\0';
      if (this->CleanLine (Start) == true)
	{
	  unsigned int Pos2 = 0;
	  while ((Start[Pos2] != '\0') && (Start[Pos2] != '='))
	    Pos2++;
	  if (Start[Pos2] == '\0')
	    {
	      char* TmpString = new char [64]; 
	      sprintf (TmpString, "syntax error at line %d\n", LineNumber);
	      this->ErrorLog += TmpString;
	      Flag = false;
	    }
	  else
	    {
	      unsigned int Pos3 = Pos2;
	      if (Pos2 > 0)
		{
		  --Pos2;
		  while ((Pos2 > 0) && ((Start[Pos2] == ' ') || (Start[Pos2] == '\t')))
		    --Pos2;
		}
	      if ((Pos2 == 0)  && ((Start[Pos2] == ' ') || (Start[Pos2] == '\t') || (Start[Pos2] == '=')))
		{
		  char* TmpString = new char [64]; 
		  sprintf (TmpString, "no parameter name defined at line %d\n", LineNumber);
		  this->ErrorLog += TmpString;
		  Flag = false;
		}
	      else
		{
		  char* TmpName = new char [Pos2 + 2];
		  strncpy (TmpName, Start, Pos2 + 1);
		  TmpName[Pos2 + 1] = '\0';
		  Pos2 = Pos3 + 1;
		  while ((Start[Pos2] != '\0') && ((Start[Pos2] == ' ') || (Start[Pos2] == '\t')))
		    ++Pos2;
		  if (Start[Pos2] == '\0')
		    {
		      char* TmpString = new char [64]; 
		      sprintf (TmpString, "no parameter value defined at line %d\n", LineNumber);
		      this->ErrorLog += TmpString;
		      delete[] TmpName;
		      Flag = false;
		    }
		  else
		    {
		      char* TmpValue = new char [strlen(Start + Pos2) + 1];
		      strcpy (TmpValue, Start + Pos2);
		      this->ParameterNames += TmpName;
		      this->ParameterValues += TmpValue;
		    }
		}
	    }
	}
      Pos++;
      LineNumber++;
    }
  delete[] TmpBuffer;
  return Flag;
}

// get the string corresponding to a configuration parameter 
//
// parameterName = string corresponding to a parameter name
// return value = string (must not be altered) corresponding to a configuration parameter, null if the parameter is not defined

char* ConfigurationParser::operator [] (char* parameterName)
{
  ListIterator<char*> IterNames(this->ParameterNames);
  ListIterator<char*> IterValues(this->ParameterValues);
  char** TmpName;
  char** TmpValue;
  while ((TmpName = IterNames()))
    {
      TmpValue = IterValues();
      if (strcmp (*TmpName, parameterName) == 0)
	{
	  return *TmpValue;
	}
    }
  return 0;
}

// dump configuration file
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& ConfigurationParser::DumpConfiguration (ostream& str)
{
  ListIterator<char*> IterNames(this->ParameterNames);
  ListIterator<char*> IterValues(this->ParameterValues);
  char** TmpName;
  char** TmpValue;
  while ((TmpName = IterNames()))
    {
      TmpValue = IterValues();
      str << *TmpName << "=" << *TmpValue << endl;
    }
  return str;
}

// print last error encountered during parsing operation
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& ConfigurationParser::PrintLastError (ostream& str)
{
  if (this->ErrorLog.GetNbrElement() > 0)
    str << this->ErrorLog[this->ErrorLog.GetNbrElement() - 1];
  return str;
}

// dump all errors encountered during parsing operation
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& ConfigurationParser::DumpErrors (ostream& str)
{
  ListIterator<char*> IterError(this->ErrorLog);
  char** TmpError;
  while ((TmpError = IterError()))
    {
      str << *TmpError;
    }
  return str;
}

// clean a line from useless comments and spaces
// 
// line = pointer to the line to clean
// return value = true if the line still contains usefull information

bool ConfigurationParser::CleanLine (char* line)
{
  int NbrCharacters = strlen(line);
  if (NbrCharacters == 0)
    return false;
  int Index = 0;
  while ((Index < NbrCharacters) && ((line[Index] == ' ') || (line[Index] == '\t')))
    Index ++;
  if (Index == NbrCharacters)
    return false;
  char* Start = line + Index;
  NbrCharacters -= Index;
  Index = 0;
  while ((Index < NbrCharacters) && ((Start[Index] != '#') || ((Index > 0) && (Start[Index - 1] == '\\'))))
    Index ++;
  if (Index == 0)
    return false;
  for (int i = 0; i < Index; ++i)
    line[i] = Start[i];
  line[Index] = '\0';
  return true;
}
