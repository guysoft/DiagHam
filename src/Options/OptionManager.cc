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


#include "config.h"
#include "GeneralTools/ListIterator.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/OptionManager.h"

#include <iostream>
#include <string.h>


using std::ostream;
using std::endl;


// constructor
//
// groupName = group full name

OptionManager::OptionManager(char* programName, char* programVersion, char* programAdditionalInformations)
{
  this->ProgramName = new char [strlen(programName) + 1];
  strcpy (this->ProgramName, programName);      
  if (programVersion == 0)
    {
      this->ProgramVersion = 0;
    }
  else
    {
      this->ProgramVersion = new char [strlen(programVersion) + 1];
      strcpy (this->ProgramVersion, programVersion);      
    }
  if (programAdditionalInformations == 0)
    {
      this->ProgramAdditionalInformations = 0;
    }
  else
    {
      this->ProgramAdditionalInformations = new char [strlen(programAdditionalInformations) + 1];
      strcpy (this->ProgramAdditionalInformations, programAdditionalInformations);      
    }
}

// destructor
//

OptionManager::~OptionManager()
{
  delete[] this->ProgramName;
  if (this->ProgramVersion != 0)
    delete[] this->ProgramVersion;
  if (this->ProgramAdditionalInformations != 0)
    delete[] this->ProgramAdditionalInformations;
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  while ((TmpGroup = IterGroup()))
    {
      delete *TmpGroup;
    }
}

// add an option group to the manager
// 
// group = pointer to the option group to add 
// return value = reference on the current option manager

OptionManager& OptionManager::operator += (OptionGroup* group)
{
  this->Groups += group;
  return *this;
}

// get option from its name
//
// optionName = string containing option name
// return value = poitner to the option if it has been found, 0 either

AbstractOption* OptionManager::operator[] (char* optionName)
{
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  AbstractOption* TmpOption = 0;
  while ((TmpOption == 0) && ((TmpGroup = IterGroup())))
    {
      TmpOption = (**TmpGroup)[optionName];
    }
  return TmpOption; 
}

// Proceed running options from command line arguments
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// output = reference on output stream used to display errors
// return value = true if proceeding succeded, false if an error occurs

bool OptionManager::ProceedOptions (char** argumentValues, int nbrArgument, ostream& output)
{
  int Pos = 1;
  int Inc;
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  while (Pos < nbrArgument)
    {
      Inc = 0;
      IterGroup.DefineList(this->Groups);
      while ((Inc == 0) && ((TmpGroup = IterGroup())) )
	{
	  Inc = (*TmpGroup)->ReadOption(argumentValues, nbrArgument, Pos);
	  if (Inc == -1)
	    {
	      (*TmpGroup)->PrintError(output);
	      return false; 
	    }
	}
      if (Inc == 0)
	{
	  output << "unknown option " <<  argumentValues[Pos] << endl;
	  return false;
	}
      Pos += Inc;
    }
  return true;
}

// print the options and their values in the current group
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// return value = reference on current output stream

ostream& OptionManager::DisplayOption (ostream& output, bool shortVersion)
{
  if (shortVersion)
    {
      output << this->ProgramName << "  ";
      ListIterator<OptionGroup*> IterGroup(this->Groups);
      OptionGroup** TmpGroup;
      while ((TmpGroup = IterGroup()))
	{
	  (*TmpGroup)->DisplayOption(output, shortVersion);
	}
      return output;
    }
  else
    {
      if ((this->ProgramVersion != 0) || (this->ProgramAdditionalInformations != 0))
	{
	  output << this->ProgramName;
	  if (this->ProgramVersion != 0)
	    {
	      output << ", version " << this->ProgramVersion << endl;
	    }
	  else
	    {
	      output << endl;
	    }
	  if (this->ProgramAdditionalInformations != 0)
	    {
	      output << this->ProgramAdditionalInformations << endl;
	    }
	}  
      output << endl << "Options:" << endl << endl;
      ListIterator<OptionGroup*> IterGroup(this->Groups);
      OptionGroup** TmpGroup;
      while ((TmpGroup = IterGroup()))
	{
	  (*TmpGroup)->DisplayOption(output, shortVersion) << endl;
	}
      return output;
    }
}

// print help concerning current option group
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& OptionManager::DisplayHelp (ostream& output)
{
  if ((this->ProgramVersion != 0) || (this->ProgramAdditionalInformations != 0))
    {
      output << this->ProgramName;
      if (this->ProgramVersion != 0)
	{
	   output << ", version " << this->ProgramVersion << endl;
	}
      else
	{
	  output << endl;
	}
      if (this->ProgramAdditionalInformations != 0)
	{
	  output << this->ProgramAdditionalInformations << endl;
	}
      output << endl;      
    }  
  output << "Usage: " << this->ProgramName << " [options] " << endl;
  output << endl << "Options:" << endl << endl;
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  while ((TmpGroup = IterGroup()))
    {
      (*TmpGroup)->DisplayHelp(output) << endl;
    }
  return output;
}

