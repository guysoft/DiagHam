////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class that provide multicolumn ascii file                //
//                                                                            //
//                        last modification : 31/12/2006                      //
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
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/ListIterator.h"

#include <string.h>
#include <fstream>
#include <iostream>


using std::ifstream;
using std::ios;
using std::endl;
using std::cout;


// default constructor
//
// separator = character which is used as separator between integer values in the string 
//             if \s (i.e \t or space) is used, then any number of consecutive \s or \t are identify as one separator

MultiColumnASCIIFile::MultiColumnASCIIFile(char separator)
{
  this->Separator = separator;
  this->Data = 0;
}

// destructor
//

MultiColumnASCIIFile::~MultiColumnASCIIFile()
{
  if (this->Data != 0)
    {
      for (int i = 0; i < this->NbrColumns; ++i)
	{
	  for (int j = 0; j < this->NbrLines; ++j)
	    delete[] this->Data[i][j];
	  delete[] this->Data[i];
	}
      delete[] this->Data;
    }
}

// parse a file
//
// filename = name of the file to parse
// return value = true if no error occurs

bool MultiColumnASCIIFile::Parse(char* filename)
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
  int MaxNbrLines = 0;
  while (Pos < Size)
    {
      if (TmpBuffer[Pos] == '\n')
	++MaxNbrLines;
      ++Pos;
    }
  Pos = 0;
  char* Start = TmpBuffer;
  int LineNumber = 1;
  bool Flag = true;
  this->NbrColumns = 0;
  this->NbrLines = 0;
  char** TmpArray = 0;
  while (Pos < Size)
    {
      Start = TmpBuffer + Pos;
      while ((Pos < Size) && (TmpBuffer[Pos] != '\n'))
	++Pos;
      TmpBuffer[Pos] = '\0';
      if (this->CleanLine (Start) == true)
	{
	  if (this->NbrColumns == 0)
	    {
	      this->NbrColumns = this->SplitLine(Start, TmpArray, this->Separator);
	      if (this->NbrColumns <= 0)
		{
		  char* TmpString = new char [strlen (filename) + 256]; 
		  sprintf (TmpString, "fatal error at line %d in file %s : can't retrieve number of columns\n", LineNumber, filename);
		  this->ErrorLog += TmpString;
		  return false;
		}
	      this->Data = new char** [this->NbrColumns];
	      for (int i = 0; i < this->NbrColumns; ++i)
		{
		  this->Data[i] = new char* [MaxNbrLines];
		  this->Data[i][0] = TmpArray[i];
		}	      
	      ++this->NbrLines;
	    }
	  else
	    {
	      if (this->FixedSplitLine(Start, TmpArray, this->NbrColumns, this->Separator) != this->NbrColumns)
		{
		  char* TmpString = new char [strlen (filename) + 256]; 
		  sprintf (TmpString, "fatal error at line %d in file %s : the number of columns is different from %d\n", LineNumber, filename, this->NbrColumns);
		  this->ErrorLog += TmpString;
		  return false;		  
		}
	      for (int i = 0; i < this->NbrColumns; ++i)
		this->Data[i][this->NbrLines] = TmpArray[i];
	      ++this->NbrLines;
	    }
	}
      ++LineNumber;
      ++Pos;
    }
  delete[] TmpArray;
  if (this->NbrLines != MaxNbrLines)
    {
      for (int i = 0; i < this->NbrColumns; ++i)
	{
	  char** TmpColumn = new char* [this->NbrLines];
	  for (int j = 0; j < this->NbrLines; ++j)
	    TmpColumn[j] = this->Data[i][j];
	  delete[] this->Data[i];
	  this->Data[i] = TmpColumn;
	}
    }
  return true;
}

// clean a line from useless comments and spaces
// 
// line = pointer to the line to clean
// return value = true if the line still contains usefull information

bool MultiColumnASCIIFile::CleanLine (char* line)
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

// split a line using a given separator
//
// string = string to split
// array = reference on the array where elements have to be stored (allocation is done by the method itself, de-allocation has to be done by hand)
// separator = character which is used as separator between columns 
//             if \s (i.e \t or space) is used, then any number of consecutive \s or \t are identify as one separator
// return value = number of elements in the line (zero if an error occured)

int MultiColumnASCIIFile::SplitLine(char* string, char**& array, char separator)
{
  char* End = string;
  int NbrElements = 1;
  if ((separator == ' ') || (separator == '\t'))
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if (((*End) == ' ') || ((*End) == '\t'))
	  {
	    ++NbrElements;
	    while ((((*End) != '\0') || ((*End) != '\n')) && (((*End) == ' ') || ((*End) == '\t')))
	      ++End;
	  }
	else
	  ++End;
      }
  else
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if ((*End) == separator)
	  ++NbrElements;
	++End;
      }
  array = new char*[NbrElements];  
  NbrElements = 0;
  long TmpLength;
  End = string;
  if ((separator == ' ') || (separator == '\t'))
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if (((*End) == ' ') || ((*End) == '\t'))
	  {
	    TmpLength = End - string;
	    array[NbrElements] = new char [TmpLength + 1];
	    strncpy (array[NbrElements], string, TmpLength);
	    array[NbrElements][TmpLength] = '\0';
	    ++NbrElements;
	    while ((((*End) != '\0') && ((*End) != '\n')) && (((*End) == ' ') || ((*End) == '\t')))
	      ++End;
	    string = End;
	  }
	else
	  ++End;
      }
  else
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if ((*End) == separator)
	  {
	    TmpLength = End - string;
	    array[NbrElements] = new char [TmpLength + 1];
	    strncpy (array[NbrElements], string, TmpLength);
	    array[NbrElements][TmpLength] = '\0';
	    ++NbrElements;
	    string = End + 1;
	  }
	++End;
      }
  TmpLength = End - string;
  array[NbrElements] = new char [TmpLength + 1];
  strncpy (array[NbrElements], string, TmpLength);
  array[NbrElements][TmpLength] = '\0';
  ++NbrElements;	
  return NbrElements;
}

// split a line using a given separator and requesting a fixed number of elements
//
// string = string to split
// array = pointer on the array to use to store elements
// nbrElements = number of elements to retrieve 
// separator = character which is used as separator between columns 
//             if \s (i.e \t or space) is used, then any number of consecutive \s or \t are identify as one separator
// return value = number of elements in the line (should be equal if no error occured)

int MultiColumnASCIIFile::FixedSplitLine(char* string, char** array, int nbrElements, char separator)
{
  char* End = string;
  int NbrElements = 0;
  long TmpLength;
  if ((separator == ' ') || (separator == '\t'))
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if (((*End) == ' ') || ((*End) == '\t'))
	  {
	    if (NbrElements < nbrElements)
	      {
		TmpLength = End - string;
		array[NbrElements] = new char [TmpLength + 1];
		strncpy (array[NbrElements], string, TmpLength);
		array[NbrElements][TmpLength] = '\0';
	      }
	    ++NbrElements;
	    while ((((*End) != '\0') && ((*End) != '\n')) && (((*End) == ' ') || ((*End) == '\t')))
	      ++End;
	    string = End;
	  }
	else
	  ++End;
      }
  else
    while (((*End) != '\0') && ((*End) != '\n'))
      {
	if ((*End) == separator)
	  {
	    if (NbrElements < nbrElements)
	      {
		TmpLength = End - string;
		array[NbrElements] = new char [TmpLength + 1];
		strncpy (array[NbrElements], string, TmpLength);
		array[NbrElements][TmpLength] = '\0';
	      }
	    ++NbrElements;
	    string = End + 1;
	  }
	++End;
      }
  if (NbrElements < nbrElements)
    {
      TmpLength = End - string;
      array[NbrElements] = new char [TmpLength + 1];
      strncpy (array[NbrElements], string, TmpLength);
      array[NbrElements][TmpLength] = '\0';
    }
  ++NbrElements;	
  return NbrElements;
}

// retrieve a given value (warning the string is not duplicated, thus should not be modified nor de-allocated)
//
// column = index of the column (zero based)
// line = index of the line (zero based)
// return value = string corresponding to the requested element (null pointer if out of range)

char* MultiColumnASCIIFile::operator ()(int column, int line)
{
  if ((column < 0) || (column >= this->NbrColumns) || (line < 0) || (line >= this->NbrLines))
    return 0;
  else
    return this->Data[column][line];
}

// get a column converting it to integer
//
// column = column index
// return value = reference on the array where the read values have to be stored (allocation is done by the method itself, de-allocation has to be done by han)

int* MultiColumnASCIIFile::GetAsIntegerArray (int column)
{
  int* TmpColumn = new int [this->NbrLines];
  char** TmpASCIIColumn = this->Data[column];
  char* TmpError;
  for (int i = 0; i < this->NbrLines; ++i)
    {
      long Tmp = strtol(TmpASCIIColumn[i], &TmpError, 0);
      if ((TmpError == TmpASCIIColumn[i]) || ((*TmpError) != '\0'))
	{
	  char* TmpString = new char [256 + strlen(TmpASCIIColumn[i])]; 
	  sprintf (TmpString, "element in column %d and line %d is an invalid integer value (%s)\n", column, i, TmpASCIIColumn[i]);
	  this->ErrorLog += TmpString;
	  delete[] TmpColumn;
	  return 0;
	}
      else 
	TmpColumn[i] = (int) Tmp;
    } 
  return TmpColumn;
}

// get a column converting it to long
//
// column = column index
// return value = reference on the array where the read values have to be stored (allocation is done by the method itself,  de-allocation has to be done by hand)

long* MultiColumnASCIIFile::GetAsLongArray (int column)
{
  long* TmpColumn = new long [this->NbrLines];
  char** TmpASCIIColumn = this->Data[column];
  char* TmpError;
  for (int i = 0; i < this->NbrLines; ++i)
    {
      TmpColumn[i] = strtol(TmpASCIIColumn[i], &TmpError, 0);
      if ((TmpError == TmpASCIIColumn[i]) || ((*TmpError) != '\0'))
	{
	  char* TmpString = new char [256 + strlen(TmpASCIIColumn[i])]; 
	  sprintf (TmpString, "element in column %d and line %d is an invalid long integer value (%s)\n", column, i, TmpASCIIColumn[i]);
	  this->ErrorLog += TmpString;
	  delete[] TmpColumn;
	  return 0;
	}
    } 
  return TmpColumn;
}

// get a column converting it to double
//
// column = column index
// return value = reference on the array where the read values have to be stored (allocation is done by the method itself, de-allocation has to be done by han)

double* MultiColumnASCIIFile::GetAsDoubleArray (int column)
{
  double* TmpColumn = new double [this->NbrLines];
  char** TmpASCIIColumn = this->Data[column];
  char* TmpError;
  for (int i = 0; i < this->NbrLines; ++i)
    {
      TmpColumn[i] = strtod(TmpASCIIColumn[i], &TmpError);
      if ((TmpError == TmpASCIIColumn[i]) || ((*TmpError) != '\0'))
	{
	  char* TmpString = new char [256 + strlen(TmpASCIIColumn[i])]; 
	  sprintf (TmpString, "element in column %d and line %d is an invalid double value (%s)\n", column, i, TmpASCIIColumn[i]);
	  this->ErrorLog += TmpString;
	  delete[] TmpColumn;
	  return 0;
	}
    } 
  return TmpColumn;
}

// print last error encountered during parsing operation
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& MultiColumnASCIIFile::PrintLastError (ostream& str)
{
  if (this->ErrorLog.GetNbrElement() > 0)
    str << this->ErrorLog[this->ErrorLog.GetNbrElement() - 1];
  return str;
}


// dump all errors encountered during parsing operation
//
// str = reference on the output stream
// return value = reference on the output stream

ostream& MultiColumnASCIIFile::DumpErrors (ostream& str)
{
  ListIterator<char*> IterError(this->ErrorLog);
  char** TmpError;
  while ((TmpError = IterError()))
    {
      str << *TmpError;
    }
  return str;
}
