#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/QuantumDot/Spectra/Spectra.h"
#include "Tools/QuantumDot/Spectra/DOSSpectra.h"

#include "GeneralTools/List.h"
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ArrayTools.h"


#include <iostream>
#include <fstream>
#ifdef __SSTREAM_STYLE__
#include <sstream>
#else
#include <strstream>
#endif
#include <string>
#include <unistd.h>
#include <dirent.h>


using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;

// list all files or directories that obey a given pattern (that can include relative/absolute path) /to/directory/patternxxxsuffix where xxx is an integer
//
// pattern = string that corresponds to the pattern  (i.e. /to/directory/pattern)
// matchedFileArray = reference on the sorted array (with respect to xxx) of files or directories names (with the optional relative/absolute path), 
//                    memory allocation isd one by the function itself
// suffix = optional suffix  to test
// return value = number of matched files
int GetAllDirectories(char* pattern, char**& matchedFileArray, char* suffix = 0);


int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("AbsorptionInQuantumWellBField" , "0.01");
  OptionGroup* AbsorptionInQuantumWellBFieldGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  Manager += AbsorptionInQuantumWellBFieldGroup;
  Manager += MiscGroup;

  (*AbsorptionInQuantumWellBFieldGroup) += new SingleStringOption('\n', "input", "name of the input file (directory mode only)", "eigenvalues");
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleStringOption('\n', "prefix", "prefix of the directories containing spectrum (directory mode only)", "");
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleIntegerOption('\n', "begin", "number of the first directory (directory mode only)", 0);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleIntegerOption('\n', "end", "number of the last directory (directory mode only)", 0);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleIntegerOption('n', "nbr-state", "number of states", 10);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleDoubleOption('\n', "min", "lower limit of the spectrum (in eV unit)", 0.0);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleDoubleOption('\n', "max", "upper limit of the spectrum (in eV unit)", 0.0);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleDoubleOption('g', "gamma", "full width at half maximum of each Lorentzian peak (in eV unit)", 0.01);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleDoubleOption('\n', "step", "length of each discretized step (in eV unit) in the spectrum", 2e-4);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleStringOption('\n', "output", "name of the output file", "AbsorptionInQuantumWellBField.txt");

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");
  (*MiscGroup) += new BooleanOption ('v', "verbose", "verbose mode", false);

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type AsorptionInQuantumWellBFiled -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrState = ((SingleIntegerOption*) Manager["nbr-state"])->GetInteger();
  double Min = ((SingleDoubleOption*) Manager["min"])->GetDouble();
  double Max = ((SingleDoubleOption*) Manager["max"])->GetDouble();
  double Gamma = ((SingleDoubleOption*) Manager["gamma"])->GetDouble();
  double Step = ((SingleDoubleOption*) Manager["step"])->GetDouble();
  char* OutputFile = ((SingleStringOption*) Manager["output"])->GetString();

  bool VerboseFlag = ((BooleanOption*) Manager["verbose"])->GetBoolean();


  char* InputFile = ((SingleStringOption*) Manager["input"])->GetString();
  char* Prefix = ((SingleStringOption*) Manager["prefix"])->GetString();
  int Begin = ((SingleIntegerOption*) Manager["begin"])->GetInteger();
  int End = ((SingleIntegerOption*) Manager["end"])->GetInteger();
  int Number = End - Begin + 1;
/*  char* Prefixbis = new char [100]; char* InputFilebis = new char [100];
  AddString (Prefixbis, Prefix, 0, "");  strcpy(InputFilebis, "/"); strcat(InputFilebis, InputFile);
  char** Files = new char* [Number]; int* State = new int[Number];
  for (int i = Begin; i <= End; ++i)
    {
      State[i - Begin] = NbrState;
      Files[i - Begin] = new char[200];
      if (i < 10)
	AddString(Files[i - Begin], Prefixbis, i, InputFilebis);
      else
	AddString(Files[i - Begin], Prefix, i, InputFilebis);
      if (VerboseFlag)
	cout << Files[i - Begin] << endl;
    }
  DOSSpectra AbsorptionInQuantumWellBField(Number, Files, State, Gamma, Min, Max, Step);
  AbsorptionInQuantumWellBField.WriteSpectra(OutputFile);*/

  char** MatchedDirectories;
  int NbrMatchedDirectories = GetAllDirectories("/home/regnault/results/quantumwell/quantumwell/v_2level2/bfield_20/run_0", MatchedDirectories);
  if (NbrMatchedDirectories == 0)
    {
      return -1;
    }
  char** MatchedSpectra = new char* [NbrMatchedDirectories];
  char*** MatchedEigenvectors = new char** [NbrMatchedDirectories];
  for (int i = 0; i < NbrMatchedDirectories; ++i)
    {
      int TmpLength = strlen(MatchedDirectories[i]);
      char* TmpPattern = new char[TmpLength + 2 + strlen(InputFile)];
      strcpy(TmpPattern, MatchedDirectories[i]);
      TmpPattern[TmpLength] = '/';
      strcpy(TmpPattern + TmpLength + 1, InputFile);
      cout << TmpPattern << endl;
      int TmpNbrEigenvectors = GetAllDirectories(TmpPattern, MatchedEigenvectors[i], ".vec"); 
      delete[] MatchedDirectories[i];
      delete[] TmpPattern;
    }
  
  return 1;
}


// list all files or directories that obey a given pattern (that can include relative/absolute path) /to/directory/patternxxxsuffix where xxx is an integer
//
// pattern = string that corresponds to the pattern  (i.e. /to/directory/pattern)
// matchedFileArray = reference on the sorted array (with respect to xxx) of files or directories names (with the optional relative/absolute path), 
//                    memory allocation isd one by the function itself
// suffix = optional suffix  to test
// return value = number of matched files

int GetAllDirectories(char* pattern, char**& matchedFileArray, char* suffix)
{
  char* Path = strrchr(pattern, '/');
  char* TmpPattern;
  long PatternLength;
  DIR* TmpDirectory;
  long PathLength = 0;
  if (Path == 0)
    {
      TmpDirectory = opendir(".");
      PatternLength = strlen(pattern);
      TmpPattern = new char [PatternLength + 1];
      strcpy(TmpPattern, pattern);
    }
  else
    {
      PathLength = (Path - pattern) + 1;
      PatternLength = strlen(Path) - 1;
      TmpPattern = new char [PatternLength + 1];
      strcpy(TmpPattern, Path + 1); 
      char* TmpPath = new char [PathLength + 1];
      strncpy (TmpPath, pattern, PathLength);
      TmpPath[PathLength] = '\0';
      TmpDirectory = opendir(TmpPath); 
      delete[] TmpPath;
    }
  dirent* DirectoryContent;
  List<char*> MatchedFiles;
  if (suffix == 0)
    {
      while ((DirectoryContent = readdir(TmpDirectory)))
	{
	  if ((strncmp(DirectoryContent->d_name, TmpPattern, PatternLength) == 0) && ((*(DirectoryContent->d_name + PatternLength)) >= '0') && 
	  ((*(DirectoryContent->d_name + PatternLength)) <= '9'))
	    {
	      char* TmpName  = new char [strlen(DirectoryContent->d_name) + 1];
	      strcpy (TmpName, DirectoryContent->d_name);
	      MatchedFiles += TmpName;
	    }
	}
    }
  else
    {
      while ((DirectoryContent = readdir(TmpDirectory)))
	{
	  if (strncmp(DirectoryContent->d_name, TmpPattern, PatternLength) == 0)
	    {
	      char* EndPos = DirectoryContent->d_name + PatternLength;
	      while (((*EndPos) != '\0') && ((*EndPos) >= '0') && ((*EndPos) <= '9'))
		++EndPos;
	      if (((*EndPos) != '\0') && (EndPos != (DirectoryContent->d_name + PatternLength)) && (strcmp(suffix, EndPos) == 0))
		{
		  char* TmpName  = new char [strlen(DirectoryContent->d_name) + 1];
		  strcpy (TmpName, DirectoryContent->d_name);
		  MatchedFiles += TmpName;		
		}
	    }
	}
    }
  closedir(TmpDirectory);
  if (MatchedFiles.GetNbrElement() == 0)
    {
      matchedFileArray = 0;
      return 0;
    }
  ListIterator<char*> MatchedFileIterator(MatchedFiles);
  char** TmpName2;
  int* FileIndices = new int [MatchedFiles.GetNbrElement()];
  matchedFileArray = new char* [MatchedFiles.GetNbrElement()];
  int Pos = 0;
  while ((TmpName2 = MatchedFileIterator()))
    {
      FileIndices[Pos] = atoi ((*TmpName2) + PatternLength); 
      if (PathLength == 0)
	{
	  matchedFileArray[Pos] = (*TmpName2);
	}
      else
	{
	  char* TmpName3 = new char[PathLength + 2 + strlen(*TmpName2)];
	  strncpy (TmpName3, pattern, PathLength);
	  TmpName3[PathLength] = '/';
	  strcpy (TmpName3 + PathLength + 1,(*TmpName2));	  
	  matchedFileArray[Pos] = TmpName3;
	  delete[] (*TmpName2);
	}
      ++Pos;
    }
  SortArrayUpOrdering<char*>(FileIndices, matchedFileArray, Pos);
  for (int i = 0; i < Pos; ++i)
    {
      cout << FileIndices[i] << " " << matchedFileArray[i] << endl;
    }
  delete[] FileIndices;
  delete[] TmpPattern;
  return Pos;
}
