#include "config.h"

#include "Vector/RealVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHENBodyQuasiHoleOverlap" , "0.01");
  OptionGroup* MainGroup = new OptionGroup ("main options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  (*MainGroup) += new SingleStringOption  ('\0', "input-file", "name of the file which contains definition of overlaps to evaluate");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHENBodyQuasiHoleOverlap -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  ConfigurationParser OverlapDefinition;
  if (OverlapDefinition.Parse(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      OverlapDefinition.DumpErrors(cout) << endl;
      return -1;
    }


  if (OverlapDefinition["Degeneracy"] == 0)
    {
      cout << "no Degeneracy defined in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  char* TmpString = OverlapDefinition["Degeneracy"];
  int MaxNbrLz = (strlen(TmpString) >> 1) + 1;
  if (MaxNbrLz == 0)
     {
      cout << "error while parsing Degeneracy in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  int* Degeneracy = new int [MaxNbrLz];
  MaxNbrLz = 0;
  char* End = TmpString;
  bool ErrorFlag = false;
  while (((*TmpString) != '\0') && (ErrorFlag == false))
    {
      while (((*End) != '\0') && ((*End) != ','))
	{
	  if ((((*End) < '0') || ((*End) > '9')) && ((*End) != ' ') && ((*End) != '\t'))
	    ErrorFlag = true;
	  ++End;
	}
      if (ErrorFlag == true)
	{
	  cout << "error while parsing Degeneracy in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
	  return -1;     
	}
      Degeneracy[MaxNbrLz] = atoi (TmpString);
      ++MaxNbrLz;
    }

  if (OverlapDefinition["InputVectors"] == 0)
    {
      cout << "no InputVectors defined in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  if (OverlapDefinition["OutputVectors"] == 0)
    {
      cout << "no OutputVectors defined in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  
  char* InputVectors = new char [strlen(OverlapDefinition["InputVectors"]) + 24];
  char* OutputVectors = new char [strlen(OverlapDefinition["OutputVectors"]) + 24];
  strcpy (InputVectors, OverlapDefinition["InputVectors"]);
  strcpy (OutputVectors, OverlapDefinition["OutputVectors"]);
  char* InputVectors2 = InputVectors + strlen(InputVectors);
  char* OutputVectors2 = OutputVectors + strlen(OutputVectors);

  for (int i = 0; i < MaxNbrLz; ++i)
    {
      int Lz = i << 1;
      cout << "Lz = " << i << endl;
      for (int j = 0; j < Degeneracy[i]; ++j)
	{
	  sprintf (InputVectors2, "%d.%d.vec", Lz, j);
	  sprintf (OutputVectors2, "%d.%d.vec", Lz, j);
	  cout << InputVectors << " " << OutputVectors << endl;
	}
      cout << "----------------------------------------------------" << endl;
    }

  delete[] Degeneracy;
  delete[] InputVectors;
  delete[] OutputVectors;

  return 0;
}
