#include "Vector/RealVector.h"
#include "Vector/RationalVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("VectorAscii2Binary" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-vector", "name of the file containing the ASCII vector");
  (*SystemGroup) += new SingleStringOption  ('o', "output-vector", "name of the file where the vector will be stored in binary");
  (*SystemGroup) += new BooleanOption  ('r', "rational", "indicate that the input vector has rational coefficients");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type VectorAscii2Binary -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-vector") == 0)
    {
      cout << "VectorAscii2Binary requires an input file" << endl << "see man page for option syntax or type VectorAscii2Binary -h" << endl;
      return -1;
    }

  if ((Manager.GetString("output-vector") == 0) && (Manager.GetBoolean("std-output") == false))
    {
      cout << "VectorAscii2Binary requires an output file" << endl << "see man page for option syntax or type VectorAscii2Binary -h" << endl;
      return -1;
    }

  MultiColumnASCIIFile AsciiVector;
  if (AsciiVector.Parse(Manager.GetString("input-vector")) == false)
    {
      AsciiVector.DumpErrors(cout) << endl;
      return -1;
    }

  if (Manager.GetBoolean("rational") == false)
    {
      double* TmpData = AsciiVector.GetAsDoubleArray(0);
      if (TmpData == 0)
	{
	  AsciiVector.DumpErrors(cout) << endl;
	  return -1;     
	}
      RealVector BinaryVector(TmpData, AsciiVector.GetNbrLines());
      BinaryVector.WriteVector(Manager.GetString("output-vector"));
    }
  else
    {
      long* TmpNumerators = AsciiVector.GetAsLongArray(0);
      long* TmpDenominators = AsciiVector.GetAsLongArray(1);
      if ((TmpNumerators == 0) || (TmpDenominators == 0))
	{
	  AsciiVector.DumpErrors(cout) << endl;
	  return -1;     
	}
      RationalVector BinaryVector(TmpNumerators, TmpDenominators, AsciiVector.GetNbrLines());
      BinaryVector.WriteVector(Manager.GetString("output-vector"));
    }

 
  return 0;
}
