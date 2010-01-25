#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "Options/Options.h"

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
  OptionManager Manager ("DiffBinaryVectors" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "names of the two vector files");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-range", "compare vectors starting from a given component (is negative, start couting from the last component) from ", 0l);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-range", "compare vectors up to a given component (0 if up to the end)", 0l);
  (*SystemGroup) += new SingleDoubleOption  ('e', "error", "rounding error", 0.0);
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type DiffBinaryVectors -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrVectors;
  char** VectorFiles = Manager.GetStrings("states",NbrVectors);

  if (NbrVectors !=2)
    {
      cout << "two vector files are required!"<<endl;
      exit(1);
    }

  RealVector State1;
  if (State1.ReadVector (VectorFiles[0]) == false)
    {
      cout << "can't open vector file " << VectorFiles[0] << endl;
      return -1;      
    }
  RealVector State2;
  if (State2.ReadVector (VectorFiles[1]) == false)
    {
      cout << "can't open vector file " << VectorFiles[1] << endl;
      return -1;      
    }
  if (State1.GetLargeVectorDimension() != State2.GetLargeVectorDimension() )
    {
      cout << "Dimension of Hilbert spaces in input files does not coincide" << endl;
      return -2;
    }
  
  long MinValue = 0l;
      if (Manager.GetInteger("min-range") > 0l)
	{
	  if (Manager.GetInteger("min-range") < State1.GetLargeVectorDimension())
	    MinValue = Manager.GetInteger("min-range");      
	}
      else
	{
	  long Tmp = State1.GetLargeVectorDimension() + Manager.GetInteger("min-range");
	  if ((Tmp >= 0) && (Tmp < State1.GetLargeVectorDimension()))
	    MinValue = Tmp;
	}
  long MaxValue = State1.GetLargeVectorDimension(); 
  if ((Manager.GetInteger("max-range") < State1.GetLargeVectorDimension()) && (Manager.GetInteger("max-range") > MinValue))
    MaxValue = Manager.GetInteger("max-range");
  double Error = Manager.GetDouble("error");

  long Count = 0l;
  if (Error == 0.0)
    {
      for (long i = MinValue; i < MaxValue; ++i)
	{
	  if (State1[i] != State2[i])
	    {
	      cout << i << " : " << State1[i] << " " << State2[i] << endl;
	      ++Count;
	    }
	}
    }
  else
    {
      for (long i = MinValue; i < MaxValue; ++i)
	{
	  if ((fabs(State1[i] - State2[i]) > Error) && (fabs(State1[i] - State2[i]) > (Error * fabs(State1[i]))))
	    {
	      cout << i << " : " << State1[i] << " " << State2[i] << endl;
	      ++Count;
	    }
	}
      
    }
  cout << "total number of different components : " << Count << " / " << State1.GetLargeVectorDimension() << endl;
  return 0;
}
