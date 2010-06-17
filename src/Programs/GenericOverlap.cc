#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("GenericOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "names of the vector files obtained using exact diagonalization");
  (*SystemGroup) += new SingleStringOption  ('\n', "state-list", "provide the names of the vector using a single column formatted text file");
  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new BooleanOption  ('s', "scalar-product", "Get the scalar product, not the overlap");
  (*SystemGroup) += new BooleanOption  ('\n', "conjugate", "Conjugate the second (complex) number");
  (*SystemGroup) += new BooleanOption  ('\n', "discard-sign", "compute sum_i |v1_i * v2_i| instead of sum_i v1_i * v2_i");
  (*SystemGroup) += new BooleanOption  ('x', "no-cross", "calculate only overlap of 1st vector with all others");
  (*SystemGroup) += new BooleanOption  ('\n', "sum", "sum all computed overlaps");
  (*SystemGroup) += new BooleanOption  ('\n', "no-square", "calculate only the scalar products");
  (*SystemGroup) += new BooleanOption  ('n', "normalize", "normalize vectors before calculating any overlaps");
  (*SystemGroup) += new BooleanOption  ('d', "dimension", "show vector dimension");
  (*SystemGroup) += new BooleanOption  ('\n', "quiet", "discard any output except the overlaps");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericOverlap -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  bool QuietFlag = Manager.GetBoolean("quiet");
  bool Scalar = Manager.GetBoolean("scalar-product");
  
  int NbrVectors;
  char** VectorFiles = 0;

  if (Manager.GetString("state-list") == 0)
    {
      VectorFiles = Manager.GetStrings("states", NbrVectors);
    }
  else
    {
      MultiColumnASCIIFile Description;
      if (Description.Parse(Manager.GetString("state-list")) == false)
	{
	  Description.DumpErrors(cout);
	  return -1;
	}      
      if (Description.GetNbrColumns() < 1)
	{
	  cout << "wrong number of columns in " << Manager.GetString("state-list") << endl;
	  return -1;
	}
      if (Manager.GetStrings("states", NbrVectors) != 0)
	{
	  int TmpNbrVectors1 = 0;
	  char** TmpVectorFiles1 =  Manager.GetStrings("states", TmpNbrVectors1);
	  char** TmpVectorFiles2 = Description.GetAsStringArray(0);
	  int TmpNbrVectors2 = Description.GetNbrLines();
	  NbrVectors = TmpNbrVectors1 + TmpNbrVectors2;
	  VectorFiles = new char* [NbrVectors];
	  for (int i = 0; i < TmpNbrVectors1; ++i)
	    VectorFiles[i] = TmpVectorFiles1[i];
	  for (int i = 0; i < TmpNbrVectors2; ++i)
	    VectorFiles[i + TmpNbrVectors1] = TmpVectorFiles2[i];
	  delete[] TmpVectorFiles2;
	}
      else
	{
	  VectorFiles = Description.GetAsStringArray(0);
	  NbrVectors = Description.GetNbrLines();
	}
    }

  if (NbrVectors < 2)
    {
      cout << "At least two vector files are required!"<<endl;
      exit(1);
    }

  if (QuietFlag == false)
    for (int i=0; i<NbrVectors; ++i)    
      cout << "File "<<i<<"  "<<VectorFiles[i]<<endl;

  Complex sp=0.0;

  int MaxVectors=(Manager.GetBoolean("no-cross")?1:NbrVectors);

  bool HaveComplex=Manager.GetBoolean("complex");
  double TotalOverlap = 0.0;
  
  if (HaveComplex)
    {      
      ComplexVector State1, State2;
      for (int i=0; i<MaxVectors; ++i)
	{
	  if (State1.ReadVector (VectorFiles[i]) == false)
	    {
	      cout << "can't open vector file " << VectorFiles[i] << endl;
	      return -1;      
	    }
	  if ((i==0)&&(Manager.GetBoolean("dimension")))
	    {
	      cout << "Vector dimension = "<<State1.GetVectorDimension() <<endl;
	    }
	  if (Manager.GetBoolean("normalize"))
	    State1/=State1.Norm();
	  for (int j=i+1; j<NbrVectors; ++j)
	    {	      
	      if (State2.ReadVector (VectorFiles[j]) == false)
		{
		  cout << "can't open vector file " << VectorFiles[j] << endl;
		  return -1;      
		}
	      if (State1.GetVectorDimension() != State2.GetVectorDimension() )
		{
		  cout << "Dimension of Hilbert spaces in input files does not coincide" << endl;
		  return -2;
		}
	      if (Manager.GetBoolean("normalize"))
		State2/=State2.Norm();
	      sp=0.0;
	      if (Manager.GetBoolean("discard-sign"))
		for (int i=0; i<State1.GetVectorDimension(); ++i)
		  sp+=Norm(State1[i]*State2[i]);
	      else
		if (Manager.GetBoolean("conjugate"))
		  for (int i=0; i<State1.GetVectorDimension(); ++i)
		    sp+= State1[i]*State2[i];
		else
		  for (int i=0; i<State1.GetVectorDimension(); ++i)
		    sp+= Conj(State1[i])*State2[i];
	      TotalOverlap += SqrNorm(sp);
	      if (Scalar==false)
		{
		  if (QuietFlag == false)
		    cout << "Overlap |<"<<i<<"|"<<j<<">|^2 = " << SqrNorm(sp) << endl;
		  else
		    cout << SqrNorm(sp) << endl;
		}
	      else
		{
		  if (QuietFlag == false)
		    cout << "Overlap |<"<<i<<"|"<<j<<">|^2 = " << sp << endl;
		  else
		    cout << sp.Re << " " << sp.Im << endl;
		}
	    }
	}
    }
  else // real vectors
    {
      RealVector State1, State2;
      for (int i=0; i<MaxVectors; ++i)
	{
	  if (State1.ReadVector (VectorFiles[i]) == false)
	    {
	      cout << "can't open vector file " << VectorFiles[i] << endl;
	      return -1;      
	    }
	  if ((i==0)&&(Manager.GetBoolean("dimension")))
	    {
	      cout << "Vector dimension = "<<State1.GetVectorDimension() <<endl;
	    }
	  if (Manager.GetBoolean("normalize"))
	    State1/=State1.Norm();
	  for (int j=i+1; j<NbrVectors; ++j)
	    {	      
	      if (State2.ReadVector (VectorFiles[j]) == false)
		{
		  cout << "can't open vector file " << VectorFiles[j] << endl;
		  return -1;      
		}
	      if (State1.GetVectorDimension() != State2.GetVectorDimension() )
		{
		  cout << "Dimension of Hilbert spaces in input files does not coincide" << endl;
		  return -2;
		}
	      if (Manager.GetBoolean("normalize"))
		State2/=State2.Norm();
	      sp=0.0;
	      if (Manager.GetBoolean("discard-sign"))
		for (int i=0; i<State1.GetVectorDimension(); ++i)
		  sp+=fabs(State1[i]*State2[i]);
	      else
		sp = State1 *State2 ;
	      TotalOverlap += SqrNorm(sp);
	      if (Scalar==false)
		{
		  if (QuietFlag == false)
		    cout << "Overlap |<"<<i<<"|"<<j<<">|^2 = " << SqrNorm(sp) << endl;
		  else
		    cout << SqrNorm(sp) << endl;
		}
	      else
		{
		  if (QuietFlag == false)
		    cout << "<"<<i<<"|"<<j<<"> = " << sp.Re << endl;
		  else
		    cout << sp.Re << endl;
		}

	    }
	}
    }
  if (Manager.GetBoolean("sum") == true)
    {
      if (QuietFlag == false)
	cout << "Total overlap = " << TotalOverlap << endl;
      else
	cout << TotalOverlap << endl;
    }
  return 0;
}
