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
  OptionManager Manager ("GenericOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "names of all vector files that should be superposed");
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new SingleDoubleOption  ('r', "random-component", "amplitude of a random component to be added",0.5);
  (*SystemGroup) += new SingleStringOption  ('o', "output", "names of output filename","superposition.vec");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrVectors, ValidVectors=0, VectorDimension=0;
  char** VectorFiles = Manager.GetStrings("states",NbrVectors);


  Vector **Vectors = new Vector*[NbrVectors], *Result, *Random;

  bool tmpB, haveVector=false;
  
  if (Manager.GetBoolean("complex"))
    {
      Result = new ComplexVector();
      for (int i=0; i<NbrVectors; ++i)
	{
	  Vectors[i]= new ComplexVector();
	  tmpB = ((ComplexVector*)Vectors[i])->ReadVector(VectorFiles[i]);
	  if (!haveVector)
	    VectorDimension=((ComplexVector*)Vectors[i])->GetVectorDimension();
	  if (haveVector && (((ComplexVector*)Vectors[i])->GetVectorDimension()!=VectorDimension))
	    {
	      cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of previous vectors!"<<endl;
	      exit(1);
	    }
	  haveVector=haveVector | tmpB;
	  if (tmpB) ++ValidVectors;
	}
      if (ValidVectors==0)
	{
	  cout <<"At least one valid vector is required!"<<endl;
	  exit(1);
	}
      Random = new ComplexVector(VectorDimension);
      for (int i = 0; i < VectorDimension; ++i)
	{
	  ((ComplexVector*)Random)->Re(i) = (rand() - 32767) * 0.5;
	  ((ComplexVector*)Random)->Im(i) = (rand() - 32767) * 0.5;
	}
      *((ComplexVector*)(Random)) /= Random->Norm();
    }
  else
    {
      Result = new RealVector();
      for (int i=0; i<NbrVectors; ++i)
	{
	  Vectors[i] = new RealVector();
	  tmpB = ((RealVector*)Vectors[i])->ReadVector(VectorFiles[i]);
	  if (!haveVector)
	    VectorDimension=((RealVector*)Vectors[i])->GetVectorDimension();
	  if (haveVector && (((RealVector*)Vectors[i])->GetVectorDimension()!=VectorDimension))
	    {
	      cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of previous vectors!"<<endl;
	      exit(1);
	    }
	  haveVector=haveVector | tmpB;
	  if (tmpB) ++ValidVectors;
	}
      if (ValidVectors==0)
	{
	  cout <<"At least one valid vector is required!"<<endl;
	  exit(1);
	}
      Random = new RealVector(VectorDimension);
      for (int i = 0; i < VectorDimension; ++i)
	(*(RealVector*)(Random))[i] = (rand() - 32767) * 0.5;
      *((RealVector*)(Random)) /= Random->Norm();
    }	

  double RandomComponent = Manager.GetDouble("random-component");  
  double GivenComponent = (1.0-RandomComponent)/NbrVectors;

  Result->Resize(VectorDimension);
  Result->ClearVector();

  for (int i=0; i<NbrVectors; ++i)
    {
      cout<<"Adding Vector "<<i<<endl;
      Result->AddLinearCombination(GivenComponent,*(Vectors[i]));
    }
  
  Result->AddLinearCombination(RandomComponent,*Random);

  char *OutputName = Manager.GetString("output");

  if (Manager.GetBoolean("complex"))
    {
      *((ComplexVector*)(Result)) /= Result->Norm();
      ((ComplexVector*)Result)->WriteVector(OutputName);
    }
  else
    {
      *((RealVector*)(Result)) /= Result->Norm();
      ((RealVector*)Result)->WriteVector(OutputName);
    }
  
}
  
