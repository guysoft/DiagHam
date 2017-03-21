#include "Options/Options.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Matrix/ComplexMatrix.h"
#include "MPSObjects/ComplexPEPSPBC.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);
  
  // some running options and help
  OptionManager Manager ("PEPSComputeState" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  
  Manager += SystemGroup;
  
  (*SystemGroup) += new SingleStringOption  ('\n', "tensor-file", "name of the file containing the eigenstate to be displayed");
  (*SystemGroup) += new SingleStringOption  ('\n', "peps-name", "name of the peps used to form the output file name");
  (*SystemGroup) += new  SingleIntegerOption ('x',"Lx", "size of the lattice in x direction", 4);
  (*SystemGroup) += new  SingleIntegerOption ('y',"Ly", "size of the lattice in y direction (only 2 or 4)", 2);
  (*SystemGroup) += new  SingleIntegerOption ('sz',"sz", "sz value", -10000);
  (*SystemGroup) += new BooleanOption ('\n', "block-tensor", "compute the 2*2 blocked tensor (instead of computing the state)");

  (*SystemGroup) += new BooleanOption ('\n', "cylinder", "compute the state on a cylinder");
  (*SystemGroup) += new SingleStringOption  ('\n', "left-vector", "name of the file containing the left boundary vector (useful only in cylinder mode)");
  (*SystemGroup) += new SingleStringOption  ('\n', "right-vector", "name of the file containing the right boundary vector (useful only in cylinder mode)");


//  (*SystemGroup) += new BooleanOption ('c', "complex", "use complex version of the code");

  (*SystemGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PEPSComputeState -h" << endl;
      return -1;
    }

  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  int Lx = Manager.GetInteger("Lx");
  int Ly = Manager.GetInteger("Ly");
  int Sz = Manager.GetInteger("sz");
  bool BlockFlag = Manager.GetBoolean("block-tensor");
  bool CylinderFlag = Manager.GetBoolean("cylinder");
  
  MultiColumnASCIIFile TensorsElementsDefinition;
  if (TensorsElementsDefinition.Parse(Manager.GetString("tensor-file")) == false)
    {
      TensorsElementsDefinition.DumpErrors(cout) << endl;
      return -1;
    } 
  
  ComplexPEPSPBC PEPS(TensorsElementsDefinition);
  
  if (BlockFlag) 
    {
      PEPS.ComputeBlockTensor();
      return 0;
    }


  char * FullOutputFileName = new char [200];
  sprintf(FullOutputFileName,"PEPS_%s_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Lx,Ly); 

  if (CylinderFlag ==false)
    {
      if (Sz ==  -10000) 
	{
	  char * FullOutputFileName = new char [200];
	  sprintf(FullOutputFileName,"PEPS_%s_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Lx,Ly); 
	  
	  if (Ly == 2 )
	    {
	      ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPS (Lx);
	      State.WriteVector(FullOutputFileName);
	    }
	  if (Ly == 4)
	    {
	      ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPS (Lx,2);
	      State.WriteVector(FullOutputFileName);
	    }
	}
      else
	{
	  char * FullOutputFileName = new char [200];
	  sprintf(FullOutputFileName,"PEPS_%s_sz_%d_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Sz,Lx,Ly); 
	  ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPSSzConstraint (Lx,Ly/2,Sz);
	  State.WriteVector(FullOutputFileName);
	}
    }
  else
    {
      ComplexVector LeftVector;
      if (LeftVector.ReadVector(Manager.GetString("left-vector")) == false)
	{
	  cout << "error while reading " << Manager.GetString("left-vector") << endl;
	  return -1;
	} 
      
      ComplexVector RightVector;
      if (RightVector.ReadVector(Manager.GetString("right-vector")) == false)
	{
	  cout << "error while reading " << Manager.GetString("right-vector") << endl;
	  return -1;
	} 
      
      char * FullOutputFileName = new char [200];
      if (Sz ==  -10000) 
	{
	  sprintf(FullOutputFileName,"PEPS_cylinder_%s_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Lx,Ly); 
	  
	  ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPS (Lx,Ly/2, ComplexDiagonalMatrix() , LeftVector, RightVector, false,false);
	  State.WriteVector(FullOutputFileName);
	}
      else
	{
	  sprintf(FullOutputFileName,"PEPS_cylinder_%s_sz_%d_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Sz,Lx,Ly); 

	  ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPSSzConstraint (Lx, Ly/2,Sz, ComplexDiagonalMatrix() , LeftVector, RightVector, false,false);
	  State.WriteVector(FullOutputFileName);
	}
    }
}
