#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"


#include "MPSObjects/MPOPeratorDefinedByFiles.h"
#include "MPSObjects/MPSSite.h"
#include "MPSObjects/DMRGFiniteSizeRealOBCMainTask.h" 


#include "Options/Options.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Matrix/RealMatrix.h"

#include "LanczosAlgorithm/LanczosManager.h"

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
  OptionManager Manager ("DMRGSixVertexModelTransferMatrixSquare" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Architecture.AddOptionGroup(&Manager);
  Manager += SystemGroup;
  Manager += MiscGroup;
  Lanczos.AddOptionGroup(&Manager);

  (*SystemGroup) += new  SingleIntegerOption ('L', "length", "length of the spin chain", 4);
  (*SystemGroup) += new  SingleIntegerOption ('D', "bond-dimension", "bond dimension", 20);
  (*SystemGroup) += new  SingleIntegerOption ('s', "sweep", "number of sweep to be performed", 4);
  (*SystemGroup) += new SingleStringOption  ('\n', "tensor-file", "name of the file containing the eigenstate to be displayed");
  (*SystemGroup) += new SingleStringOption  ('\n', "vector-file", "name of the file containing the eigenstate to be displayed");

  (*MiscGroup) += new  BooleanOption ('\n', "print-tensor", "print the tensor elements", false);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericOverlap -h" << endl;
      return -1;
    }

  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  

  MultiColumnASCIIFile TensorsElementsDefinition;
  if (TensorsElementsDefinition.Parse(Manager.GetString("tensor-file")) == false)
	{
	  TensorsElementsDefinition.DumpErrors(cout) << endl;
	  return -1;
	} 

  if (TensorsElementsDefinition.GetNbrColumns() != 5)
{
  cout <<" The tensor file should have 5 columnns"<<endl;
}

  MultiColumnASCIIFile BoundaryVectorsDefinition;
  if (BoundaryVectorsDefinition.Parse(Manager.GetString("vector-file")) == false)
	{
	  BoundaryVectorsDefinition.DumpErrors(cout) << endl;
	  return -1;
	} 

  if ((BoundaryVectorsDefinition.GetNbrColumns() > 2) || ( BoundaryVectorsDefinition.GetNbrColumns() < 1 ) )
  {
      cout <<" The vector should be defined and cannot be more than 2"<<endl;
  }


  if(Manager.GetBoolean("test-idmrg") == true)
  {
  
   return 0;
  }

 
  int NbrSites = Manager.GetInteger("length");
  MPOPeratorDefinedByFiles TransferMatrix(NbrSites,TensorsElementsDefinition,BoundaryVectorsDefinition);
  int PhysicalDimension = TransferMatrix.GetPhysicalDimension();
  int MaxBondDimension = Manager.GetInteger("bond-dimension");;

 
if(Manager.GetBoolean("print-tensor") == true)
  {
  TransferMatrix.PrintTensorElements();
   return 0;
  }

  MPSSite * Lattice = new MPSSite[NbrSites];
  Lattice[0] = MPSSite(0, PhysicalDimension, 0, &Lattice[1], MaxBondDimension,&TransferMatrix);
  for(int i = 1 ; i < NbrSites - 1 ; i++ )
    {
      Lattice[i] = MPSSite(i, PhysicalDimension, &Lattice[i-1], &Lattice[i+1], MaxBondDimension, &TransferMatrix);
    }
  Lattice[NbrSites-1] = MPSSite(NbrSites-1, PhysicalDimension, &Lattice[NbrSites-2], 0, MaxBondDimension,&TransferMatrix);
  int CurrentDimension = 1;
  int NextCurrentDimension = PhysicalDimension;
  for(int i = 0;  i < (NbrSites>>1) ; i++)
  {
    Lattice[i].SetBondDimension(CurrentDimension,NextCurrentDimension);
    Lattice[NbrSites - i - 1].SetBondDimension(NextCurrentDimension,CurrentDimension);

   CurrentDimension = NextCurrentDimension;
   NextCurrentDimension *=  PhysicalDimension;
   if(NextCurrentDimension > MaxBondDimension)
     { 
        NextCurrentDimension =   MaxBondDimension;
     }	
  }
  int NbrSweep = Manager.GetInteger("sweep");
  DMRGFiniteSizeRealOBCMainTask Algorithm (Lattice, &TransferMatrix, NbrSites, NbrSweep, MaxBondDimension, Architecture.GetArchitecture(), &Lanczos);
  Algorithm.RunAlgorithm();
  delete [] Lattice;
  return 0;
}
