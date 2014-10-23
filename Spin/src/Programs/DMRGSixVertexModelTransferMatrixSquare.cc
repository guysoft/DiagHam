#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"


#include "MPSObjects/MPOPeratorSixVertexModelTransferMatrixSquare.h"
#include "MPSObjects/MPSSite.h"
#include "MPSObjects/DMRGFiniteSizeRealOBCMainTask.h" 
#include "Options/Options.h"
#include "Matrix/RealMatrix.h"


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
  Architecture.AddOptionGroup(&Manager);
  Manager += SystemGroup;
  Manager += MiscGroup;

  
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
  
  
  int NbrSites = 4;
  int PhysicalDimension = 2;
  int MaxBondDimension = 10;
  MPOPeratorSixVertexModelTransferMatrixSquare TransferMatrix(NbrSites);  
  TransferMatrix.PrintTensorElements();

  MPSSite * Lattice = new MPSSite[NbrSites];
  Lattice[0] = MPSSite(0, PhysicalDimension, 0, Lattice + 1, MaxBondDimension,&TransferMatrix);
  for(int i = 1 ; i < NbrSites - 1 ; i++ )
    {
      Lattice[i] = MPSSite(i, PhysicalDimension, &Lattice[i-1], &Lattice[i+1], MaxBondDimension, &TransferMatrix);
    }
  Lattice[NbrSites-1] = MPSSite(NbrSites-1, PhysicalDimension, &Lattice[NbrSites-2], 0, 10,&TransferMatrix);
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
        NextCurrentDimension=   MaxBondDimension;
     }	
  }
  int NbrSweep = 1;
  DMRGFiniteSizeRealOBCMainTask Algorithm(Lattice, &TransferMatrix, NbrSites, NbrSweep, MaxBondDimension, Architecture.GetArchitecture());
  Algorithm.RunAlgorithm();
  delete [] Lattice;
  return 0;
}
