#include "HilbertSpace/FermionOnS2xS2.h"

#include "Hamiltonian/ParticleOnS2xS2DeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/StringTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteY = number of sites along the y direction
// return value = Hilbert space dimension
long FQHES2xS2FermionsEvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, 
						    int kxMomentum, int kyMomentum, int nbrSiteY);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHES2xS2FermionsTwoBodyGeneric" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux1", "number of flux quanta for the first sphere", 0);
  (*SystemGroup) += new SingleIntegerOption  ('k', "nbr-flux2", "number of flux quanta for the second sphere", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHES2xS2FermionsTwoBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrFermions = Manager.GetInteger("nbr-particles");
  int NbrFluxQuanta1 = Manager.GetInteger("nbr-flux1");
  int NbrFluxQuanta2 = Manager.GetInteger("nbr-flux2");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  
  bool FirstRun = true;

  char* OutputName = new char [256];
  sprintf (OutputName, "fermions_s2xs2_delta_n_%d_2s1_%d_2s2_%d.dat", NbrFermions, NbrFluxQuanta1, NbrFluxQuanta2);

//   int MaxTotalLz = (NbrFluxQuanta1 * NbrFermions) - (NbrFermions * (NbrFermions - 1));
//   int MaxTotalKz = (NbrFluxQuanta2 * NbrFermions) - (NbrFermions * (NbrFermions - 1));
//   int MinTotalLz = -MaxTotalLz;
//   int MinTotalKz = -MaxTotalKz;

  int MaxTotalLz = ((NbrFluxQuanta1 * (NbrFermions / (NbrFluxQuanta2 + 1))) - ((NbrFermions / (NbrFluxQuanta2 + 1)) * ((NbrFermions / (NbrFluxQuanta2 + 1)) - 1))) * (NbrFluxQuanta2 + 1);
  int MaxTotalKz = 0;//(NbrFluxQuanta2 * NbrFermions) - (NbrFermions * (NbrFermions - 1));
  int MinTotalLz = -MaxTotalLz;
  int MinTotalKz = -MaxTotalKz;
  
  for (int TotalLz = MinTotalLz; TotalLz <= MaxTotalLz; TotalLz += 2)
    {
      for (int TotalKz = MinTotalKz; TotalKz <= MaxTotalKz; TotalKz += 2)
	{
	  int KxMomentum = (TotalLz + (NbrFermions * NbrFluxQuanta1)) >> 1;
	  int KyMomentum = (TotalKz + (NbrFermions * NbrFluxQuanta2)) >> 1;
 	  cout << ((TotalLz - MinTotalLz) >> 1) << " " << ((TotalKz - MinTotalKz) >> 1) << " " << FQHES2xS2FermionsEvaluateHilbertSpaceDimension(NbrFermions, NbrFluxQuanta1, NbrFluxQuanta2, 0, 0, 
																		 KxMomentum, KyMomentum, NbrFluxQuanta2 + 1) << endl;
	}
    }



  for (int TotalLz = MinTotalLz; TotalLz <= MaxTotalLz; TotalLz += 2)
    {
      for (int TotalKz = MinTotalKz; TotalKz <= MaxTotalKz; TotalKz += 2)
	{
// 	  ParticleOnSphere* Space = 0;
// 	  Space = new FermionOnS2xS2(NbrFermions, NbrFluxQuanta1, NbrFluxQuanta2, TotalLz, TotalKz);
// 	  cout << (MaxTotalLz - TotalLz) << " " << (MaxTotalKz - TotalKz) << " " << Space->GetHilbertSpaceDimension() << endl;
	  
// 	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
// 	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
// 	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
// 	  //       for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
// 	  // 	Space->PrintState(cout, i);
	  
// 	  AbstractQHEHamiltonian* Hamiltonian = 0;
// 	  Hamiltonian = new ParticleOnS2xS2DeltaHamiltonian(Space, NbrFermions, NbrFluxQuanta1, NbrFluxQuanta2, Architecture.GetArchitecture(), Memory);
      
// 	  char* EigenvectorName = 0;
// 	  if (Manager.GetBoolean("eigenstate") == true)	
// 	    {
// 	      char* TmpVectorExtension = new char [64];
// 	      sprintf (TmpVectorExtension, "_lz_%d_kz_%d", TotalLz, TotalKz);
// 	      EigenvectorName = ReplaceString(OutputName, ".dat", TmpVectorExtension);
// 	    }
	
// 	  char* ContentPrefix = new char[256];
// 	  sprintf (ContentPrefix, "%d %d", TotalLz, TotalKz);	  
// 	  char* SubspaceLegend = new char[256];
// 	  sprintf (SubspaceLegend, "Lz Kz");

// 	  GenericRealMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, ContentPrefix, SubspaceLegend, 0, OutputName, FirstRun, EigenvectorName);
// 	  MainTaskOperation TaskOperation (&Task);	  	 
// 	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
// 	  delete Hamiltonian;
// 	  delete Space;
// 	  if (EigenvectorName != 0)
// 	    {
// 	      delete[] EigenvectorName;
// 	    }
// 	  if (FirstRun == true)
// 	    FirstRun = false;
	}       
    }
  return 0;
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteY = number of sites along the y direction
// return value = Hilbert space dimension

long FQHES2xS2FermionsEvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, 
						    int kxMomentum, int kyMomentum, int nbrSiteY)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if ((currentTotalKx == kxMomentum) && (currentTotalKy == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  Count += FQHES2xS2FermionsEvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, kxMomentum, kyMomentum, nbrSiteY);
  Count += FQHES2xS2FermionsEvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, kxMomentum, kyMomentum, nbrSiteY);
  return Count;
}

