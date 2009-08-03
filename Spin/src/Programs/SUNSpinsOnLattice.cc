#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"

#include "HilbertSpace/GenericSUNSpinCollection.h"
#include "Hamiltonian/SUNSpinOnLatticeQuadraticHamiltonian.h"
#include "Tools/LatticeConnections.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <iostream>

using std::cout;
using std::endl;

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("SpinOnLatticeWithSUN" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");      
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);  
  Manager += SystemGroup;
  LatticeConnections::AddOptionGroup(&Manager);  
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('n', "level-n", "level of SU(N) symmetry group", 3);
  (*SystemGroup) += new MultipleIntegerOption  ('t', "cartan", "eigenvalues of the generators of the cartan algebra",',');

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout); 

  int LevelN=Manager.GetInteger("level-n");
  int NbrCartan;
  int *CartanQuantumNumbers;
  CartanQuantumNumbers = Manager.GetIntegers("cartan",NbrCartan);
  unsigned long MemorySpace = ((unsigned long)Manager.GetInteger("fast-search")) << 20;
  if (ULONG_MAX>>20 < (unsigned long)Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  unsigned long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  char *LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  
  if (NbrCartan<LevelN-1)
    {
      cout << "Please select the desired subspace!"<<endl
	   << "Use option -t to indicate the first N-1 eigenvalues of the Cartan operators"<<endl;
      exit(-1);
    }
  if (NbrCartan>LevelN-1)
    cout << "More values of Cartan operators than needed. Will ignore all beyond N-1"<<endl;  
 
  LatticeConnections *Lattice = new LatticeConnections();

  int NbrSpins=Lattice->GetNbrSites();

  GenericSUNSpinCollection *Space = new GenericSUNSpinCollection(LevelN, NbrSpins, CartanQuantumNumbers, MemorySpace);
  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  
  AbstractSUNSpinHamiltonian *Hamiltonian =
    new SUNSpinOnLatticeQuadraticHamiltonian(Space, Lattice, Architecture.GetArchitecture(),
					     Memory, LoadPrecalculationFileName);

  delete Space;
  delete Lattice;
  delete Hamiltonian;
  delete [] CartanQuantumNumbers;
}
