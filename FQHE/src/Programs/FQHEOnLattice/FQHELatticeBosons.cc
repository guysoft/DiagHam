#include "HilbertSpace/BosonOnLattice.h"
#include "Hamiltonian/ParticleOnLatticeDeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnLatticeMainTask.h"

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


#include "Matrix/ComplexMatrix.h"

// store imaginary Hamiltonian into an complex matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

ComplexMatrix& GetHamiltonianIm (AbstractHamiltonian *H, ComplexMatrix& M)
{
  ComplexVector TmpV1 (H->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV2 (H->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < H->GetHilbertSpaceDimension(); i++)
    {
      TmpV1.Im(i) = 1.0;
      H->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < H->GetHilbertSpaceDimension(); j++)
	{
	  M.SetMatrixElement(i, j, TmpV2[j]);
	}
      TmpV1.Im(i) = 0.0;
    }
  return M;  
}

// store imaginary Hamiltonian into an complex matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

ComplexMatrix& GetHamiltonian (AbstractHamiltonian *H, ComplexMatrix& M)
{
  ComplexVector TmpV1 (H->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV2 (H->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < H->GetHilbertSpaceDimension(); i++)
    {
      TmpV1.Re(i) = 1.0;
      H->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < H->GetHilbertSpaceDimension(); j++)
	{
	  M.SetMatrixElement(i, j, TmpV2[j]);
	}
      TmpV1.Re(i) = 0.0;
    }
  return M;  
}


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeBosons" , "0.01");  
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  QHEOnLatticeMainTask::AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 8);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 5);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 1);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice (-1=all)", -1);
  (*SystemGroup) += new SingleDoubleOption  ('u', "contactU", "prefactor U of the contact interaction (kinetic term ~ 1)", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('d', "deltaPotential", "Introduce a delta-potential at the origin", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "negative-hopping", "reverse sign of hopping terms", false);
  
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  int NbrSites = Lx*Ly;  
  bool ReverseHopping = Manager.GetBoolean("negative-hopping");
  double ContactU = Manager.GetDouble("contactU");
  double Delta = Manager.GetDouble("deltaPotential");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  bool FirstRun = true;

  int NbrFluxValues = 1;
  if (NbrFluxQuanta == -1)
    {
      NbrFluxQuanta = 0;
      NbrFluxValues = NbrSites;
    }

  char* OutputName;
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      OutputName = new char [256];
      char reverseHoppingString[4]="";
      char deltaString[20]="";
      if (ReverseHopping)
	sprintf(reverseHoppingString,"rh_");
      if (Delta!=0.0)
	sprintf(deltaString,"d_%g_",Delta);  
      if (NbrFluxValues == 1)
	sprintf (OutputName, "bosons_lattice_n_%d_x_%d_y_%d_u_%g_%s%sq_%d.dat", NbrBosons, Lx, Ly, ContactU, reverseHoppingString, deltaString, NbrFluxQuanta);
      else
	sprintf (OutputName, "bosons_lattice_n_%d_x_%d_y_%d_u_%g_%s%sq.dat", NbrBosons, Lx, Ly, ContactU, reverseHoppingString, deltaString);
    }
  ParticleOnLattice* Space=new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  
  AbstractQHEOnLatticeHamiltonian* Hamiltonian;
  Hamiltonian = new ParticleOnLatticeDeltaHamiltonian(Space, NbrBosons, Lx, Ly, NbrFluxQuanta, ContactU,
						      ReverseHopping, Delta, Architecture.GetArchitecture(), Memory, LoadPrecalculationFileName);

  // testing Hamiltonian:
  
  ComplexMatrix HRe(Hamiltonian->GetHilbertSpaceDimension(),Hamiltonian->GetHilbertSpaceDimension());
  ComplexMatrix HIm(Hamiltonian->GetHilbertSpaceDimension(),Hamiltonian->GetHilbertSpaceDimension());
  Complex one, two, M_I(0.0,1.0);
  for (int i=0; i<Hamiltonian->GetHilbertSpaceDimension(); ++i)
    for (int j=0; j<Hamiltonian->GetHilbertSpaceDimension(); ++j)
      {
	HRe.GetMatrixElement(i,j,one);
	HIm.GetMatrixElement(i,j,two);
	one= one*M_I;
	if (Norm(one-two)>1e-10)
	  cout << "Discrepancy in "<<i<<", "<<j<<": "<<one << " vs " << two << endl;
      }

  ComplexVector TmpV1a (Hamiltonian->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV1b (Hamiltonian->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV2a (Hamiltonian->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV2b (Hamiltonian->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < Hamiltonian->GetHilbertSpaceDimension(); i++)
    {
      TmpV1a.Re(i) = (rand() - 32767) * 0.5;
      TmpV1a.Im(i) = (rand() - 32767) * 0.5;
    }
  TmpV1a /= TmpV1a.Norm();
  TmpV1b = TmpV1a*M_I;
  Hamiltonian->LowLevelMultiply(TmpV1a, TmpV2a);
  Hamiltonian->LowLevelMultiply(TmpV1b, TmpV2b);
  for (int j=0; j<Hamiltonian->GetHilbertSpaceDimension(); ++j)
      {
	one = TmpV2a[j];
	two = TmpV2b[j];
	one = one*M_I;
	if (Norm(one-two)>1e-10)
	  cout << "Discrepancy in "<<j<<": "<<one << " vs " << two << endl;
      }  


  for (int iter=0; iter<NbrFluxValues; ++iter, ++NbrFluxQuanta)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << "NbrFluxQuanta="<<NbrFluxQuanta<<endl;
      
      if (!FirstRun) Hamiltonian->SetNbrFluxQuanta(NbrFluxQuanta);
  
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate"))	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "bosons_lattice_n_%d_x_%d_y_%d_u_%g_q_%d", NbrBosons, Lx, Ly, ContactU, NbrFluxQuanta);
	}
      QHEOnLatticeMainTask Task (&Manager, Space, Hamiltonian, NbrFluxQuanta, 0.0, OutputName, FirstRun, EigenvectorName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      if (FirstRun == true)
	FirstRun = false;
    }
  
  delete Hamiltonian;
  delete Space;  

  return 0;
}
