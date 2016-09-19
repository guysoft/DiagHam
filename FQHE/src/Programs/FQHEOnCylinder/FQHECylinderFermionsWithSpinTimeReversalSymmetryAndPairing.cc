#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinAndPairing.h"
#include "HilbertSpace/FermionOnSphereWithSpinAndPairingAllLz.h"

#include "Hamiltonian/ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairing.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHECylinderFermionsWithSpinTimeReversalSymmetryAndPairing" , "0.01");
  OptionGroup* SystemGroup  = new OptionGroup("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);
  

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 8);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-ky", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-ky", "number of ky value to evaluate", -1);
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "if positive, fix the total number of particles", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "force-negativeky", "manually force to compute the negative ky sectors if ky is conserved");
  (*SystemGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the cylinder", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new SingleStringOption ('\n', "confining-file", "file describing the confining potential");
  (*SystemGroup) += new SingleStringOption ('\n', "superconducting-file", "file describing the superconducting order paarameter");
  (*SystemGroup) += new SingleDoubleOption ('\n', "charging-energy", "factor in front of the charging energy (i.e 1/(2C))", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "average-nbrparticles", "average number of particles", 0.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderFermionsWithSpinTimeReversalSymmetryAndPairing -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int LzMax = Manager.GetInteger("lzmax");
  int TotalSz = Manager.GetInteger("total-sz");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int InitialLz = Manager.GetInteger("initial-ky");
  int NbrLz = Manager.GetInteger("nbr-ky");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  double Ratio = Manager.GetDouble("aspect-ratio");
  double Perimeter = Manager.GetDouble("cylinder-perimeter");
  if (Perimeter != 0.0)
    {
      Ratio = 2.0 * M_PI * (LzMax + 1) / (Perimeter * Perimeter);
    }
  bool FirstRun = true;

  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;
  double* OneBodyPotentialUpDown = 0;
  double* OneBodyPotentialPairing = 0;
  double** PseudoPotentials  = new double*[3];
  for (int i = 0; i < 3; ++i)
    {
      PseudoPotentials[i] = new double[LzMax + 1];
    }
  Complex** OffDiagonalOneBodyPotentialUpUp = 0;
  Complex** OffDiagonalOneBodyPotentialDownDown = 0;
  Complex* ComplexOneBodyPotentialPairing = 0;
  Complex** OffDiagonalComplexOneBodyPotentialPairing = 0;
  bool PreserveKySymmetryFlag = false;
  int MaximumMomentumTransfer = 0;

  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      if (FQHESphereSU2GetPseudopotentialsWithPairing(Manager.GetString("interaction-file"), LzMax, PseudoPotentials, 
						      OneBodyPotentialUpUp, OneBodyPotentialDownDown, OneBodyPotentialUpDown, OneBodyPotentialPairing) == false)
	{
	  return -1;
	}
    }
  if (OneBodyPotentialUpDown != 0)
    {
      cout << "warning, OneBodyPotentialUpDown is not supported" << endl;
    }

  if (Manager.GetString("confining-file") != 0)
    {
      MultiColumnASCIIFile ConfiningFile;
      if (ConfiningFile.Parse(Manager.GetString("confining-file")) == false)
	{
	  ConfiningFile.DumpErrors(cout);
	  return -1;
	}
      if ((ConfiningFile.GetNbrColumns() < 4) || (ConfiningFile.GetNbrLines() == 0))
	{
	  cout << "error, " << Manager.GetString("confining-file") << " has an invalid format" << endl;
	  return -1;
	}
      int* CreationIndices = ConfiningFile.GetAsIntegerArray(0);
      int* AnnihilationIndices = ConfiningFile.GetAsIntegerArray(1);
      double* TmpUpUpPotential = ConfiningFile.GetAsDoubleArray(2);
      double* TmpDownDownPotential = ConfiningFile.GetAsDoubleArray(3);
      double* TmpUpUpPhasePotential;
      if (ConfiningFile.GetNbrColumns() > 4)
	{
	  TmpUpUpPhasePotential = ConfiningFile.GetAsDoubleArray(4);
	}
      else
	{
	  TmpUpUpPhasePotential = new double[ConfiningFile.GetNbrLines()];
	  for (int i = 0; i < ConfiningFile.GetNbrLines(); ++i)
	    {
	      TmpUpUpPhasePotential[i] = 0.0;
	    }
	}
      double* TmpDownDownPhasePotential;
      if (ConfiningFile.GetNbrColumns() > 5)
	{
	  TmpDownDownPhasePotential = ConfiningFile.GetAsDoubleArray(5);
	}
      else
	{
	  TmpDownDownPhasePotential = new double[ConfiningFile.GetNbrLines()];
	  for (int i = 0; i < ConfiningFile.GetNbrLines(); ++i)
	    {
	      TmpDownDownPhasePotential[i] = 0.0;
	    }
	}
      for (int i = 0; i < ConfiningFile.GetNbrLines(); ++i)
	{
	  if ((CreationIndices[i] < 0) || (AnnihilationIndices[i] < 0) || (CreationIndices[i] > LzMax) || (AnnihilationIndices[i] > LzMax))
	    {
	      cout << "invalid indices line " << (i+1) << endl;
	      return -1;
	    }
	  if (abs(CreationIndices[i] - AnnihilationIndices[i]) > MaximumMomentumTransfer)
	    {
	      MaximumMomentumTransfer = abs(CreationIndices[i] - AnnihilationIndices[i]);
	    }
	}
      OneBodyPotentialUpUp = new double[LzMax + 1];
      OneBodyPotentialDownDown = new double[LzMax + 1];
      for (int i = 0; i <= LzMax; ++i)
	{
	  OneBodyPotentialUpUp[i] = 0.0;
	  OneBodyPotentialDownDown[i] = 0.0;	  
	}
      if (MaximumMomentumTransfer > 0)
	{
	  OffDiagonalOneBodyPotentialUpUp = new Complex*[LzMax + 1];
	  OffDiagonalOneBodyPotentialDownDown = new Complex*[LzMax + 1];
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      OffDiagonalOneBodyPotentialUpUp[i] = new Complex[MaximumMomentumTransfer];
	      OffDiagonalOneBodyPotentialDownDown[i] = new Complex[MaximumMomentumTransfer];
	      for (int j = 0; j < MaximumMomentumTransfer; ++j)	  
		{
		  OffDiagonalOneBodyPotentialUpUp[i][j] = 0.0;
		  OffDiagonalOneBodyPotentialDownDown[i][j] = 0.0;
		}
	    }	  
	}
      for (int i = 0; i < ConfiningFile.GetNbrLines(); ++i)
	{
	  int TmpMomentumTransfer = CreationIndices[i] - AnnihilationIndices[i];
	  if (TmpMomentumTransfer == 0)
	    {
	      OneBodyPotentialUpUp[CreationIndices[i]] = TmpUpUpPotential[i];
	      OneBodyPotentialDownDown[CreationIndices[i]] = TmpDownDownPotential[i];
	    }
	  else
	    {
	      if (TmpMomentumTransfer > 0)
		{
		  OffDiagonalOneBodyPotentialUpUp[CreationIndices[i]][TmpMomentumTransfer - 1] = (TmpUpUpPotential[i] * Phase(M_PI * TmpUpUpPhasePotential[i]));
		  OffDiagonalOneBodyPotentialDownDown[CreationIndices[i]][TmpMomentumTransfer - 1] = (TmpDownDownPotential[i] * Phase(M_PI * TmpDownDownPhasePotential[i]));
		}
	      else
		{
		  OffDiagonalOneBodyPotentialUpUp[CreationIndices[i]][-TmpMomentumTransfer - 1] = (TmpUpUpPotential[i] * Phase(-M_PI * TmpUpUpPhasePotential[i]));
		  OffDiagonalOneBodyPotentialDownDown[CreationIndices[i]][-TmpMomentumTransfer - 1] = (TmpDownDownPotential[i] * Phase(-M_PI * TmpDownDownPhasePotential[i]));
		}
	    }
	}
    }
  if (Manager.GetString("superconducting-file") != 0)
    {
      MultiColumnASCIIFile SuperconductingFile;
      if (SuperconductingFile.Parse(Manager.GetString("superconducting-file")) == false)
	{
	  SuperconductingFile.DumpErrors(cout);
	  return -1;
	}
      if ((SuperconductingFile.GetNbrColumns() < 3) || (SuperconductingFile.GetNbrLines() == 0))
	{
	  cout << "error, " << Manager.GetString("confining-file") << " has an invalid format" << endl;
	  return -1;
	}
      int* CreationIndices = SuperconductingFile.GetAsIntegerArray(0);
      int* AnnihilationIndices = SuperconductingFile.GetAsIntegerArray(1);
      double* TmpPairingPotential = SuperconductingFile.GetAsDoubleArray(2);
      double* TmpPairingPhasePotential;
      if (SuperconductingFile.GetNbrColumns() > 3)
	{
	  TmpPairingPhasePotential = SuperconductingFile.GetAsDoubleArray(3);
	}
      else
	{
	  TmpPairingPhasePotential = new double[SuperconductingFile.GetNbrLines()];
	  for (int i = 0; i < SuperconductingFile.GetNbrLines(); ++i)
	    {
	      TmpPairingPhasePotential[i] = 0.0;
	    }
	}
      for (int i = 0; i < SuperconductingFile.GetNbrLines(); ++i)
	{
	  if ((CreationIndices[i] < 0) || (AnnihilationIndices[i] < 0) || (CreationIndices[i] > LzMax) || (AnnihilationIndices[i] > LzMax))
	    {
	      cout << "invalid indices line " << (i+1) << endl;
	      return -1;
	    }
	  if (abs(CreationIndices[i] - AnnihilationIndices[i]) > MaximumMomentumTransfer)
	    {
	      MaximumMomentumTransfer = abs(CreationIndices[i] - AnnihilationIndices[i]);
	    }
	}
      ComplexOneBodyPotentialPairing = new Complex[LzMax + 1];
      for (int i = 0; i <= LzMax; ++i)
	{
	  ComplexOneBodyPotentialPairing[i] = 0.0;
	}
      if (MaximumMomentumTransfer > 0)
	{
	  OffDiagonalComplexOneBodyPotentialPairing = new Complex*[LzMax + 1];
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      OffDiagonalComplexOneBodyPotentialPairing[i] = new Complex[(2 * MaximumMomentumTransfer) + 1];
	      for (int j = 0; j < MaximumMomentumTransfer; ++j)	  
		{
		  OffDiagonalComplexOneBodyPotentialPairing[i][j] = 0.0;
		}
	    }	  
	}
      for (int i = 0; i < SuperconductingFile.GetNbrLines(); ++i)
	{
	  int TmpMomentumTransfer = CreationIndices[i] - AnnihilationIndices[i];
	  if (TmpMomentumTransfer == 0)
	    {
	      ComplexOneBodyPotentialPairing[CreationIndices[i]] = TmpPairingPotential[i] * Phase(M_PI * TmpPairingPhasePotential[i]);
	    }
	  else
	    {
	      OffDiagonalComplexOneBodyPotentialPairing[CreationIndices[i]][MaximumMomentumTransfer + TmpMomentumTransfer] = TmpPairingPotential[i] * Phase(M_PI * TmpPairingPhasePotential[i]);
	    }
	}	  
    }

  if (MaximumMomentumTransfer == 0)
    {
      PreserveKySymmetryFlag = true;
    }


  char* GeometryName = new char[128];
  if (Perimeter > 0.0)	
    {
      sprintf (GeometryName, "cylinder_perimeter_%.6f", Perimeter);
    }
  else
    {
      sprintf (GeometryName, "cylinder_ratio_%.6f", Ratio);
    }
  char* OutputFileName = new char [256 + strlen(Manager.GetString("interaction-name")) + strlen(GeometryName)];
  if (Manager.GetInteger("nbr-particles") >= 0)
    {
      if (PreserveKySymmetryFlag == true)
	{
	  sprintf (OutputFileName, "fermions_%s_su2_%s_cenergy_%.6f_n0_%.6f_pairing_n_%ld_2s_%d_sz_%d_lz", GeometryName, 
		   Manager.GetString("interaction-name"), 
		   Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), 
		   Manager.GetInteger("nbr-particles"), LzMax, TotalSz);
	}
      else
	{
	  sprintf (OutputFileName, "fermions_%s_su2_%s_cenergy_%.6f_n0_%.6f_pairing_n_%ld_2s_%d_sz_%d", GeometryName, 
		   Manager.GetString("interaction-name"), 
		   Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), 
		   Manager.GetInteger("nbr-particles"), LzMax, TotalSz);
	}
    }
  else
    {
      if (PreserveKySymmetryFlag == true)
	{
	  sprintf (OutputFileName, "fermions_%s_su2_%s_cenergy_%.6f_n0_%.6f_pairing_n_0_2s_%d_sz_%d_lz", GeometryName, 
		   Manager.GetString("interaction-name"), 
		   Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), LzMax, TotalSz);
	}
      else
	{
	  sprintf (OutputFileName, "fermions_%s_su2_%s_cenergy_%.6f_n0_%.6f_pairing_n_0_2s_%d_sz_%d.", GeometryName, 
		   Manager.GetString("interaction-name"), 
		   Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), LzMax, TotalSz);
	}
    }
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  int MinNbrParticles = abs(TotalSz);
  int MaxNbrParticles = (2 * (LzMax + 1)) - abs(TotalSz);
  int MaxL = 0;
  if (Manager.GetInteger("nbr-particles") >= 0)
    {
      MinNbrParticles = Manager.GetInteger("nbr-particles");
      MaxNbrParticles = MinNbrParticles;
    }

  for (int TmpNbrParticles = MinNbrParticles; TmpNbrParticles <= MaxNbrParticles; TmpNbrParticles += 2)
    {
      int NbrUp = (TmpNbrParticles + TotalSz) / 2;
      int NbrDown = (TmpNbrParticles - TotalSz) / 2;
      if ((NbrUp <= (LzMax + 1)) && (NbrUp >= 0) && (NbrDown <= (LzMax + 1)) && (NbrDown >= 0))
	{
	  int TmpMaxL = (((LzMax - NbrUp + 1) * NbrUp) + ((LzMax - NbrDown + 1) * NbrDown));
	  if (TmpMaxL > MaxL)
	    MaxL = TmpMaxL;
	}
    }

  if (PreserveKySymmetryFlag == true)
    {
      Lanczos.SetRealAlgorithms();
      char* CommentLine = new char [512];
      sprintf (CommentLine, " Lz E");
      int  L = 0;
      if (InitialLz >= 0)
	{
	  L = InitialLz;
	  if ((abs(MaxL) & 1) != 0)
	    L |= 1;
	  else
	    L &= ~0x1;
	}
      if (NbrLz > 0)
	{
	  if (L + (2 * (NbrLz - 1)) < MaxL)
	    MaxL = L + (2 * (NbrLz - 1));
	}
      else
	{
	  if (Manager.GetBoolean("force-negativeky") == true)
	    {
	      L = -MaxL;
	    }
	}
      
      for (; L <= MaxL; L += 2)
	{
	  ParticleOnSphereWithSpin* Space;
	  if (Manager.GetInteger("nbr-particles") >= 0)
	    {
	      Space = new FermionOnSphereWithSpin(Manager.GetInteger("nbr-particles"), L, LzMax, TotalSz);
	    }
	  else
	    {
	      Space = new FermionOnSphereWithSpinAndPairing(L, LzMax, TotalSz);
	    }
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian = 0;
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Hamiltonian = new ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairing (Space, LzMax, Ratio,
													 PseudoPotentials, OneBodyPotentialUpUp, 
													 OneBodyPotentialDownDown,
													 OneBodyPotentialPairing,
													 Manager.GetDouble("charging-energy"), 
													 Manager.GetDouble("average-nbrparticles"),
													 Architecture.GetArchitecture(), 
													 Memory, DiskCacheFlag,
													 LoadPrecalculationFileName);
	  double Shift = - 10.0;
	  Hamiltonian->ShiftHamiltonian(Shift);
	  char* EigenvectorName = 0;
	  if (Manager.GetBoolean("eigenstate") == true)	
	    {
	      EigenvectorName = new char [256 + strlen(OutputFileName)];
	      sprintf (EigenvectorName, "%s_lz_%d", OutputFileName, L);
	    }
	  
	  char* TmpString = new char[64];
	  sprintf (TmpString, "%d ", L);
 	  GenericRealMainTask Task(&Manager, Space, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, EigenvectorName);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  if (EigenvectorName != 0)
	    {
	      delete[] EigenvectorName;
	    }
	  delete Hamiltonian;
	  if (FirstRun == true)
	    FirstRun = false;
	  delete Space;
	}
    }
  else
    {
      char* CommentLine = new char [512];
      sprintf (CommentLine, " E");
      ParticleOnSphereWithSpin* Space;
      Space = new FermionOnSphereWithSpinAndPairingAllLz(LzMax, TotalSz);
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();

      double Shift = - 10.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{	  
	  EigenvectorName = new char [256 + strlen(OutputFileName)];
	  sprintf (EigenvectorName, "%s", OutputFileName);
	}
      
      char* TmpString = new char[64];
      sprintf (TmpString, "");
      GenericComplexMainTask Task(&Manager, Space, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
				  FirstRun, EigenvectorName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      delete Hamiltonian;
      delete Space;      
    }
  return 0;
}
