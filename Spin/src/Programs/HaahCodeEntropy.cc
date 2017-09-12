#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "GeneralTools/ArrayTools.h"

#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;


// build the binary mask for for a Hamiltonian Z term 
// 
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = binary mask for for a Hamiltonian Z term 
unsigned long BuildHamiltonianZTermMask (int x, int y, int z, int nbrSiteX, int nbrSiteY, int nbrSiteZ);

// build the linearized index from the spin coordinates
//
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// spinIndex = index of the spin we want to address (either 0 or 1)
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = linearized index
int GetHaahCodeLinearizedIndex (int x, int y, int z, int spinIndex, int nbrSiteX, int nbrSiteY, int nbrSiteZ);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("HaahCodeEntropy" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('z', "nbr-sitez", "number of sites along the z direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbra-sitex", "number of sites along the x direction for the region A", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbra-sitey", "number of sites along the y direction for the region A", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbra-sitez", "number of sites along the z direction for the region A", 2);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HaahCodeEntropy -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSitesX = Manager.GetInteger("nbr-sitex");
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSitesZ = Manager.GetInteger("nbr-sitez");
  int NbrSitesXA = Manager.GetInteger("nbra-sitex");
  int NbrSitesYA = Manager.GetInteger("nbra-sitey");
  int NbrSitesZA = Manager.GetInteger("nbra-sitez");

  int TotalNbrSites = NbrSitesX * NbrSitesY * NbrSitesZ;
  int TotalNbrSpins = 2 * TotalNbrSites;
  long GroundStateDimension = 2l << TotalNbrSites;

  unsigned long* GroundState = new unsigned long[GroundStateDimension];

  GroundState[0l] = 0x0ul;
  GroundStateDimension = 1l;

  // avoid the double counting of configurations when having PBC and only two sites
  int EffectiveNbrSitesX = NbrSitesX;
  if (EffectiveNbrSitesX == 2)
    {
      EffectiveNbrSitesX = 1;
    }
  int EffectiveNbrSitesY = NbrSitesY;
  if (EffectiveNbrSitesY == 2)
    {
      EffectiveNbrSitesY = 1;
    }
  int EffectiveNbrSitesZ = NbrSitesZ;
  if (EffectiveNbrSitesZ == 2)
    {
      EffectiveNbrSitesZ = 1;
    }
  for (int x = 0; x < EffectiveNbrSitesX; ++x)
    {
      for (int y = 0; y < EffectiveNbrSitesY; ++y)
	{
	  for (int z = 0; z < EffectiveNbrSitesZ; ++z)
	    {
	      unsigned long TmpMask = BuildHamiltonianZTermMask(x, y, z, NbrSitesX, NbrSitesY, NbrSitesZ);
	      for (long i = 0l; i < GroundStateDimension; ++i)
		{
		  GroundState[i + GroundStateDimension] = GroundState[i] ^ TmpMask;
		}
	      GroundStateDimension *= 2l;
	    }	  
	}      
    }
  double GroundStateNormalizationCoefficient = 1.0 / sqrt((double) GroundStateDimension);
  double GroundStateSqrNormalizationCoefficient = 1.0 / ((double) GroundStateDimension);
  cout << GroundStateDimension << " components generated for the groundstate" << endl;

  int RegionANbrSites = NbrSitesXA * NbrSitesYA * NbrSitesZA;
  unsigned long RegionAMask = 0x0ul;  
  for (int x = 0; x < NbrSitesXA; ++x)
    {
      for (int y = 0; y < NbrSitesYA; ++y)
	{
	  for (int z = 0; z < NbrSitesZA; ++z)
	    {
	      RegionAMask |= 0x1ul << GetHaahCodeLinearizedIndex(x, y, z, 0, NbrSitesX, NbrSitesY, NbrSitesZ);
	      RegionAMask |= 0x1ul << GetHaahCodeLinearizedIndex(x, y, z, 1, NbrSitesX, NbrSitesY, NbrSitesZ);
	    }
	}
    }
  unsigned long RegionBMask = ((0x1ul << (2 * TotalNbrSites)) - 0x1ul) & (~RegionAMask);

  unsigned long* RegionAHilbertSpace = new unsigned long[GroundStateDimension];
  long RegionAHilbertSpaceDimension = 0l;
  for (long i = 0l; i < GroundStateDimension; ++i)
    {
      RegionAHilbertSpaceDimension += SearchInSortedArrayAndInsert(GroundState[i] & RegionAMask, RegionAHilbertSpace, RegionAHilbertSpaceDimension);
    }
  cout << "Hilbert space dimension for the A region = " << RegionAHilbertSpaceDimension << endl;

  cout << "Building the entanglement matrix" << endl;
  unsigned long* RegionBConfigurations = new unsigned long[GroundStateDimension];
  int* RegionAIndices = new int[GroundStateDimension];
  for (long i = 0l; i < GroundStateDimension; ++i)
    {
      RegionBConfigurations[i] = GroundState[i] & RegionBMask;
      RegionAIndices[i] = SearchInArray(GroundState[i] & RegionAMask, RegionAHilbertSpace, RegionAHilbertSpaceDimension);
    }
  SortArrayDownOrdering(RegionBConfigurations, RegionAIndices, GroundStateDimension);

  cout << "Building the reduced density matrix" << endl;
  RealSymmetricMatrix ReducedDensityMatrix (RegionAHilbertSpaceDimension, true);
  long TmpIndex = 0l;
  while (TmpIndex < GroundStateDimension)
    {
      long TmpIndex2 = TmpIndex + 1l;
      while ((TmpIndex2 < GroundStateDimension) && (RegionBConfigurations[TmpIndex] == RegionBConfigurations[TmpIndex2]))
	{
	  ++TmpIndex2;
	}
      for (long i = TmpIndex; i < TmpIndex2; ++i)
	{
	  for (long j= TmpIndex; j < TmpIndex2; ++j)
	    {
	      if (RegionAIndices[i] <= RegionAIndices[j])
		{
		  ReducedDensityMatrix.AddToMatrixElement(RegionAIndices[i], RegionAIndices[j], GroundStateSqrNormalizationCoefficient);
		}
	    }
	}
      TmpIndex = TmpIndex2;
    }
  cout << "Diagonalizing the reduced density matrix" << endl;
  RealDiagonalMatrix ReducedDensityMatrixEigenvalues(RegionAHilbertSpaceDimension, true);
  ReducedDensityMatrix.LapackDiagonalize(ReducedDensityMatrixEigenvalues);

  double ReducedDensityMatrixTrace = 0.0;
  double ReducedDensityMatrixEntanglementEntropy = 0.0;
  long ReducedDensityMatrixNbrNonZeroEigenvalues = 0l;
  for (long i = 0l; i < RegionAHilbertSpaceDimension; ++i)
    {
      ReducedDensityMatrixTrace += ReducedDensityMatrixEigenvalues[i];
      if (ReducedDensityMatrixEigenvalues[i] > 0.0)
	{
	  ReducedDensityMatrixEntanglementEntropy -= ReducedDensityMatrixEigenvalues[i] * log(ReducedDensityMatrixEigenvalues[i]);
	  ReducedDensityMatrixNbrNonZeroEigenvalues++;
	}
    }
  cout << "Trace of the reduced density matrix = " << ReducedDensityMatrixTrace << endl;
  cout << "Number of non zero eigenvalues for the reduced density matrix = " << ReducedDensityMatrixNbrNonZeroEigenvalues << endl;
  cout << "Entangement entropy = " << (ReducedDensityMatrixEntanglementEntropy / log(2.0)) << " * log 2" << endl;
   
  delete[] GroundState;

  return 0;
}


// build the linearized index from the spin coordinates
//
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// spinIndex = index of the spin we want to address (either 0 or 1)
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = linearized index

int GetHaahCodeLinearizedIndex (int x, int y, int z, int spinIndex, int nbrSiteX, int nbrSiteY, int nbrSiteZ)
{
  return (spinIndex + 2 * ((z % nbrSiteZ) + nbrSiteZ * ((y % nbrSiteY) + nbrSiteY * (x % nbrSiteX))));  
}

// build the binary mask for for a Hamiltonian Z term 
// 
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = binary mask for for a Hamiltonian Z term 

unsigned long BuildHamiltonianZTermMask (int x, int y, int z, int nbrSiteX, int nbrSiteY, int nbrSiteZ)
{
  unsigned long TmpMask = 0x0ul;
  TmpMask |= 0x1ul << GetHaahCodeLinearizedIndex(x,     y,     z,     1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask |= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z,     0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask |= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y + 1, z,     1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask |= 0x1ul << GetHaahCodeLinearizedIndex(x,     y,     z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask |= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask |= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z + 1, 1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask |= 0x1ul << GetHaahCodeLinearizedIndex(x,     y + 1, z + 1, 1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask |= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y + 1, z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  return TmpMask;
}

