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

// build the binary mask for for a Hamiltonian Z term  (unsigned long long version)
// 
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = binary mask for for a Hamiltonian Z term 
ULONGLONG BuildHamiltonianZTermMaskLongLong (int x, int y, int z, int nbrSiteX, int nbrSiteY, int nbrSiteZ);

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

// get the parity of spin 0 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)
unsigned long GetSpin0Parity(unsigned long state);

// get the parity of spin 1 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)
unsigned long GetSpin1Parity(unsigned long state);

// get the parity of spin 0 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)
ULONGLONG GetSpin0Parity(ULONGLONG state);

// get the parity of spin 1 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)
ULONGLONG GetSpin1Parity(ULONGLONG state);


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
  double* FullReducedDensityMatrixEigenvalues;
  long FullReducedDensityMatrixNbrEigenvalues;


  if (TotalNbrSpins <= 64)
    {
      unsigned long* GroundState = new unsigned long[GroundStateDimension];
      
      GroundState[0l] = 0x0ul;
      GroundStateDimension = 1l;
      
      for (int x = 0; x < NbrSitesX; ++x)
	{
	  for (int y = 0; y < NbrSitesY; ++y)
	    {
	      for (int z = 0; z < NbrSitesZ; ++z)
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
      SortArrayUpOrderingAndRemoveDuplicates(GroundState, GroundStateDimension);
      
      
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
      unsigned long* RegionBConfigurations = new unsigned long[GroundStateDimension];
      int* RegionAIndices = new int[GroundStateDimension];
      FullReducedDensityMatrixEigenvalues = new double[GroundStateDimension];
      FullReducedDensityMatrixNbrEigenvalues = 0;
      unsigned long* RegionAHilbertSpaceFixedParities = new unsigned long[RegionAHilbertSpaceDimension];
      
      for (unsigned long TmpParity = 0x0ul; TmpParity <= 0x3ul; ++TmpParity)
	{
	  long RegionAHilbertSpaceDimensionFixedParities = 0l;
	  for  (long i = 0l; i < RegionAHilbertSpaceDimension; ++i)
	    {
	      unsigned long Tmp = RegionAHilbertSpace[i];
	      if (((GetSpin1Parity(Tmp) << 1) | GetSpin0Parity(Tmp)) == TmpParity)
		{
		  RegionAHilbertSpaceFixedParities[RegionAHilbertSpaceDimensionFixedParities] = Tmp;
		  ++RegionAHilbertSpaceDimensionFixedParities;
		}
	    }
	  
	  if (RegionAHilbertSpaceDimensionFixedParities > 0)
	    {
	      cout << "Building the entanglement matrix for the parity sector " << TmpParity << endl;
	      long RegionBHilbertSpaceDimensionFixedParities = 0l;
	      for (long i = 0l; i < GroundStateDimension; ++i)
		{
		  int Tmp = SearchInArray(GroundState[i] & RegionAMask, RegionAHilbertSpaceFixedParities, RegionAHilbertSpaceDimensionFixedParities);
		  if (Tmp >= 0)
		    {
		      RegionBConfigurations[RegionBHilbertSpaceDimensionFixedParities] = GroundState[i] & RegionBMask;
		      RegionAIndices[RegionBHilbertSpaceDimensionFixedParities] = Tmp;
		      ++RegionBHilbertSpaceDimensionFixedParities;
		    }
		}
	      SortArrayDownOrdering(RegionBConfigurations, RegionAIndices, RegionBHilbertSpaceDimensionFixedParities);
	      
	      cout << "Building the reduced density matrix" << endl;
	      RealSymmetricMatrix ReducedDensityMatrix (RegionAHilbertSpaceDimensionFixedParities, true);
	      long TmpIndex = 0l;
	      while (TmpIndex < RegionBHilbertSpaceDimensionFixedParities)
		{
		  long TmpIndex2 = TmpIndex + 1l;
		  while ((TmpIndex2 < RegionBHilbertSpaceDimensionFixedParities) && (RegionBConfigurations[TmpIndex] == RegionBConfigurations[TmpIndex2]))
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
	      cout << "Diagonalizing the reduced density matrix (" << RegionAHilbertSpaceDimensionFixedParities << "x" << RegionAHilbertSpaceDimensionFixedParities << ")" << endl;
	      RealDiagonalMatrix ReducedDensityMatrixEigenvalues(ReducedDensityMatrix.GetNbrRow(), true);
	      ReducedDensityMatrix.LapackDiagonalize(ReducedDensityMatrixEigenvalues);
	      for (long i = 0l; i < ReducedDensityMatrixEigenvalues.GetNbrRow(); ++i)
		{
		  FullReducedDensityMatrixEigenvalues[FullReducedDensityMatrixNbrEigenvalues] = ReducedDensityMatrixEigenvalues[i];
		  ++FullReducedDensityMatrixNbrEigenvalues;
		}
	    }
	}
      delete[] RegionBConfigurations;
      delete[] RegionAIndices;
      delete[] GroundState;
    }      
  else
    {
      timeval StartingTime;
      timeval EndingTime;
      gettimeofday (&(StartingTime), 0);
      ULONGLONG* GroundState = new ULONGLONG[GroundStateDimension];
      
      GroundState[0l] = ((ULONGLONG) 0x0ul);
      GroundStateDimension = 1l;
      
      for (int x = 0; x < NbrSitesX; ++x)
	{
	  for (int y = 0; y < NbrSitesY; ++y)
	    {
	      for (int z = 0; z < NbrSitesZ; ++z)
		{
		  ULONGLONG TmpMask = BuildHamiltonianZTermMaskLongLong(x, y, z, NbrSitesX, NbrSitesY, NbrSitesZ);
		  for (long i = 0l; i < GroundStateDimension; ++i)
		    {
		      GroundState[i + GroundStateDimension] = GroundState[i] ^ TmpMask;
		    }
		  GroundStateDimension *= 2l;
		}	  
	    }      
	}
      gettimeofday (&(EndingTime), 0);
      double DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			   ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
      cout << GroundStateDimension << " intermediate components generated for the groundstate (done in " << DeltaTime << "s)" << endl;
      gettimeofday (&(StartingTime), 0);
      SortArrayUpOrderingAndRemoveDuplicates(GroundState, GroundStateDimension);
      gettimeofday (&(EndingTime), 0);
      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
		    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
      cout << GroundStateDimension << " distinct components generated for the groundstate (done in " << DeltaTime << "s)" << endl;

      double GroundStateNormalizationCoefficient = 1.0 / sqrt((double) GroundStateDimension);
      double GroundStateSqrNormalizationCoefficient = 1.0 / ((double) GroundStateDimension);
      
      gettimeofday (&(StartingTime), 0);
      int RegionANbrSites = NbrSitesXA * NbrSitesYA * NbrSitesZA;
      ULONGLONG RegionAMask = ((ULONGLONG) 0x0ul);  
      for (int x = 0; x < NbrSitesXA; ++x)
	{
	  for (int y = 0; y < NbrSitesYA; ++y)
	    {
	      for (int z = 0; z < NbrSitesZA; ++z)
		{
		  RegionAMask |= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x, y, z, 0, NbrSitesX, NbrSitesY, NbrSitesZ);
		  RegionAMask |= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x, y, z, 1, NbrSitesX, NbrSitesY, NbrSitesZ);
		}
	    }
	}
      ULONGLONG RegionBMask = ((((ULONGLONG) 0x1ul) << (2 * TotalNbrSites)) - ((ULONGLONG) 0x1ul)) & (~RegionAMask);
      
      ULONGLONG* RegionAHilbertSpace = new ULONGLONG[GroundStateDimension];
      long RegionAHilbertSpaceDimension = 0l;
      for (long i = 0l; i < GroundStateDimension; ++i)
	{
	  RegionAHilbertSpaceDimension += SearchInSortedArrayAndInsert(GroundState[i] & RegionAMask, RegionAHilbertSpace, RegionAHilbertSpaceDimension);
	}
      ULONGLONG* TmpArray = new ULONGLONG[RegionAHilbertSpaceDimension]; 
      for (long i = 0l; i < RegionAHilbertSpaceDimension; ++i)
	{
	  TmpArray[i] = RegionAHilbertSpace[i];
	}
      delete[] RegionAHilbertSpace;
      RegionAHilbertSpace = TmpArray;
      gettimeofday (&(EndingTime), 0);
      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
		    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
      cout << "Hilbert space dimension for the A region = " << RegionAHilbertSpaceDimension << " (done in " << DeltaTime << "s)" << endl;
      ULONGLONG* RegionBConfigurations = new ULONGLONG[GroundStateDimension];
      int* RegionAIndices = new int[GroundStateDimension];
      FullReducedDensityMatrixEigenvalues = new double[RegionAHilbertSpaceDimension];
      FullReducedDensityMatrixNbrEigenvalues = 0;
      ULONGLONG* RegionAHilbertSpaceFixedParities = new ULONGLONG[RegionAHilbertSpaceDimension];
      
      for (ULONGLONG TmpParity = ((ULONGLONG) 0x0ul); TmpParity <= ((ULONGLONG) 0x3ul); ++TmpParity)
	{
	  gettimeofday (&(StartingTime), 0);
	  long RegionAHilbertSpaceDimensionFixedParities = 0l;
	  for  (long i = 0l; i < RegionAHilbertSpaceDimension; ++i)
	    {
	      ULONGLONG Tmp = RegionAHilbertSpace[i];
	      if (((GetSpin1Parity(Tmp) << 1) | GetSpin0Parity(Tmp)) == TmpParity)
		{
		  RegionAHilbertSpaceFixedParities[RegionAHilbertSpaceDimensionFixedParities] = Tmp;
		  ++RegionAHilbertSpaceDimensionFixedParities;
		}
	    }
	  gettimeofday (&(EndingTime), 0);	      
	  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	  cout << "Hilbert space dimension for the region A in  the parity sector = " << ((unsigned long) TmpParity) << " (done in " << DeltaTime << "s)"  << endl;
	  
	  if (RegionAHilbertSpaceDimensionFixedParities > 0)
	    {
	      cout << "Building the entanglement matrix" << endl;
	      gettimeofday (&(StartingTime), 0);
	      long RegionBHilbertSpaceDimensionFixedParities = 0l;
	      for (long i = 0l; i < GroundStateDimension; ++i)
		{
		  int Tmp = SearchInArray(GroundState[i] & RegionAMask, RegionAHilbertSpaceFixedParities, RegionAHilbertSpaceDimensionFixedParities);
		  if (Tmp >= 0)
		    {
		      RegionBConfigurations[RegionBHilbertSpaceDimensionFixedParities] = GroundState[i] & RegionBMask;
		      RegionAIndices[RegionBHilbertSpaceDimensionFixedParities] = Tmp;
		      ++RegionBHilbertSpaceDimensionFixedParities;
		    }
		}
	      SortArrayDownOrdering(RegionBConfigurations, RegionAIndices, RegionBHilbertSpaceDimensionFixedParities);
	      gettimeofday (&(EndingTime), 0);	      
	      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	      cout << "done in " << DeltaTime << "s"  << endl;
	      gettimeofday (&(StartingTime), 0);
	      cout << "Building the reduced density matrix" << endl;

	      RealSymmetricMatrix ReducedDensityMatrix (RegionAHilbertSpaceDimensionFixedParities, true);
	      long TmpIndex = 0l;
	      while (TmpIndex < RegionBHilbertSpaceDimensionFixedParities)
		{
		  long TmpIndex2 = TmpIndex + 1l;
		  while ((TmpIndex2 < RegionBHilbertSpaceDimensionFixedParities) && (RegionBConfigurations[TmpIndex] == RegionBConfigurations[TmpIndex2]))
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
	      gettimeofday (&(EndingTime), 0);	      
	      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	      cout << "done in " << DeltaTime << "s"  << endl;
	      gettimeofday (&(StartingTime), 0);
	      cout << "Diagonalizing the reduced density matrix (" << RegionAHilbertSpaceDimensionFixedParities << "x" << RegionAHilbertSpaceDimensionFixedParities << ")" << endl;
	      RealDiagonalMatrix ReducedDensityMatrixEigenvalues(ReducedDensityMatrix.GetNbrRow(), true);
	      ReducedDensityMatrix.LapackDiagonalize(ReducedDensityMatrixEigenvalues);
	      gettimeofday (&(EndingTime), 0);	      
	      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	      cout << "done in " << DeltaTime << "s"  << endl;
	      for (long i = 0l; i < ReducedDensityMatrixEigenvalues.GetNbrRow(); ++i)
		{
		  FullReducedDensityMatrixEigenvalues[FullReducedDensityMatrixNbrEigenvalues] = ReducedDensityMatrixEigenvalues[i];
		  ++FullReducedDensityMatrixNbrEigenvalues;
		}
	    }
	}
      delete[] RegionBConfigurations;
      delete[] RegionAIndices;
      delete[] GroundState;
    }

  double ReducedDensityMatrixTrace = 0.0;
  double ReducedDensityMatrixEntanglementEntropy = 0.0;
  long ReducedDensityMatrixNbrNonZeroEigenvalues = 0l;
  for (long i = 0l; i < FullReducedDensityMatrixNbrEigenvalues; ++i)
    {
      ReducedDensityMatrixTrace += FullReducedDensityMatrixEigenvalues[i];
      if (FullReducedDensityMatrixEigenvalues[i] > 0.0)
	{
	  ReducedDensityMatrixEntanglementEntropy -= FullReducedDensityMatrixEigenvalues[i] * log(FullReducedDensityMatrixEigenvalues[i]);
	  ReducedDensityMatrixNbrNonZeroEigenvalues++;
	}
    }
  cout << "Trace of the reduced density matrix = " << ReducedDensityMatrixTrace << endl;
  cout << "Number of non zero eigenvalues for the reduced density matrix = " << ReducedDensityMatrixNbrNonZeroEigenvalues << endl;
  cout << "Entangement entropy = " << (ReducedDensityMatrixEntanglementEntropy / log(2.0)) << " * log 2" << endl;
   

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

inline int GetHaahCodeLinearizedIndex (int x, int y, int z, int spinIndex, int nbrSiteX, int nbrSiteY, int nbrSiteZ)
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
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x,     y,     z,     1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z,     0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y + 1, z,     1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x,     y,     z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z + 1, 1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x,     y + 1, z + 1, 1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y + 1, z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  return TmpMask;
}

// build the binary mask for for a Hamiltonian Z term  (unsigned long long version)
// 
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = binary mask for for a Hamiltonian Z term 

ULONGLONG BuildHamiltonianZTermMaskLongLong (int x, int y, int z, int nbrSiteX, int nbrSiteY, int nbrSiteZ)
{
  ULONGLONG TmpMask = 0x0ul;
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x,     y,     z,     1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z,     0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y + 1, z,     1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x,     y,     z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z + 1, 1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x,     y + 1, z + 1, 1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y + 1, z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  return TmpMask;
}

// get the parity of spin 0 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)

inline unsigned long GetSpin0Parity(unsigned long state)
{
  state &= 0x5555555555555555ul;
  state ^= state >> 32;
  state ^= state >> 16;
  state ^= state >> 8;
  state ^= state >> 4;
  state ^= state >> 2;
  state ^= state >> 1;
  return (state &0x1ul);
}

// get the parity of spin 1 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)

inline unsigned long GetSpin1Parity(unsigned long state)
{
  state &= 0xaaaaaaaaaaaaaaaaul;
  state ^= state >> 32;
  state ^= state >> 16;
  state ^= state >> 8;
  state ^= state >> 4;
  state ^= state >> 2;
  state ^= state >> 1;
  return (state &0x1ul);
}


// get the parity of spin 0 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)

inline ULONGLONG GetSpin0Parity(ULONGLONG state)
{
  state &= (((ULONGLONG) 0x5555555555555555ul) << 32) | ((ULONGLONG) 0x5555555555555555ul);
  state ^= state >> 64;
  state ^= state >> 32;
  state ^= state >> 16;
  state ^= state >> 8;
  state ^= state >> 4;
  state ^= state >> 2;
  state ^= state >> 1;
  return (state &0x1ul);
}

// get the parity of spin 1 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)

inline ULONGLONG GetSpin1Parity(ULONGLONG state)
{
  state &= (((ULONGLONG) 0xaaaaaaaaaaaaaaaaul) << 32) | ((ULONGLONG) 0xaaaaaaaaaaaaaaaaul);
  state ^= state >> 64;
  state ^= state >> 32;
  state ^= state >> 16;
  state ^= state >> 8;
  state ^= state >> 4;
  state ^= state >> 2;
  state ^= state >> 1;
  return (state &0x1ul);
}

