#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Tools/FQHESpectrum/QHEOnSphereLzSortedSpectrum.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


// subdivide a basis into different L orthonormal eigenstate
//
// vectors = reference on the matrix whose columns are the vectors that span the basis
// oper = pointer to operator whivh allow to evaluate total L matrix elements
// totalMaxLz = maximum L value that can be reach by the system (-1/2 if L is half integer)
// subspaceSize = reference on the array that will be filled  with the dimension of each fixed L subspace (index equal to L if L is integer, L-1/2 if L is half integer)
// subspacePositions = reference on the array that will be filled  with the position of the first occurence of a vector with a given fixed L subspace 
//                     (index equal to L if L is integer, L-1/2 if L is half integer)
int ReshuffleVectors (RealVector* vectors, int nbrVectors, double componentError);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereExtractRootPartitions" , "0.01");
  OptionGroup* MainGroup = new OptionGroup ("main options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += MiscGroup;
  Manager += MainGroup;
 
  (*MainGroup) += new SingleStringOption  ('\0', "input-file", "name of the file which contains the spectrum");
  (*MainGroup) += new SingleDoubleOption  ('\n', "energy-error", "energy below which a state is considered as a zero energy states", 1e-10);
  (*MainGroup) += new SingleDoubleOption  ('\n', "component-error", "value below which norm of a vector component is assumed to be zero", 1e-8);
  (*MainGroup) += new SingleIntegerOption  ('\n', "output-precision", "numerical display precision", 14, true, 2, true, 14);
  (*MainGroup) += new SingleStringOption  ('o', "output", "output name for the root partition list (default name uses input-file, replacing dat extension with root)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereExtractRootPartitions -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(Manager.GetString("input-file") == 0)
    {
      cout << "no input spectrum" << endl << "see man page for option syntax or type FQHESphereExtractRootPartitions -h" << endl;
      return -1;
    }

  double EnergyError = Manager.GetDouble("energy-error");
  double VectorError = Manager.GetDouble("component-error");
  int NbrFluxQuanta = 0;
  int NbrParticles = 0;
  bool FermionFlag = false;
  bool SU2SpinFlag = false;
  bool SU3SpinFlag = false;
  bool SU4SpinFlag = false;
  int TotalSz = 0;
  int TotalTz = 0;
  int TotalY = 0;
  int TotalIz = 0;
  int TotalPz = 0;
  if (FQHEOnSphereFindSystemInfoFromFileName(Manager.GetString("input-file"), NbrParticles, NbrFluxQuanta, FermionFlag) == false)
    {
      cout << "can't retrieve system information form file name " << Manager.GetString("input-file") << endl;
	return -1;
    }
  FQHEOnSphereFindInternalSymmetryGroupFromFileName(Manager.GetString("input-file"), SU2SpinFlag, SU3SpinFlag, SU4SpinFlag);

  FQHEOnSphereLzSortedSpectrum Spectrum (Manager.GetString("input-file"), EnergyError);
  if (Spectrum.IsSpectrumValid() == false)
    {
      cout << "Spectrum " << Manager.GetString("input-file") << " is unreadable or is not valid" << endl;
      return -1;           
    }
  
  int TotalMaxLz = Spectrum.GetMaxLzValue();
  int MaxLzValue = 0;
  if ((TotalMaxLz & 1) != 0)
    MaxLzValue = 1;
  int* LzDegeneracy = new int[TotalMaxLz + 1];
  long TotalNbrZeroEnergyStates = 0l;
  while ((MaxLzValue <= TotalMaxLz) && (fabs(Spectrum.GetEnergy(MaxLzValue, 0)) < EnergyError))
    {
//      cout << MaxLzValue << " " << Spectrum.GetEnergy(MaxLzValue, 0) << endl;
      LzDegeneracy[MaxLzValue / 2] = Spectrum.GetDegeneracy(MaxLzValue, 0);
      if (MaxLzValue == 0)
	TotalNbrZeroEnergyStates += LzDegeneracy[MaxLzValue / 2];
      else
	TotalNbrZeroEnergyStates += 2l * LzDegeneracy[MaxLzValue / 2];
      MaxLzValue += 2;
    }
  MaxLzValue -= 2;
  if (TotalNbrZeroEnergyStates == 0l)
    {
      cout << "Spectrum " << Manager.GetString("input-file") << " does not contain any zero energy state" << endl;      
      return -1;
    }

  char* OutputFileName = Manager.GetString("output");
  if (OutputFileName == 0)
    {
      OutputFileName = ReplaceExtensionToFileName (Manager.GetString("input-file"), "dat" , "root");
    }
  ofstream File;
  File.open(OutputFileName, ios::out);
  cout << "total nbr of root partitions : " << TotalNbrZeroEnergyStates << endl;
  cout << "nbr of root partitions per Lz sector : " << endl;;
  File << "total nbr of root partitions : " << TotalNbrZeroEnergyStates << endl;
  File << "nbr of root partitions per Lz sector : " << endl;

  for (int i = (MaxLzValue & 1); i <= MaxLzValue; i += 2)
    {
      cout << LzDegeneracy[i >> 1] << " ";
      File << LzDegeneracy[i >> 1] << " ";
    }
  File << endl << endl << "-----------------------------------------" << endl;
  cout << endl << endl << "-----------------------------------------" << endl;

  char* BaseFileName = RemoveExtensionFromFileName(Manager.GetString("input-file"), "dat");
  char* TmpFileName = new char [strlen(BaseFileName) + 256]; 
  BaseFileName[strlen(BaseFileName) - 1] = '\0';
  

  for (int TotalLz = (MaxLzValue & 1); TotalLz <= MaxLzValue; TotalLz += 2)
    {
      if ((TotalMaxLz & 1) != 0)
	{
	  File << "Lz = " << TotalLz << "/2 : " << endl;
	}
      else
	{
	  File << "Lz = " << (TotalLz  >> 1) << " : " << endl;
	}      
      int TmpNbrStates = LzDegeneracy[TotalLz >> 1];
      RealVector* TmpVectors = new RealVector[TmpNbrStates];
      for (int j = 0; j < TmpNbrStates; ++j)
	{
	  sprintf (TmpFileName, "%s_%d.%d.vec", BaseFileName, TotalLz, j);
 	  if (TmpVectors[j].ReadVector(TmpFileName) == false)
	    {
	      cout << "error while reading " << TmpFileName << endl;
	      return -1;
	    }
	}
      int* RootPositions = new int [TmpNbrStates];
      for (int j = 0; j < TmpNbrStates; ++j)
	RootPositions[j] = 0;

      int MinTmpPos = ReshuffleVectors(TmpVectors, TmpNbrStates, VectorError);
      RootPositions[0] = MinTmpPos;

      for (int k = 1; k < TmpNbrStates; ++k)
	{
	  for (int j = k; j < TmpNbrStates; ++j)
	    TmpVectors[j].AddLinearCombination(-TmpVectors[j][MinTmpPos] / TmpVectors[k - 1][MinTmpPos], TmpVectors[k -1]);
	  MinTmpPos = ReshuffleVectors(TmpVectors + k, TmpNbrStates - k, VectorError);
	  RootPositions[k] = MinTmpPos;
	}
      
      ParticleOnSphere* Space = 0;
      if (FermionFlag == false)
	{
	  if (SU2SpinFlag == false)
	    {
	      Space = new BosonOnSphereShort(NbrParticles, TotalLz, NbrFluxQuanta);
	    }
	  else
	    {
	      Space = new BosonOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
	    }
	}
      else
	{
	  if ((SU2SpinFlag == false) && (SU3SpinFlag == false) && (SU4SpinFlag == false))
	    {
#ifdef __64_BITS__
	      if (NbrFluxQuanta <= 63)
		Space = new FermionOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, NbrFluxQuanta);
#else
	      if (NbrFluxQuanta <= 31)
		Space = new FermionOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, NbrFluxQuanta);
#endif
	    }
	  else
	    if (SU2SpinFlag == true)
	      Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
	    else 
	      if (SU3SpinFlag == true)
		Space = new FermionOnSphereWithSU3Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalTz, TotalY);
	      else
		if (SU4SpinFlag == true)
		  Space = new FermionOnSphereWithSU4Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, TotalIz, TotalPz);	    
	}
      
     for (int j = 0; j < TmpNbrStates; ++j)
       Space->PrintState(File, RootPositions[j]) << endl;
      File << "nbr states = " << TmpNbrStates << endl;
      File << "-----------------------------------------" << endl;

      delete Space;
      delete[] RootPositions;
      delete[] TmpVectors;
    }
  File.close();

  return 0;
}


int ReshuffleVectors (RealVector* vectors, int nbrVectors, double componentError)
{
  int MinTmpPos = vectors[0].GetVectorDimension();
  int TmpVectorPos = 0;
  for (int j = 0; j < nbrVectors; ++j)      
    {
      int TmpPos = 0;
      while ((TmpPos < vectors[0].GetVectorDimension()) && (fabs(vectors[0][TmpPos]) < componentError))
	{
	  TmpPos++;
	}
      if (TmpPos < MinTmpPos)
	{
	  TmpVectorPos = j;
	  MinTmpPos = TmpPos;
	}
    }
  if (TmpVectorPos != 0)
    {
      RealVector TmpVector = vectors[TmpVectorPos];
      vectors[TmpVectorPos] = vectors[0];
       vectors[0] = TmpVector;
    }
  TmpVectorPos = 0;
  double MaxComponent = fabs(vectors[0][MinTmpPos]);
  for (int j = 1; j < nbrVectors; ++j)      
    {
      double TmpComponent = fabs(vectors[j][MinTmpPos]);
      if (TmpComponent > MaxComponent)
	{
	  TmpVectorPos = j;
	  MaxComponent = TmpComponent;
	}
    }
  if (TmpVectorPos != 0)
    {
      RealVector TmpVector = vectors[TmpVectorPos];
      vectors[TmpVectorPos] = vectors[0];
       vectors[0] = TmpVector;
    }
  return MinTmpPos;
}
