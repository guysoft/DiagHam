#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FTIRealSpaceEntanglementEntropyParticlePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*SystemGroup) += new BooleanOption  ('\n', "decoupled", "assume that the total spin is a good quantum number of the problem");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "particles have a SU(2) spin");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIRealSpaceEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
    
    
  int NbrSpaces = 1;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  ComplexVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  int* TotalKx = 0;
  int* TotalKy = 0;
  int NbrParticles = 0;
  int NbrSites = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  bool Statistics = true;
  double* Coefficients = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  int TotalSpin = 0;
  bool TwoDTranslationFlag = false;
  bool SU2SpinFlag = Manager.GetBoolean("su2-spin");
  bool GutzwillerFlag = false;

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FTIEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) && 
      (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  if ((Manager.GetString("degenerated-groundstate") != 0) && 
      (IsFile(Manager.GetString("degenerated-groundstate")) == false))
    {
      cout << "can't open file " << Manager.GetString("degenerated-groundstate") << endl;
      return -1;
    }

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalKx = new int[1];
      TotalKy = new int[1];
      Coefficients = new double[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));
      Coefficients[0] = 1.0;
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       NbrSpaces = DegeneratedFile.GetNbrLines();
       GroundStateFiles = new char* [NbrSpaces];
       TotalKx = new int[NbrSpaces];
       TotalKy = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));		   
	 }
       if (DegeneratedFile.GetNbrColumns() == 1)
	 {
	   Coefficients = new double[NbrSpaces];
	   for (int i = 0; i < NbrSpaces; ++i)
	     Coefficients[i] = 1.0 / ((double) NbrSpaces);
	 }
       else
	 {
	   double TmpSum = 0.0;
	   Coefficients = DegeneratedFile.GetAsDoubleArray(1);
	   for (int i = 0; i < NbrSpaces; ++i)
	     TmpSum += Coefficients[i];
	   TmpSum = 1.0 / TmpSum;
	   for (int i = 0; i < NbrSpaces; ++i)
	     Coefficients[i] *= TmpSum;
	 }
    }
    
    if (FTIHubbardModelFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, NbrSites, Statistics, GutzwillerFlag) == false)
      {
	cout << "error while retrieving system parameters from file name " << GroundStateFiles[0] << endl;
	return -1;
	
      }
    TwoDTranslationFlag = FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(GroundStateFiles[0],
											   NbrParticles, NbrSites, TotalKx[0], TotalKy[0], NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag);
    bool TotalSpinConservedFlag = FTIHubbardModelWithSzFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, NbrSites, TotalSpin, Statistics, GutzwillerFlag);
    
    if (TwoDTranslationFlag == true)
      { 
	for (int i = 0; i < NbrSpaces; ++i)
	  {
	    TotalKx[i] = 0;
	    TotalKy[i] = 0;
	    if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(GroundStateFiles[i],
										 NbrParticles, NbrSites, TotalKx[i], TotalKy[i], NbrSiteX, NbrSiteY, Statistics, GutzwillerFlag) == false)
	      {
		cout << "error while retrieving 2D translation parameters from file name " << GroundStateFiles[i] << endl;
		return -1;
	      }
	  } 
      }

  GroundStates = new ComplexVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i] << endl;
	  return -1;      
	}
    }
    
    if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      if (TwoDTranslationFlag == false)
	{
	  if ((Manager.GetBoolean("decoupled") == false) || (SU2SpinFlag == false))
	    {
	      DensityMatrixFile << "#  N    lambda";
	    }
	  else
	    {
	      DensityMatrixFile << "#  N    Sz    lambda";
	    }
	}
      else
	{
	  if ((Manager.GetBoolean("decoupled") == false) || (SU2SpinFlag == false))
	    {
	      DensityMatrixFile << "#  N    Kx    Ky    lambda";
	    }
	  else
	    {
	      DensityMatrixFile << "#  N    Kx    Ky    Sz    lambda";
	    }
      }
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }
  cout << "NbrParticles = " << NbrParticles << " NbrSites = "  << NbrSites << endl;

   
  
  int TotalNbrSites = NbrSiteX * NbrSiteY;
  if (TwoDTranslationFlag == false)
    TotalNbrSites = 1;
  int* NbrGroundStatePerMomentumSector = new int[TotalNbrSites];
  ComplexVector** GroundStatePerMomentumSector = new ComplexVector*[TotalNbrSites];
  double** CoefficientPerMomentumSector = new double*[TotalNbrSites];
  for (int i = 0; i < TotalNbrSites; ++i)
    {
      NbrGroundStatePerMomentumSector[i] = 0;
      GroundStatePerMomentumSector[i] = 0;
      CoefficientPerMomentumSector[i] = 0;
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i] << endl;
	  return -1;      
	}
      int TmpIndex;
      if(TwoDTranslationFlag == true)
	TmpIndex = (TotalKx[i] * NbrSiteY) + TotalKy[i];
      else
	TmpIndex = 0;
      NbrGroundStatePerMomentumSector[TmpIndex]++; 
    }
  for (int i = 0; i < TotalNbrSites; ++i)
    {
      if (NbrGroundStatePerMomentumSector[i] > 0)
	{
	  GroundStatePerMomentumSector[i] = new ComplexVector[NbrGroundStatePerMomentumSector[i]];
	  CoefficientPerMomentumSector[i] = new double[NbrGroundStatePerMomentumSector[i]];
	}
      NbrGroundStatePerMomentumSector[i] = 0;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex;
      if(TwoDTranslationFlag == true )
	TmpIndex = (TotalKx[i] * NbrSiteY) + TotalKy[i];
      else
	TmpIndex = 0;
      GroundStatePerMomentumSector[TmpIndex][NbrGroundStatePerMomentumSector[TmpIndex]] = GroundStates[i];
      CoefficientPerMomentumSector[TmpIndex][NbrGroundStatePerMomentumSector[TmpIndex]] = Coefficients[i];
      NbrGroundStatePerMomentumSector[TmpIndex]++;
    }  


  int MaxNbrSpaces;
  if (TwoDTranslationFlag == false)
    MaxNbrSpaces = 1;
  else
    MaxNbrSpaces = NbrSiteX * NbrSiteY;
  ParticleOnSphere** Spaces = new ParticleOnSphere*[MaxNbrSpaces];
  for (int i = 0; i < MaxNbrSpaces; ++i)
    {
      Spaces[i] = 0;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex;
      if (TwoDTranslationFlag == false)
	TmpIndex = 0;
      else
	TmpIndex = TotalKx[i] * NbrSiteY + TotalKy[i];
      
      if (Spaces[TmpIndex] == 0)
	{
	  if (Statistics == true)
	    {
	      if (SU2SpinFlag == false)
		{
		  if (TwoDTranslationFlag == false)
		    Spaces[TmpIndex] = new FermionOnLatticeRealSpace (NbrParticles, NbrSites);
		  else
		    Spaces[TmpIndex] = new FermionOnLatticeRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
		}
	      else
		{
		  if (GutzwillerFlag == false)
		    {
		      if (TotalSpinConservedFlag == false)
			{
			  if (TwoDTranslationFlag == false)
			    { 
			      Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
			    }
			  else
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
			}
		      else
			{
			  if (TwoDTranslationFlag == false)
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpace (NbrParticles, TotalSpin, NbrSites, 10000000);
			  else
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, TotalSpin, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
			}
		    }
		  else
		    {
		      if (TotalSpinConservedFlag == false)
			{
			  if (TwoDTranslationFlag == false)
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
			  else
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY);
			}
		      else
			{
			  if (TwoDTranslationFlag == false)
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, TotalSpin, NbrSites, 10000000);
			  else
			    Spaces[TmpIndex] = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, TotalSpin, NbrSites, TotalKx[i], NbrSiteX, TotalKy[i], NbrSiteY, 10000000);
			}
		      
		    }
		}
	    }
	  else
	    {
	      cout << " Bosonic statistics not implemented" << endl;
	    }
	}
      
      if (Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	  cout << Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() << " " << GroundStates[i].GetLargeVectorDimension() << endl;
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }
  
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);
  
  int MaxSubsystemNbrParticles = (NbrParticles >> 1) + (NbrParticles & 1);
  if (Manager.GetInteger("max-na") > 0)
    MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
  int SubsystemNbrParticles = Manager.GetInteger("min-na");
  
  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      int ComplementaryNbrParticles = NbrParticles - SubsystemNbrParticles;
      int MinSz = -SubsystemNbrParticles;
      if ((TotalSpin - ComplementaryNbrParticles) > MinSz)
	MinSz = (TotalSpin - ComplementaryNbrParticles);
      int MaxSz = SubsystemNbrParticles;
      if ((TotalSpin + ComplementaryNbrParticles) < MaxSz)
	MaxSz = (TotalSpin + ComplementaryNbrParticles);
      if ((Manager.GetBoolean("decoupled") == false) || (Manager.GetBoolean("su2-spin")) == false)
	{
	  MinSz = 0;
	  MaxSz = 0;
	}
      int SubsystemTotalKxMin = 0;
      int SubsystemTotalKyMin = 0;
      int SubsystemTotalKxMax = 1;
      int SubsystemTotalKyMax = 1;
      if(TwoDTranslationFlag == true)
      {
	SubsystemTotalKxMax = NbrSiteX;
	SubsystemTotalKyMax = NbrSiteY;
      }
      
      for (int SubsystemTotalSz = MinSz; SubsystemTotalSz <= MaxSz; SubsystemTotalSz += 2)
	{
	  for (int SubsystemTotalKx = SubsystemTotalKxMin; SubsystemTotalKx < SubsystemTotalKxMax; ++SubsystemTotalKx)
	    {
	      for (int SubsystemTotalKy = SubsystemTotalKyMin; SubsystemTotalKy < SubsystemTotalKyMax; ++SubsystemTotalKy)
		{
		  if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
		    {      
		      if(TwoDTranslationFlag == false)
			cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << endl;
		      else
			cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy << endl;
		    }
		  else
		    {      
		      if(TwoDTranslationFlag == false)
			cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " Sz = " << SubsystemTotalSz << endl;
		      else
			cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy 
			     << " Sz = " << SubsystemTotalSz << endl;
		    }
		  timeval TotalStartingTime;
		  timeval TotalEndingTime;
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
		  int TmpIndex = 0;
		  while (NbrGroundStatePerMomentumSector[TmpIndex] == 0)
		    ++TmpIndex;
		  HermitianMatrix PartialDensityMatrix;
		  if (Statistics == true)
		    {
		      if (SU2SpinFlag == false)
			{      
			  if (TwoDTranslationFlag == false)
			    {
			      PartialDensityMatrix = ((FermionOnLatticeRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				{
				  HermitianMatrix TmpMatrix = ((FermionOnLatticeRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				  PartialDensityMatrix += TmpMatrix;
				}
			    }
			  else
			    {
			      PartialDensityMatrix = ((FermionOnLatticeRealSpaceAnd2DTranslation*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
				{
				  HermitianMatrix TmpMatrix = ((FermionOnLatticeRealSpaceAnd2DTranslation*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
				  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
				  PartialDensityMatrix += TmpMatrix;
				}
			    }
			}
		      else
			{
			  if (GutzwillerFlag == false)
			    {
			      if (Manager.GetBoolean("decoupled") == false)
				{
				  if (TwoDTranslationFlag == false)
				    {
				      PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				  else
				    {
				      cout << "Error: 2d translations not yet implemented" << endl;
				      return -1;
				    }
				}
			      else
				{
				  if (TwoDTranslationFlag == false)
				    {
				      PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				  else
				    {
				      PartialDensityMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				}
			    }
			  else
			    {
			      if (Manager.GetBoolean("decoupled") == false)
				{
				  if (TwoDTranslationFlag == false)
				    {
				      PartialDensityMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				  else
				    {
				      cout << "Error: 2d translations not yet implemented" << endl;
				      return -1;
				    }
				}
			      else
				{
				  if (TwoDTranslationFlag == false)
				    {
				      PartialDensityMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      for (int i = 1; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				  else
				    {
				      cout << "Error: 2d translations not yet implemented" << endl;
				      return -1;
				    }
				}
			    }
			}
		    }
		  else
		    {
		      cout << "Error: Bosonic statistics not implemented" << endl;
		      return -1;
		    }
		  ++TmpIndex;
		  
		  for (; TmpIndex < TotalNbrSites; ++TmpIndex)
		    {
		      if (NbrGroundStatePerMomentumSector[TmpIndex] != 0)
			{
			  if (Statistics == true)
			    {
			      if (SU2SpinFlag == false)
				{      
				  if (TwoDTranslationFlag == false)
				    {
				      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				  else
				    {				      
				      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
					{
					  HermitianMatrix TmpMatrix = ((FermionOnLatticeRealSpaceAnd2DTranslation*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
					  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
					  PartialDensityMatrix += TmpMatrix;
					}
				    }
				}
			      else
				{
				  if (GutzwillerFlag == false)
				    {
				      if (Manager.GetBoolean("decoupled") == false)
					{
					  if (TwoDTranslationFlag == false)
					    {
					      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						{
						  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						  PartialDensityMatrix += TmpMatrix;
						}
					    }
					  else
					    {
					      cout << "Error: 2d translations not yet implemented" << endl;
					      return -1;
					    }
					}
				      else
					{
					  if (TwoDTranslationFlag == false)
					    {
					      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						{
						  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						  PartialDensityMatrix += TmpMatrix;
						}
					    }
					  else
					    {
					      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						{
						  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						  PartialDensityMatrix += TmpMatrix;
						}
					    }
					}
				    }
				  else
				    {
				      if (Manager.GetBoolean("decoupled") == false)
					{
					  if (TwoDTranslationFlag == false)
					    {
					      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						{
						  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						  PartialDensityMatrix += TmpMatrix;
						}
					    }
					  else
					    {
					      cout << "Error: 2d translations not yet implemented" << endl;
					      return -1;
					    }
					}
				      else
					{
					  if (TwoDTranslationFlag == false)
					    {
					      for (int i = 0; i < NbrGroundStatePerMomentumSector[TmpIndex]; ++i)
						{
						  HermitianMatrix TmpMatrix = ((FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) Spaces [TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalSz, GroundStatePerMomentumSector[TmpIndex][i], Architecture.GetArchitecture());
						  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][i];
						  PartialDensityMatrix += TmpMatrix;
						}
					    }
					  else
					    {
					      cout << "Error: 2d translations not yet implemented" << endl;
					      return -1;
					    }
					}
				    }
				}
			    }
			  else
			    {
			      cout << "Error: Bosonic statistics not implemented" << endl;
			      return -1;
			    }			  
			}
		    }
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
		      
		    }
		  if (PartialDensityMatrix.GetNbrRow() > 1)
		    {
		      if (ShowTimeFlag == true)
			{
			  gettimeofday (&(TotalStartingTime), 0);
			  
			}
		      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
		      if (LapackFlag == true)
			{
			  PartialDensityMatrix.LapackDiagonalize(TmpDiag);
			}
		      else
			{
			  PartialDensityMatrix.Diagonalize(TmpDiag);
			}
#else
		      PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
		      TmpDiag.SortMatrixDownOrder();
		      
		      if (DensityMatrixFileName != 0)
			{ 
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 		      DensityMatrixFile.precision(14);
			  double Trace = 0.0;
			  if (TwoDTranslationFlag == false)
			    {
			      if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
				{
				  
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    {
				      DensityMatrixFile << SubsystemNbrParticles << " " << TmpDiag[i] << endl;
				      Trace += TmpDiag[i];
				    }
				  cout << "Trace = " << Trace << endl;
				}
			      else
				{
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    {
				      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << TmpDiag[i] << endl;
				      Trace += TmpDiag[i];
				    }
				  cout << "Trace = " << Trace << endl;
				}
			    }
			  else
			    {
			      if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
				{				  
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpDiag[i] << endl;
				}
			      else
				{
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalSz << " " << TmpDiag[i] << endl;
				}
			    }
			  DensityMatrixFile.close();
			}
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			{
			  if (TmpDiag[i] > 1e-14)
			    {
			      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
			      DensitySum +=TmpDiag[i];
			    }
			}
		      if (ShowTimeFlag == true)
			{
			  gettimeofday (&(TotalEndingTime), 0);
			  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
			  cout << "diagonalization done in " << Dt << "s" << endl;
			}
		    }
		  else
		    if (PartialDensityMatrix.GetNbrRow() == 1)
		      {
			double TmpValue = PartialDensityMatrix(0,0);
			if (DensityMatrixFileName != 0)
			  {
			    ofstream DensityMatrixFile;
			    DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 		DensityMatrixFile.precision(14);
			    if (TwoDTranslationFlag == false)
			      {
				if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
				  {
				    for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				      DensityMatrixFile << SubsystemNbrParticles << " " << TmpValue << endl;
				  }
				else
				  {
				    for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << TmpValue << endl;
				  }
			      }
			    else
			      {
				if ((SU2SpinFlag == false) || (Manager.GetBoolean("decoupled") == false))
				  {
				    for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpValue << endl;
				  }
				else
				  {
				    for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalSz 
							<< " " <<  TmpValue << endl;
				  }
			      }
			    DensityMatrixFile.close();
			  }		  
			if (TmpValue > 1e-14)
			  {
			    EntanglementEntropy += TmpValue * log(TmpValue);
			    DensitySum += TmpValue;
			  }
		      }
		}
	    }
	}
      File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
    }
  File.close();
  
  return 0;
}
