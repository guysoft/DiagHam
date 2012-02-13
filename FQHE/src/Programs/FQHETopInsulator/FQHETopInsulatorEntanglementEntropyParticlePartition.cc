#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSpinMomentumSpace.h"

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
  OptionManager Manager ("FQHETopInsulatorEntanglementEntropyParticlePartition" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "3d", "consider a 3d model instead of a 2d model");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETopInsulatorEntanglementEntropyParticlePartition -h" << endl;
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
  int* TotalKz = 0;
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  int NbrSiteZ = 0;
  bool Statistics = true;
  double* Coefficients = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  bool Flag3d = Manager.GetBoolean("3d");

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHETopInsulatorEntanglementEntropyParticlePartition -h" << endl;
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
      TotalKz = new int[1];
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
       TotalKz = new int[NbrSpaces];
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

  if (Flag3d == false)
    {
      NbrSiteZ = 1;
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  TotalKx[i] = 0;
	  TotalKy[i] = 0;
	  TotalKz[i] = 0;
	  double Mass = 0.0;
	  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i],
								  NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], Mass, Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
    }
  else
    {
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  TotalKx[i] = 0;
	  TotalKy[i] = 0;
	  TotalKz[i] = 0;
	  if (FQHEOnCubicLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i],
								 NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	  cout << GroundStateFiles[i] << " " << NbrParticles << " " << NbrSiteX << " " << NbrSiteY << " " << NbrSiteZ << " " << TotalKx[i] << " " << TotalKy[i] << " " << TotalKz[i] << endl;
	}
    }


  GroundStates = new ComplexVector [NbrSpaces];  
  int TotalNbrSites = NbrSiteX * NbrSiteY * NbrSiteZ;
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
      int TmpIndex = (((TotalKx[i] * NbrSiteY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
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
      int TmpIndex = (((TotalKx[i] * NbrSiteY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
      GroundStatePerMomentumSector[TmpIndex][NbrGroundStatePerMomentumSector[TmpIndex]] = GroundStates[i];
      CoefficientPerMomentumSector[TmpIndex][NbrGroundStatePerMomentumSector[TmpIndex]] = Coefficients[i];
      NbrGroundStatePerMomentumSector[TmpIndex]++;
    }  

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      if (Flag3d == false)
	{
	  DensityMatrixFile << "#  N    Kx    Ky    lambda";
	}
      else
	{
	  DensityMatrixFile << "#  N    Kx    Ky    Kz    lambda";
	}
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  int MaxNbrSpaces = NbrSiteX * NbrSiteY * NbrSiteZ;
  ParticleOnSphere** Spaces = new ParticleOnSphere*[MaxNbrSpaces];
  for (int i = 0; i < MaxNbrSpaces; ++i)
    {
      Spaces[i] = 0;
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpIndex = (((TotalKx[i] * NbrSiteY) + TotalKy[i]) * NbrSiteZ) + TotalKz[i];
      if (Spaces[TmpIndex] == 0)
	{
	  if (Flag3d == false)
	    {
	      if (Statistics == true)
		Spaces[TmpIndex] = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
	      else
		Spaces[TmpIndex] = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
	    }
	  else
	    {
	      Spaces[TmpIndex] = new FermionOnCubicLatticeWithSpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx[i], TotalKy[i], TotalKz[i]);
	    }
	}
      if (Spaces[TmpIndex]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
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
      for (int SubsystemTotalKx = 0; SubsystemTotalKx < NbrSiteX; ++SubsystemTotalKx)
	{
	  for (int SubsystemTotalKy = 0; SubsystemTotalKy < NbrSiteY; ++SubsystemTotalKy)
	    {
	      for (int SubsystemTotalKz = 0; SubsystemTotalKz < NbrSiteZ; ++SubsystemTotalKz)
		{
		  if (Flag3d == false)
		    cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy << endl;
		  else
		    cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy << " Kz=" << SubsystemTotalKz << endl;
		  
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
		  if (Flag3d == false)
		    {
		      if (Statistics == true)
			{
			  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
			    {
			      PartialDensityMatrix = ((FermionOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			    }
			  else
			    {
			      PartialDensityMatrix = ((FermionOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
			    }
			}
		      else
			{
			  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
			    {
			      PartialDensityMatrix = ((BosonOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			      PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			    }
			  else
			    {
 			      PartialDensityMatrix = ((BosonOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
			    }
			}
		    }
		  else
		    {
		      if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
			{
			  PartialDensityMatrix = ((FermionOnCubicLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
			  PartialDensityMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
			}
		      else
			{
			  PartialDensityMatrix = ((FermionOnCubicLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
			}
		    }
		  
		  ++TmpIndex;
		  
		  for (; TmpIndex < TotalNbrSites; ++TmpIndex)
		    {
		      if (NbrGroundStatePerMomentumSector[TmpIndex] != 0)
			{
			  if (Flag3d == false)
			    {
			      if (Statistics == true)
				{
				  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
				    {
				      HermitianMatrix TmpMatrix = ((FermionOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      PartialDensityMatrix += TmpMatrix;
				    }
				  else
				    {
				      HermitianMatrix TmpMatrix = ((FermionOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
				      PartialDensityMatrix += TmpMatrix;
				    }
				}
			      else
				{
				  if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
				    {
				      HermitianMatrix TmpMatrix = ((BosonOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				      TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				      PartialDensityMatrix += TmpMatrix;
				    }
				  else
				    {
// 				      HermitianMatrix TmpMatrix = ((BosonOnSquareLatticeMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
// 				      PartialDensityMatrix += TmpMatrix;
				    }
				}
			    }
			  else
			    {
			      if (NbrGroundStatePerMomentumSector[TmpIndex] == 1)
				{
				  HermitianMatrix TmpMatrix = ((FermionOnCubicLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, GroundStatePerMomentumSector[TmpIndex][0], Architecture.GetArchitecture());
				  TmpMatrix *= CoefficientPerMomentumSector[TmpIndex][0];
				  PartialDensityMatrix += TmpMatrix;
				}
			      else
				{
				  HermitianMatrix TmpMatrix = ((FermionOnCubicLatticeWithSpinMomentumSpace*) Spaces[TmpIndex])->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKx, SubsystemTotalKy, SubsystemTotalKz, NbrGroundStatePerMomentumSector[TmpIndex], GroundStatePerMomentumSector[TmpIndex], CoefficientPerMomentumSector[TmpIndex], Architecture.GetArchitecture());
				  PartialDensityMatrix += TmpMatrix;
				}
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
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  if (Flag3d == false)
			    {
			      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpDiag[i] << endl;
			    }
			  else
			    {
			      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalKz << " " << TmpDiag[i] << endl;
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
			  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
						((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
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
			    DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			    DensityMatrixFile.precision(14);
			    if (Flag3d == false)
			      {
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpValue << endl;
			      }
			    else
			      {
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << SubsystemTotalKz << " " << TmpValue << endl;
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
