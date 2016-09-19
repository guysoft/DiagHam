#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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
  OptionManager Manager ("FQHESphereWithSU2SpinEntanglementEntropyParticlePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new BooleanOption ('\n', "no-sz", "indicates that the input states are not Sz eigenstates");
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*SystemGroup) += new SingleStringOption  ('\n', "realspace-cut", "use real space partition instead of particle partition, providing the orbital weights");
  (*SystemGroup) += new BooleanOption  ('\n', "use-alt", "use alternative Hilbert space for  bosonic states");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "lza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Lz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "sza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Sz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSU2SpinEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereWithSU2SpinEntanglementEntropyParticlePartition -h" << endl;
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


  int NbrParticles = 0; 
  int LzMax = 0; 
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  bool SVDFlag = Manager.GetBoolean("use-svd");
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterLza = Manager.GetInteger("lza-eigenstate");
  int FilterSza = Manager.GetInteger("sza-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  bool NoSzFlag = Manager.GetBoolean("no-sz");
  bool ComplexFlag = Manager.GetBoolean("complex");
  bool RealSpaceCut = false;
  if (Manager.GetString("realspace-cut") != 0)
    {
      RealSpaceCut = true;
    }
  int* TotalLz = 0;
  int* TotalSz = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  ParticleOnSphereWithSpin** Spaces = 0;
  RealVector* GroundStates = 0;
  ComplexVector* ComplexGroundStates = 0;
  char** GroundStateFiles = 0;

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalLz = new int[1];
      TotalSz = new int[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));      
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
       TotalLz = new int[NbrSpaces];
       TotalSz = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalLz[i] = 0;
      TotalSz[i] = 0;
      if (NoSzFlag == false)
	{
	  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(GroundStateFiles[i],
								   NbrParticles, LzMax, TotalLz[i], TotalSz[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
      else
	{
	  if (FQHEOnSphereFindSystemInfoFromVectorFileName(GroundStateFiles[i],
							   NbrParticles, LzMax, TotalLz[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
    }

  if (ComplexFlag == false)
    {
      GroundStates = new RealVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	  {
	    cout << "can't open vector file " << GroundStateFiles[i] << endl;
	    return -1;      
	  }
    }
  else
    {
      ComplexGroundStates = new ComplexVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	if (ComplexGroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	  {
	    cout << "can't open vector file " << GroundStateFiles[i] << endl;
	    return -1;      
	  }
    }
  
  Spaces = new ParticleOnSphereWithSpin* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (Statistics == true)
	{
	  if (NoSzFlag == false)
	    {
	      Spaces[i] = new FermionOnSphereWithSpin (NbrParticles, TotalLz[i], LzMax, TotalSz[i]);
	    }
	  else
	    {
	      cout << "error : bosons without Sz are not yet supported" << endl;
	      return 0;
//	      Spaces[i] = new FermionOnSphereWithSpin (NbrParticles, TotalLz[i], LzMax);
	    }
	  
	}
      else
	{
	  if (NoSzFlag == false)
	    {
	      if (Manager.GetBoolean("use-alt") == false)
		{
		  Spaces[i] = new BosonOnSphereWithSpin (NbrParticles, TotalLz[i], LzMax, TotalSz[i]);
		}
	      else
		{
		  Spaces[i] = new BosonOnSphereWithSU2Spin (NbrParticles, TotalLz[i], LzMax, TotalSz[i]);
		}
	    }
	  else
	    {
	      cout << "error : bosons without Sz are not yet supported" << endl;
	      return 0;
//	      Spaces[i] = new BosonOnSphereWithSpin (NbrParticles, TotalLz[i], LzMax);
	    }
	}
      if(ComplexFlag == false)
	{
	  if (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and ground state" << endl;
	      return 0;
	    }
	}
      else
	{
	  if (Spaces[i]->GetLargeHilbertSpaceDimension() != ComplexGroundStates[i].GetLargeVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and ground state" << endl;
	      return 0;
	    }
	}
    }
  
  int* NbrNUp = new int[NbrSpaces];
  int* NbrNDown = new int[NbrSpaces];
  int MaxNbrNUp = 0;
  int MaxNbrNDown = 0;
  for (int i = 0; i < NbrSpaces; ++i)
    { 
      NbrNUp[i] = (NbrParticles + TotalSz[i]);
      NbrNDown[i] = (NbrParticles - TotalSz[i]);
      NbrNUp[i] >>= 1;
      NbrNDown[i] >>= 1;
      if (NbrNUp[i] > MaxNbrNUp)
	MaxNbrNUp = NbrNUp[i];
      if (NbrNDown[i] > MaxNbrNDown)
	MaxNbrNDown = NbrNDown[i];
    }
  if (NoSzFlag == true)
    {
      MaxNbrNUp = 0;
      MaxNbrNDown = 0;
   }

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      if (NoSzFlag == false)
	{
	  DensityMatrixFile << "#  N    Sz    Nup    Ndown    Lz    lambda";
	}
      else
	{
	  DensityMatrixFile << "#  N    Lz    lambda";
	}
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  ofstream File;
  if (Manager.GetString("output-file") != 0)
    {
      File.open(Manager.GetString("output-file"), ios::binary | ios::out);
    }
  else
    {
      char* TmpFileName = 0;
      if (RealSpaceCut == false)
	{
	  TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
	}
      else
	{
	  TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "realent");	  
	}
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
    {
      MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
    }
  else
    {
      if (RealSpaceCut == true)
        MaxSubsystemNbrParticles = NbrParticles;
    }
  int SubsystemNbrParticles = Manager.GetInteger("min-na");
  if ((RealSpaceCut == true) && (SubsystemNbrParticles == 1))
    {
      SubsystemNbrParticles = 0;
    }
  double TotalTrace = 0;

  double* WeightAOrbitalsUp = 0;
  double* WeightBOrbitalsUp = 0;
  double* WeightAOrbitalsDown = 0;
  double* WeightBOrbitalsDown = 0;
  int NbrAOrbitals = LzMax + 1;
  int NbrBOrbitals = LzMax + 1;
  if (Manager.GetString("realspace-cut") != 0)
    {
      ConfigurationParser RealSpaceWeights;
      if (RealSpaceWeights.Parse(Manager.GetString("realspace-cut")) == false)
	{
	  RealSpaceWeights.DumpErrors(cout) << endl;
	  return -1;
	}
      double* TmpSquareWeights = 0;
      int TmpNbrOrbitals = 0;
      if (RealSpaceWeights.GetAsDoubleArray("OrbitalSquareWeights", ' ', TmpSquareWeights, TmpNbrOrbitals) == false)
	{
	  cout << "OrbitalSquareWeights is not defined or as a wrong value" << endl;
	  return -1;
	}
      if (TmpNbrOrbitals > (LzMax + 1))
	{
	  cout << "error, the number of weights (" << TmpNbrOrbitals << ") cannot exceed the number of orbitals (" << (LzMax + 1) << ")" << endl;
	  return -1;
	}
      NbrAOrbitals = (LzMax + 1 + TmpNbrOrbitals) / 2;
      WeightAOrbitalsUp = new double [NbrAOrbitals];
      WeightAOrbitalsDown = new double [NbrAOrbitals];
      for (int i = 0; i < (NbrAOrbitals - TmpNbrOrbitals); ++i)
	{
	  WeightAOrbitalsUp[i] = 1.0;
	  WeightAOrbitalsDown[i] = 1.0;
	}
      for (int i = NbrAOrbitals - TmpNbrOrbitals; i < NbrAOrbitals; ++i)
	{
	  WeightAOrbitalsUp[i] = sqrt(TmpSquareWeights[i - NbrAOrbitals + TmpNbrOrbitals]);
	  WeightAOrbitalsDown[i] = sqrt(TmpSquareWeights[i - NbrAOrbitals + TmpNbrOrbitals]);
	}
      NbrBOrbitals = (LzMax + 1 + TmpNbrOrbitals) / 2;
      WeightBOrbitalsUp = new double [NbrBOrbitals];
      WeightBOrbitalsDown = new double [NbrBOrbitals];
      for (int i = 0; i < TmpNbrOrbitals; ++i)
	{
	  WeightBOrbitalsUp[i] = sqrt(1.0 - TmpSquareWeights[i]);
	  WeightBOrbitalsDown[i] = sqrt(1.0 - TmpSquareWeights[i]);
	}
      for (int i = TmpNbrOrbitals; i < NbrBOrbitals; ++i)
	{
	  WeightBOrbitalsUp[i] = 1.0;
	  WeightBOrbitalsDown[i] = 1.0;
	}      
    }  
  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      for (int SubsystemNbrNUp = 0; SubsystemNbrNUp <= MaxNbrNUp; ++SubsystemNbrNUp)
	{
	  int SubsystemNbrNDown = SubsystemNbrParticles - SubsystemNbrNUp;
	  if (((SubsystemNbrNDown >= 0) && (SubsystemNbrNDown <= MaxNbrNDown)) || (NoSzFlag == true))
	    {
	      int SubsystemTotalSz = SubsystemNbrNUp - SubsystemNbrNDown;
	      
	      int SubsystemMaxTotalLz = 0;
	      if (Statistics == true)
		{
		  SubsystemMaxTotalLz = (((SubsystemNbrNUp * LzMax) - (SubsystemNbrNUp * (SubsystemNbrNUp - 1)))
					 + ((SubsystemNbrNDown * LzMax) - (SubsystemNbrNDown * (SubsystemNbrNDown - 1))));
		}
	      else
		{
		  SubsystemMaxTotalLz = SubsystemNbrParticles * LzMax;
		}
	      
	      int SubsystemTotalLz = -SubsystemMaxTotalLz; 
	      for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
		{
		  timeval TotalStartingTime;
		  timeval TotalEndingTime;
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
		  if (NoSzFlag == false)
		    {
		      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Nup=" << SubsystemNbrNUp 
			   << " subsystem total Ndown=" << SubsystemNbrNDown << " subsystem total Lz=" << SubsystemTotalLz << endl;
		    }
		  else
		    {
		      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
		    }
		  RealSymmetricMatrix PartialDensityMatrix;
		  HermitianMatrix ComplexPartialDensityMatrix;
		  RealMatrix PartialEntanglementMatrix; 
		  ComplexMatrix ComplexPartialEntanglementMatrix; 
		  
		  for (int i = 0; i < NbrSpaces; ++i)
		    {
		      int ComplementaryNbrNUp = NbrNUp[i] - SubsystemNbrNUp;
		      int ComplementaryNbrNDown = NbrNDown[i] - SubsystemNbrNDown;
		      int ComplementaryMaxTotalLz = 0;
		      if (Statistics == true)
			{
			  ComplementaryMaxTotalLz = (((ComplementaryNbrNUp * LzMax) - (ComplementaryNbrNUp * (ComplementaryNbrNUp - 1)))
						     + ((ComplementaryNbrNDown * LzMax) - ( ComplementaryNbrNDown * ( ComplementaryNbrNDown - 1))));
			}
		      else
			{
			  ComplementaryMaxTotalLz = (ComplementaryNbrNUp + ComplementaryNbrNDown) * LzMax;
			}
		      
		      if ((SubsystemNbrNUp <= NbrNUp[i]) && (SubsystemNbrNDown <= NbrNDown[i]) && (abs(TotalLz[i] - SubsystemTotalLz) <= ComplementaryMaxTotalLz ))
			{
			  RealSymmetricMatrix TmpMatrix;
			  HermitianMatrix ComplexTmpMatrix;
			  
			  if (Statistics == true)
			    {
			      if (ComplexFlag == false)
				{
				  if (SVDFlag == true)
				    {
				      if (RealSpaceCut == true)
					{
					  PartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz,  SubsystemNbrNUp - SubsystemNbrNDown, GroundStates[i] , true);
					  if(PartialEntanglementMatrix.GetNbrRow() != 0)
					    {
					      Spaces[i]->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, SubsystemNbrNUp - SubsystemNbrNDown ,Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), PartialEntanglementMatrix);
					    }
					}
				      else
					{
					  PartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz , 
																    SubsystemNbrNUp - SubsystemNbrNDown, GroundStates[i] , false);
					}
				    }
				  else
				    {
				      TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, 
													   SubsystemNbrNUp - SubsystemNbrNDown, GroundStates[i]);
				    }
				}
			    }
			  else
			    {
			      if (ComplexFlag == false)
				{
				  if (SVDFlag == true)
				    {
				      if (RealSpaceCut == true)
					{
					  PartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz,  SubsystemNbrNUp - SubsystemNbrNDown, GroundStates[i] , true);
					  if(PartialEntanglementMatrix.GetNbrRow() != 0)
					    {
					      Spaces[i]->EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, SubsystemNbrNUp - SubsystemNbrNDown,
																	   NbrAOrbitals, WeightAOrbitalsUp, WeightAOrbitalsDown,
																	   NbrBOrbitals, WeightBOrbitalsUp, WeightBOrbitalsDown,
																	   PartialEntanglementMatrix);
					    }
					}
				      else
					{
					  PartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz , 
																    SubsystemNbrNUp - SubsystemNbrNDown, GroundStates[i] , false);
					}
				    }
				  else
				    {
				      if (NoSzFlag == false)
					{
					  if (ComplexFlag == false)
					    TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, 
														 SubsystemNbrNUp, SubsystemNbrNDown, GroundStates[i], Architecture.GetArchitecture());
					  //else
					  //ComplexTmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz,   SubsystemNbrNUp, SubsystemNbrNDown, ComplexGroundStates[i], Architecture.GetArchitecture());
					}
				      else
					{
					  if (ComplexFlag == false)
					    TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, 
														 GroundStates[i], Architecture.GetArchitecture());
					  else
					    ComplexTmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz,  ComplexGroundStates[i], Architecture.GetArchitecture());
					  
					}
				    }
				}
			    }
			  if (SVDFlag == false)
			    {
			      if (ComplexFlag == false)
				{
				  if (PartialDensityMatrix.GetNbrRow() > 0)
				    PartialDensityMatrix += TmpMatrix;
				  else
				    PartialDensityMatrix = TmpMatrix;
				}
			      else
				{
				  if (ComplexPartialDensityMatrix.GetNbrRow() > 0)
				    ComplexPartialDensityMatrix += ComplexTmpMatrix;
				  else
				    ComplexPartialDensityMatrix = ComplexTmpMatrix;
				}
			    }
			}
		    }
		  if (ComplexFlag == false)
		    {
		      if (NbrSpaces > 1)
			{ 
			  if (SVDFlag == true)
			    {
			      cout <<"SVD is not possible for more than one state"<<endl;
			      return -1;
			    }
			  else
			    PartialDensityMatrix /= ((double) NbrSpaces);
			}
		    }
		  else
		    {
		      if (NbrSpaces > 1)
			ComplexPartialDensityMatrix /= ((double) NbrSpaces);
		    }
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
		    }
		  if (ComplexFlag == false)
		    {
		      if ((PartialDensityMatrix.GetNbrRow() > 1)||(PartialEntanglementMatrix.GetNbrRow() >= 1))
			{
			  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
			  if (SVDFlag == false)
			    {
#ifdef __LAPACK__
			      if (LapackFlag == true)
				{
				  if ((EigenstateFlag == true) && (FilterLza == SubsystemTotalLz) && (FilterSza == SubsystemTotalSz))
				    {
				      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
								PartialDensityMatrix.GetNbrRow(), true);
				      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
					TmpEigenstates[i][i] = 1.0;
				      PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
				      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				      char* TmpEigenstateName = new char[512];
				      int MaxNbrEigenstates = NbrEigenstates;
				      if (NbrEigenstates == 0)
					MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				      for (int i = 0; i < MaxNbrEigenstates; ++i)
					{
					  if (TmpDiag[i] > 1e-14)
					    {
					      if (NoSzFlag == false)
						{
						  if (Statistics == true)
						    {
						      sprintf (TmpEigenstateName,
							       "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
							       NbrParticles, LzMax, TotalLz[0], 
							       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
						    }
						  else
						    {
						      sprintf (TmpEigenstateName,
							       "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
							       NbrParticles, LzMax, TotalLz[0], 
							       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
						    }
						}
					      else
						{
						  if (Statistics == true)
						    {
						      sprintf (TmpEigenstateName,
							       "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
							       NbrParticles, LzMax, TotalLz[0], 
							       SubsystemNbrParticles, SubsystemTotalLz, i);
						    }
						  else
						    {
						      sprintf (TmpEigenstateName,
							       "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
							       NbrParticles, LzMax, TotalLz[0], 
							       SubsystemNbrParticles, SubsystemTotalLz, i);
						    }
						}
					      TmpEigenstates[i].WriteVector(TmpEigenstateName);
					    }
					}
				      delete[] TmpEigenstateName;
				    }
				  else
				    {
				      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
				    }
				}
			      else
				{
				  if ((EigenstateFlag == true) && (FilterLza == SubsystemTotalLz) && (FilterSza == SubsystemTotalSz))
				    {
				      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
								PartialDensityMatrix.GetNbrRow(), true);
				      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
					TmpEigenstates[i][i] = 1.0;
				      PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
				      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				      char* TmpEigenstateName = new char[512];
				      int MaxNbrEigenstates = NbrEigenstates;
				      if (NbrEigenstates == 0)
					MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				      for (int i = 0; i < MaxNbrEigenstates; ++i)
					{
					  if (TmpDiag[i] > 1e-14)
					    {
					      if (NoSzFlag == false)
						{
						  if (Statistics == true)
						    {
						      sprintf (TmpEigenstateName,
							       "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
							       NbrParticles, LzMax, TotalLz[0], 
							       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
						    }
						  else
						    {
						      sprintf (TmpEigenstateName,
							       "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
							       NbrParticles, LzMax, TotalLz[0], 
							       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
						    }
						}
					      else
						{
						  if (Statistics == true)
						    {
						      sprintf (TmpEigenstateName,
							       "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
							       NbrParticles, LzMax, TotalLz[0], 
							       SubsystemNbrParticles, SubsystemTotalLz, i);
						    }
						  else
						    {
						      sprintf (TmpEigenstateName,
							       "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
							       NbrParticles, LzMax, TotalLz[0], 
							       SubsystemNbrParticles, SubsystemTotalLz, i);
						    }
						}
					      TmpEigenstates[i].WriteVector(TmpEigenstateName);
					    }
					}
				      delete[] TmpEigenstateName;
				    }
				  else
				    {
				      PartialDensityMatrix.Diagonalize(TmpDiag);
				    }
				}
#else
			      PartialDensityMatrix.Diagonalize(TmpDiag);
#endif 
			    }
			  else
			    {
			      double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
			      int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
			      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
				{
				  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
				}
			      for (int i = 0; i < TmpDimension; ++i)
				{
				  TmpValues[i] *= TmpValues[i];
				}
			      TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);  
			    }
			  
			  TmpDiag.SortMatrixDownOrder();
			  if (DensityMatrixFileName != 0)
			    {
			      ofstream DensityMatrixFile;
			      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			      DensityMatrixFile.precision(14);
			      if (NoSzFlag == false)
				{
				  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
				    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << SubsystemNbrNUp << " " 
						      <<  SubsystemNbrNDown << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
				}
			      else
				{
				  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
				    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
				}
			      DensityMatrixFile.close();
			    }
			  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
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
				if (NoSzFlag == false)
				  {
				    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << SubsystemNbrNUp << " " 
						      << SubsystemNbrNDown << " " << SubsystemTotalLz << " " << TmpValue << endl;
				  }
				else
				  {
				    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
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
		  else
		    {
		      if (ComplexPartialDensityMatrix.GetNbrRow() > 1)
			{
			  RealDiagonalMatrix TmpDiag (ComplexPartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
			  if (LapackFlag == true)
			    {
			      if ((EigenstateFlag == true) && (FilterLza == SubsystemTotalLz) && (FilterSza == SubsystemTotalSz))
				{
				  ComplexMatrix TmpEigenstates(ComplexPartialDensityMatrix.GetNbrRow(),
							       ComplexPartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < ComplexPartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  ComplexPartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[512];
				  int MaxNbrEigenstates = NbrEigenstates;
			      if (NbrEigenstates == 0)
				MaxNbrEigenstates = ComplexPartialDensityMatrix.GetNbrRow();
			      for (int i = 0; i < MaxNbrEigenstates; ++i)
				{
				  if (TmpDiag[i] > 1e-14)
				    {
				      if (NoSzFlag == false)
					{
					  if (Statistics == true)
					    {
					      sprintf (TmpEigenstateName,
						       "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
						       NbrParticles, LzMax, TotalLz[0], 
						       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
					    }
					  else
					    {
					      sprintf (TmpEigenstateName,
						       "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
						       NbrParticles, LzMax, TotalLz[0], 
						       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
					    }
					}
				      else
					{
					  if (Statistics == true)
					    {
					      sprintf (TmpEigenstateName,
						       "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
						       NbrParticles, LzMax, TotalLz[0], 
						       SubsystemNbrParticles, SubsystemTotalLz, i);
					    }
					  else
					    {
					      sprintf (TmpEigenstateName,
						       "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
						       NbrParticles, LzMax, TotalLz[0], 
						       SubsystemNbrParticles, SubsystemTotalLz, i);
					    }
					}
				      TmpEigenstates[i].WriteVector(TmpEigenstateName);
				    }
				}
			      delete[] TmpEigenstateName;
			    }
			  else
			    {
			      ComplexPartialDensityMatrix.LapackDiagonalize(TmpDiag);
			    }
			}
		      else
			{
			  if ((EigenstateFlag == true) && (FilterLza == SubsystemTotalLz) && (FilterSza == SubsystemTotalSz))
			    {
			      ComplexMatrix TmpEigenstates(ComplexPartialDensityMatrix.GetNbrRow(),
							   ComplexPartialDensityMatrix.GetNbrRow(), true);
			      for (int i = 0; i < ComplexPartialDensityMatrix.GetNbrRow(); ++i)
				TmpEigenstates[i][i] = 1.0;
			      ComplexPartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
			      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
			      char* TmpEigenstateName = new char[512];
			      int MaxNbrEigenstates = NbrEigenstates;
			      if (NbrEigenstates == 0)
				MaxNbrEigenstates = ComplexPartialDensityMatrix.GetNbrRow();
			      for (int i = 0; i < MaxNbrEigenstates; ++i)
				{
				  if (TmpDiag[i] > 1e-14)
				    {
				      if (NoSzFlag == false)
					{
					  if (Statistics == true)
					    {
					      sprintf (TmpEigenstateName,
						       "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
						       NbrParticles, LzMax, TotalLz[0], 
						       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
					    }
					  else
					    {
					      sprintf (TmpEigenstateName,
						       "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
						       NbrParticles, LzMax, TotalLz[0], 
						       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
					    }
					}
				      else
					{
					  if (Statistics == true)
					    {
					      sprintf (TmpEigenstateName,
						       "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
						       NbrParticles, LzMax, TotalLz[0], 
						       SubsystemNbrParticles, SubsystemTotalLz, i);
					    }
					  else
					    {
					      sprintf (TmpEigenstateName,
						       "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
						       NbrParticles, LzMax, TotalLz[0], 
						       SubsystemNbrParticles, SubsystemTotalLz, i);
					    }
					}
				      TmpEigenstates[i].WriteVector(TmpEigenstateName);
				    }
				}
			      delete[] TmpEigenstateName;
			    }
			  else
			    {
			      ComplexPartialDensityMatrix.Diagonalize(TmpDiag);
			    }
			}
#else
		      ComplexPartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
		      TmpDiag.SortMatrixDownOrder();
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  if (NoSzFlag == false)
			    {
			      for (int i = 0; i < ComplexPartialDensityMatrix.GetNbrRow(); ++i)
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << SubsystemNbrNUp << " " 
						  <<  SubsystemNbrNDown << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			    }
			  else
			    {
			      for (int i = 0; i < ComplexPartialDensityMatrix.GetNbrRow(); ++i)
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			    }
			  DensityMatrixFile.close();
			}
		      for (int i = 0; i < ComplexPartialDensityMatrix.GetNbrRow(); ++i)
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
		    if (ComplexPartialDensityMatrix.GetNbrRow() == 1)
		      {
			double TmpValue = ComplexPartialDensityMatrix(0,0);
			if (DensityMatrixFileName != 0)
			  {
			    ofstream DensityMatrixFile;
			    DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			    DensityMatrixFile.precision(14);
			    if (NoSzFlag == false)
			      {
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << SubsystemNbrNUp << " " 
						  << SubsystemNbrNDown << " " << SubsystemTotalLz << " " << TmpValue << endl;
			      }
			    else
			      {
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
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
	    }
      if (RealSpaceCut == true)
	TotalTrace += DensitySum;

	  File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
	}
  if (RealSpaceCut == true)
    {
      cout <<" Total density matrix trace is equal to "<< TotalTrace<<endl;
      File << "# Total density matrix trace is equal to "<< TotalTrace<<endl;
    }
  File.close();
}

