#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/BosonOnSphereLong.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisLong.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// evaluate Hilbert space dimension for bosons
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension
long BosonEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for bosons
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension
long BosonShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereEntanglementEntropyParticlePartition" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new BooleanOption  ('\n', "compute-lvalue", "compute the L value of each reduced density matrix eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n', "largest-lz", "only compute the largest block of the reduced density matrix (Lz=0 or 1/2)");
  (*SystemGroup) += new BooleanOption  ('\n', "positive-lz", "only compute the positive Lz sectors");
  (*SystemGroup) += new BooleanOption  ('\n', "realspace-cut", "use real space partition instead of particle partition");
  (*SystemGroup) += new SingleDoubleOption ('\n', "realspace-theta-top", "inclination angle that defines the top of the real space parition (in degrees)", 0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "realspace-theta-bot", "inclination angle that defines the bottom of the real space parition (in degrees)", 90);
  (*SystemGroup) += new SingleDoubleOption ('\n', "realspace-phi-range", "angle between the 2 longitudes that defines the real space parition (in degrees)", 360);
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*SystemGroup) += new BooleanOption  ('\n', "entanglement-matrix", "compute the entanglement matrix instead of the density matrix");
  (*SystemGroup) += new BooleanOption  ('\n', "check-size", "check the amount of storage required for the entanglement matrix calculation");  
  (*SystemGroup) += new BooleanOption  ('\n', "matrix-realcut", "compute the entanglement matrix for the real space partition");
  (*SystemGroup) += new BooleanOption  ('\n', "use-entanglement", "use previously computed entanglement matrix to compute the density matrix");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new BooleanOption ('\n', "lza-filter", "compute the eigenstates of the reduced density matrix only for a given value of Lza");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "lza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Lz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereEntanglementEntropyParticlePartition -h" << endl;
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


  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int LzMax = Manager.GetInteger("lzmax"); 
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  bool ComputeLValueFlag = Manager.GetBoolean("compute-lvalue");
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  bool LargestLSector = Manager.GetBoolean("largest-lz");
  bool PositiveLzSectors = Manager.GetBoolean("positive-lz");
  bool RealSpaceCut = Manager.GetBoolean("realspace-cut");
  int FilterLza = Manager.GetInteger("lza-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  bool UseEntanglementMatrixFlag = Manager.GetBoolean("use-entanglement");
  bool EntanglementMatrixForRealCut = Manager.GetBoolean("matrix-realcut");
  int* TotalLz = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  ParticleOnSphere** Spaces = 0;
  RealVector* GroundStates = 0;
  char** GroundStateFiles = 0;

  if ((ComputeLValueFlag == true) && (DensityMatrixFileName == 0))
    {
      cout << "compute-lvalue only valid when density-matrix is activated" << endl;
      return - 1;
    }

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalLz = new int[1];
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
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalLz[i] = 0;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(GroundStateFiles[i],
						       NbrParticles, LzMax, TotalLz[i], Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
      if (Statistics == true)
	{
	  cout << GroundStateFiles[i] << " is not a bosonic state" << endl;
	  return -1;
	}
      if (((NbrParticles * LzMax) & 1) != (TotalLz[i] & 1))
	{
	  cout << "incompatible values for nbr-particles, nbr-flux and total-lz for ground state file " << GroundStateFiles[i] << endl;
	  return -1;
	}
    }


  GroundStates = new RealVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[i] << endl;
	return -1;      
      }


  Spaces = new ParticleOnSphere* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
#ifdef  __64_BITS__
      if ((LzMax + NbrParticles - 1) < 63)
#else
	if ((LzMax + NbrParticles - 1) < 31)	
#endif
	  {
	    if (Manager.GetBoolean("haldane") == false)
	      {
		Spaces[i] = new BosonOnSphereShort (NbrParticles, TotalLz[i], LzMax);
	      }
	    else
	      {
		int* ReferenceState = 0;
		if (Manager.GetString("reference-file") == 0)
		  {
		    cout << "error, a reference file is needed" << endl;
		    return 0;
		  }
		ConfigurationParser ReferenceStateDefinition;
		if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
		  {
		    ReferenceStateDefinition.DumpErrors(cout) << endl;
		    return 0;
		  }
		if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
		  {
		    cout << "NbrParticles is not defined or as a wrong value" << endl;
		    return 0;
		  }
		if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax < 0))
		  {
		    cout << "LzMax is not defined or as a wrong value" << endl;
		    return 0;
		  }
		int MaxNbrLz;
		if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		  {
		    cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
		    return 0;     
		  }
		if (MaxNbrLz != (LzMax + 1))
		  {
		    cout << "wrong LzMax value in ReferenceState" << endl;
		    return 0;     
		  }
		if (Manager.GetString("load-hilbert") != 0)
		  Spaces[i] = new BosonOnSphereHaldaneBasisShort(Manager.GetString("load-hilbert"));
		else
		  {
		    Spaces[i] = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz[i], LzMax, ReferenceState);	  
		  }
	      }
	  }
	else
	  {
	    Spaces[i] = new BosonOnSphereLong (NbrParticles, TotalLz[i], LzMax);
	  }
      
      if (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }


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
  if (((RealSpaceCut == true) || (EntanglementMatrixForRealCut == true)) && (SubsystemNbrParticles == 1))
    {
      SubsystemNbrParticles = 0;
    }

  if (Manager.GetBoolean("entanglement-matrix") == false)
    {
      if (DensityMatrixFileName != 0)
	{
	  ofstream DensityMatrixFile;
	  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
	  DensityMatrixFile << "#  N    Lz    lambda";
	  if (ComputeLValueFlag == true)
	    DensityMatrixFile << " L^2 L";
	  DensityMatrixFile << endl;
	  DensityMatrixFile.close();
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
      
      for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  double EntanglementEntropy = 0.0;
	  double DensitySum = 0.0;
	  
	  int SubsystemMaxTotalLz = SubsystemNbrParticles * LzMax;
	  int SubsystemTotalLz = -SubsystemMaxTotalLz; 
	  if (LargestLSector == true)
	    {
	      if (((LzMax * NbrParticles) & 1) == 0)
		{
		  SubsystemTotalLz = 0;
		  SubsystemMaxTotalLz = 0;
		}
	      else
		{
		  SubsystemTotalLz = 1;
		  SubsystemMaxTotalLz = 1;
		}
	    }
	  if (PositiveLzSectors == true)
	    {
	      if (((LzMax * NbrParticles) & 1) == 0)
		{
		  SubsystemTotalLz = 0;
		}
	      else
		{
		  SubsystemTotalLz = 1;
		}
	    }
	  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
	    {
	      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
	      timeval TotalStartingTime;
	      timeval TotalEndingTime;
	      if (ShowTimeFlag == true)
		{
		  gettimeofday (&(TotalStartingTime), 0);
		}
	      RealSymmetricMatrix PartialDensityMatrix;
	      if (UseEntanglementMatrixFlag == false)
		{
		  if (RealSpaceCut == false)	    
		    PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[0], Architecture.GetArchitecture());
		  else
		    {
		      if ((2 * SubsystemNbrParticles) <= NbrParticles)
			{
			  PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz,Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), GroundStates[0], Architecture.GetArchitecture());
			}
		      else
			{
			  if (Manager.GetDouble("realspace-theta-top") != 0.0)
			    {
			      PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz,Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), GroundStates[0], Architecture.GetArchitecture());
			    }
			  else
			    {
			      PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixRealSpacePartition(NbrParticles - SubsystemNbrParticles, TotalLz[0] - SubsystemTotalLz, Manager.GetDouble("realspace-theta-bot"), 180.0, Manager.GetDouble("realspace-phi-range"), GroundStates[0], Architecture.GetArchitecture());		    
			    }
			}
		    }
		  for (int i = 1; i < NbrSpaces; ++i)
		    {
		      RealSymmetricMatrix TmpMatrix;
		      if (RealSpaceCut == false)
			{
			  TmpMatrix =  Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[i], Architecture.GetArchitecture());
			}
		      else
			{
			  if ((2 * SubsystemNbrParticles) <= NbrParticles)
			    {
			      TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz, Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), GroundStates[i], Architecture.GetArchitecture());
			    }
			  else
			    {
			      if (Manager.GetDouble("realspace-theta-top") != 0.0)
				{
				  TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz, Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), GroundStates[i], Architecture.GetArchitecture());
				}
			      else
				{
				  TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixRealSpacePartition(NbrParticles - SubsystemNbrParticles, TotalLz[i] - SubsystemTotalLz, Manager.GetDouble("realspace-theta-bot"), 180.0, Manager.GetDouble("realspace-phi-range"), GroundStates[i], Architecture.GetArchitecture());
				}
			    }
			}
		      PartialDensityMatrix += TmpMatrix;
		    }
		  if (NbrSpaces > 1)
		    PartialDensityMatrix /= ((double) NbrSpaces);
		}
	      else
		{
		  char* TmpExtension = new char[256];
		  if (RealSpaceCut == true)
		    {
		      if ((2 * SubsystemNbrParticles) <= NbrParticles)
			sprintf (TmpExtension, "pem_realcut_na_%d_lza_%d.mat", SubsystemNbrParticles, SubsystemTotalLz);
		      else
			sprintf (TmpExtension, "pem_realcut_na_%d_lza_%d.mat", NbrParticles - SubsystemNbrParticles, TotalLz[0] - SubsystemTotalLz);
		    }
		  else
		    {
		      sprintf (TmpExtension, "pem_na_%d_lza_%d.mat", SubsystemNbrParticles, SubsystemTotalLz);
		}
		  char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", TmpExtension);
		  if (TmpFileName == 0)
		    {
		      cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
		      return 0;
		    }
		  if (IsFile(TmpFileName) == true)
		    {
		      RealMatrix TmpEntanglementMatrix;
		      if (TmpEntanglementMatrix.ReadMatrix(TmpFileName) == false)
			{
			  cout << "can't read " <<  TmpFileName << endl;
			  return -1;
			}
		      if (RealSpaceCut == false)	    
			{
			  PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixParticlePartitionFromEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, TmpEntanglementMatrix);
			}
		      else
			{
			  if ((2 * SubsystemNbrParticles) <= NbrParticles)
			    {
			      if (SubsystemNbrParticles != 0)
				PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixRealSpacePartitionFromEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz,Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), TmpEntanglementMatrix);
			      else
				PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz,Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), GroundStates[0]);
			    }
			  else
			    {
			      if (SubsystemNbrParticles != NbrParticles)
				PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixRealSpacePartitionFromEntanglementMatrix(NbrParticles - SubsystemNbrParticles, TotalLz[0] - SubsystemTotalLz, Manager.GetDouble("realspace-theta-bot"), 180.0, Manager.GetDouble("realspace-phi-range"), TmpEntanglementMatrix);
			      else
				PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixRealSpacePartition(NbrParticles - SubsystemNbrParticles, TotalLz[0] - SubsystemTotalLz, Manager.GetDouble("realspace-theta-bot"), 180.0, Manager.GetDouble("realspace-phi-range"), GroundStates[0]);
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
		  if (ComputeLValueFlag == false)
		    {
#ifdef __LAPACK__
		      if (LapackFlag == true)
			PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		      else
			PartialDensityMatrix.Diagonalize(TmpDiag);
#else
		      PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
		      TmpDiag.SortMatrixDownOrder();
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			  DensityMatrixFile.close();
			}
		    }
		  else
		    {
		      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
						PartialDensityMatrix.GetNbrRow(), true);
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			TmpEigenstates[i][i] = 1.0;
#ifdef __LAPACK__
		      if (LapackFlag == true)
			PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
		      else
			PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates);
#else
		      PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates);
#endif
		      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
		      BosonOnSphereShort TmpDestinationHilbertSpace(SubsystemNbrParticles, SubsystemTotalLz, LzMax);
		      ParticleOnSphereSquareTotalMomentumOperator OperMomentum (&TmpDestinationHilbertSpace, LzMax);
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      char* TmpEigenstateName = new char[512];
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			{
			  double TmpSqrMomentum = (OperMomentum.MatrixElement(TmpEigenstates[i], TmpEigenstates[i])).Re;
			  double TmpMomentum = 0.0;
			  if (TmpSqrMomentum > 0.0)
			    TmpMomentum = (0.5 * (sqrt ((4.0 * TmpSqrMomentum) + 1.0) - 1.0));
			  DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << " " << TmpSqrMomentum << " " << TmpMomentum << endl;
			  if ((EigenstateFlag == true) && ((Manager.GetBoolean("lza-filter") == false) || (FilterLza == SubsystemTotalLz)) && ((NbrEigenstates == 0) || (NbrEigenstates > i)))
			    {
			      sprintf (TmpEigenstateName,
				       "bosons_particlereduceddensity_na_%d_n_%d_2s_%d_lz_%d.%d.vec",
				       SubsystemNbrParticles, NbrParticles, LzMax, SubsystemTotalLz, i);
			      TmpEigenstates[i].WriteVector(TmpEigenstateName);
			    }
			}
		      delete[] TmpEigenstateName;
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
			if (ComputeLValueFlag == false)
			  {
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
			  }
			else		      
			  {
			    if (SubsystemNbrParticles == 1)
			      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << " " << ((LzMax * (LzMax + 2)) / 4.0) << " " << (LzMax / 2.0) << endl;
			    else
			      {
				BosonOnSphereShort TmpDestinationHilbertSpace(SubsystemNbrParticles, SubsystemTotalLz, LzMax);
				ParticleOnSphereSquareTotalMomentumOperator OperMomentum (&TmpDestinationHilbertSpace, LzMax);
				RealVector TmpEigenstate(1);
				TmpEigenstate[0] = 1.0;
				double TmpSqrMomentum = (OperMomentum.MatrixElement(TmpEigenstate, TmpEigenstate)).Re;
				double TmpMomentum = 0.0;
				if (TmpSqrMomentum > 0.0)
				  TmpMomentum = (0.5 * (sqrt ((4.0 * TmpSqrMomentum) + 1.0) - 1.0));
				DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << " " << TmpSqrMomentum << " " << TmpMomentum << endl;
				if ((EigenstateFlag == true) && ((Manager.GetBoolean("lza-filter") == false) || (FilterLza == SubsystemTotalLz)) && ((NbrEigenstates == 0) || (NbrEigenstates > 0)))
				  {
				    char* TmpEigenstateName = new char[512];
				    RealVector TmpEigenstate(1);
				    TmpEigenstate[0] = 1.0;
				    sprintf (TmpEigenstateName,
					     "bosons_particlereduceddensity_na_%d_n_%d_2s_%d_lz_%d.0.vec",
					     SubsystemNbrParticles, NbrParticles, LzMax, SubsystemTotalLz);
				    TmpEigenstate.WriteVector(TmpEigenstateName);
				    delete[] TmpEigenstateName;
				  }
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
	  File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
	  cout << "trace = " << DensitySum << endl;
	}
      File.close();
    }
  else
    {
      if (Manager.GetBoolean("check-size") == true)
	{
	  long MaxEntanglementMatrix = 0;
	  long MaxEntanglementMatrixRow = 0;
	  long MaxEntanglementMatrixColumn = 0;
	  long MaxDensityMatrix = 0;
	  long Size = 0l;
	  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	    {
	      int SubsystemMaxTotalLz = SubsystemNbrParticles * LzMax;
	      int SubsystemTotalLz = -SubsystemMaxTotalLz; 
	      if (LargestLSector == true)
		{
		  if (((LzMax * NbrParticles) & 1) == 0)
		    {
		      SubsystemTotalLz = 0;
		      SubsystemMaxTotalLz = 0;
		    }
		  else
		    {
		      SubsystemTotalLz = 1;
		      SubsystemMaxTotalLz = 1;
		    }
		}
	      if (PositiveLzSectors == true)
		{
		  if (((LzMax * NbrParticles) & 1) == 0)
		    {
		      SubsystemTotalLz = 0;
		    }
		  else
		    {
		      SubsystemTotalLz = 1;
		    }
		}
	      for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
		{
		  long SizeA = BosonEvaluateHilbertSpaceDimension(SubsystemNbrParticles, LzMax, SubsystemTotalLz);
		  long SizeB = BosonEvaluateHilbertSpaceDimension(NbrParticles - SubsystemNbrParticles, LzMax, TotalLz[0] - SubsystemTotalLz);
		  Size += SizeA * SizeB;
		  if (MaxEntanglementMatrix < (SizeA * SizeB))
		    {
		      MaxEntanglementMatrix = (SizeA * SizeB);
		      MaxEntanglementMatrixRow = SizeA;
		      MaxEntanglementMatrixColumn = SizeB;		      
		    }
		  if (MaxDensityMatrix < SizeA)
		    {
		      MaxDensityMatrix = SizeA;
		    }
		}
	    }
	  Size <<= 3;
	  cout << "total required storage = ";
	  PrintMemorySize(cout, Size) << endl;    
	  cout << "largest file size = ";
	  PrintMemorySize(cout, MaxEntanglementMatrix << 3) << "(matrix  " << MaxEntanglementMatrixRow << "x" << MaxEntanglementMatrixColumn << ")" << endl;
	  cout << "largest density matrix = " << MaxDensityMatrix << "x" << MaxDensityMatrix << endl;
	  return 0;
	}
      for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  int SubsystemMaxTotalLz = SubsystemNbrParticles * LzMax;
	  int SubsystemTotalLz = -SubsystemMaxTotalLz; 
	  if (LargestLSector == true)
	    {
	      if (((LzMax * NbrParticles) & 1) == 0)
		{
		  SubsystemTotalLz = 0;
		  SubsystemMaxTotalLz = 0;
		}
	      else
		{
		  SubsystemTotalLz = 1;
		  SubsystemMaxTotalLz = 1;
		}
	    }
	  if (PositiveLzSectors == true)
	    {
	      if (((LzMax * NbrParticles) & 1) == 0)
		{
		  SubsystemTotalLz = 0;
		}
	      else
		{
		  SubsystemTotalLz = 1;
		}
	    }
	  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
	    {
	      RealMatrix TmpEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[0], EntanglementMatrixForRealCut);
	      if (TmpEntanglementMatrix.GetNbrRow() > 0)
		{
		  char* TmpExtension = new char[256];
		  if (EntanglementMatrixForRealCut == true)
		    {
		      sprintf (TmpExtension, "pem_realcut_na_%d_lza_%d.mat", SubsystemNbrParticles, SubsystemTotalLz);
		    }
		  else
		    {
		      sprintf (TmpExtension, "pem_na_%d_lza_%d.mat", SubsystemNbrParticles, SubsystemTotalLz);
		    }
		  char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", TmpExtension);
		  if (TmpFileName == 0)
		    {
		      cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
		      return 0;
		    }
		  TmpEntanglementMatrix.WriteMatrix(TmpFileName);
		  delete[] TmpFileName;
		  delete[] TmpExtension;
		}
	    }
	}
    }
}

// evaluate Hilbert space dimension for bosons
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

long BosonEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  return BosonShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, (totalLz + lzMax * nbrBosons) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for bosons
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long BosonShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons == 0) || ((nbrBosons * lzMax) < totalLz))
    return 0l;
  if (((nbrBosons * lzMax) == totalLz) || (lzMax == 0) || (totalLz == 0))
    {
      return 1l;
    }
  long TmpDim = 0;
  while ((totalLz >= 0) && (nbrBosons > 0))
    {
      TmpDim += BosonShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz);
      --nbrBosons;
      totalLz -= lzMax;
    }
  return TmpDim;
}
