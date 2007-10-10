#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
// #include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
// #include "HilbertSpace/FermionOnSphereUnlimited.h"
// #include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
// #include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
// #include "HilbertSpace/FermionOnSphereLong.h"
// #include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
// #include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
// #include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereEntanglementEntropy" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-la", "minimum size of the subsystem whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-la", "maximum size of the subsystem whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonEntanglementEntropy -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((SingleStringOption*) Manager["ground-file"])->GetString() == 0)
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereEntanglementEntropy -h" << endl;
      return -1;
    }
  if (IsFile(((SingleStringOption*) Manager["ground-file"])->GetString()) == false)
    {
      cout << "can't open file " << ((SingleStringOption*) Manager["ground-file"])->GetString() << endl;
    }

  bool HaldaneBasisFlag = ((BooleanOption*) Manager["haldane"])->GetBoolean();
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger(); 
#ifdef __LAPACK__
  bool LapackFlag = ((BooleanOption*) Manager["use-lapack"])->GetBoolean();
#endif
  char* DensityMatrixFileName = ((SingleStringOption*) Manager["density-matrix"])->GetString();
  int TotalLz = 0;
  bool Statistics = true;
  if (QHEOnSphereFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["ground-file"])->GetString(),
						  NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << ((SingleStringOption*) Manager["ground-file"])->GetString() << endl;
      return -1;
    }
  if (Statistics == true)
    {
      cout << ((SingleStringOption*) Manager["ground-file"])->GetString() << " is not a bosonic state" << endl;
      return -1;
    }
  if (((SingleIntegerOption*) Manager["total-lz"])->GetInteger() >= 0)
    TotalLz = ((SingleIntegerOption*) Manager["total-lz"])->GetInteger(); 

  if (((NbrParticles * LzMax) & 1) != (TotalLz & 1))
    {
      cout << "incompatible values for nbr-particles, nbr-flux and total-lz" << endl;
      return -1;
    }

  RealVector GroundState;
  if (GroundState.ReadVector (((SingleStringOption*) Manager["ground-file"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["ground-file"])->GetString() << endl;
      return -1;      
    }


  ParticleOnSphere* Space = 0;
  if ((SymmetrizedBasis == false) || (TotalLz != 0))
    Space = new BosonOnSphere (NbrParticles, TotalLz, LzMax);
  else
    {
      Space = new BosonOnSphere (NbrParticles, TotalLz, LzMax);
      BosonOnSphereSymmetricBasis TmpSpace(NbrParticles, LzMax);
      RealVector OutputState = TmpSpace.ConvertToNbodyBasis(GroundState, *((BosonOnSphere*) Space));
      GroundState = OutputState;
    }

  if (Space->GetHilbertSpaceDimension() != GroundState.GetVectorDimension())
    {
      cout << "dimension mismatch between Hilbert space and ground state" << endl;
      return 0;
    }


  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "# l_a    N    Lz    lambda" << endl;
      DensityMatrixFile.close();
    }

  ofstream File;
  if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
    File.open(((SingleStringOption*) Manager["output-file"])->GetString(), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(((SingleStringOption*) Manager["ground-file"])->GetString(), "vec", "ent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << ((SingleStringOption*) Manager["ground-file"])->GetString() << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);
  int MeanSubsystemSize = LzMax >> 1;
  if ((LzMax & 1) != 0)
    ++MeanSubsystemSize;
  if (((SingleIntegerOption*) Manager["max-la"])->GetInteger() > 0)
    {
      MeanSubsystemSize = ((SingleIntegerOption*) Manager["max-la"])->GetInteger();
      if (MeanSubsystemSize > LzMax)
	MeanSubsystemSize = LzMax;
    }
  int SubsystemSize = ((SingleIntegerOption*) Manager["min-la"])->GetInteger();
  if (SubsystemSize < 1)
    SubsystemSize = 1;
  BinomialCoefficients Coefs(MeanSubsystemSize + NbrParticles - 1);
  for (; SubsystemSize <= MeanSubsystemSize; ++SubsystemSize)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      int MaxSubsystemNbrParticles = NbrParticles;
      if (MaxSubsystemNbrParticles > SubsystemSize)
	MaxSubsystemNbrParticles = SubsystemSize;
      long MaximumSize = 0;
      for (int i = 0; i <= NbrParticles; ++i)
	MaximumSize += Coefs(SubsystemSize + i - 1, i);
      double* TmpDensityMatrixEigenvalues = new double [MaximumSize];
      long TmpDensityMatrixEigenvaluePosition = 0;
      for (int SubsystemNbrParticles = 0; SubsystemNbrParticles <= NbrParticles; ++SubsystemNbrParticles)
	{
	  int SubsystemTotalLz = 0;
	  int SubsystemLzMax = SubsystemSize - 1;
	  int SubsystemMaxTotalLz = SubsystemNbrParticles * SubsystemLzMax;
	  SubsystemTotalLz = -SubsystemMaxTotalLz; 
	  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
// 	    if ((TotalLz - SubsystemTotalLz) <= (((LzMax + 1) * (NbrParticles - 2 *SubsystemNbrParticles)) + (SubsystemSize * SubsystemNbrParticles) - 
// 						 ((NbrParticles - SubsystemNbrParticles) * (NbrParticles - SubsystemNbrParticles))))
	      {
		cout << "processing subsystem size=" << SubsystemSize << "  subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
		RealSymmetricMatrix PartialDensityMatrix = Space->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, GroundState);
		if (PartialDensityMatrix.GetNbrRow() > 1)
		  {
		    RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
		    if (LapackFlag == true)
		      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		    else
		      PartialDensityMatrix.Diagonalize(TmpDiag);
#else
		    PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
		    TmpDiag.SortMatrixDownOrder();
		    for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		      TmpDensityMatrixEigenvalues[TmpDensityMatrixEigenvaluePosition++] = TmpDiag[i];
		    if (DensityMatrixFileName != 0)
		      {
			ofstream DensityMatrixFile;
			DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			DensityMatrixFile.precision(14);
			for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			  DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			DensityMatrixFile.close();
		      }
		  }
		else
		  if (PartialDensityMatrix.GetNbrRow() == 1)
		    {
		      double TmpValue = PartialDensityMatrix(0,0);
		      TmpDensityMatrixEigenvalues[TmpDensityMatrixEigenvaluePosition++] = TmpValue;
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
			  DensityMatrixFile.close();
			}		  
		    }
	      }
	}
      EntanglementEntropy = 0.0;
      DensitySum = 0.0;
      cout << "sorting density matrix eigenvalues and computing entanglement entropy" << endl;
      SortArrayDownOrdering(TmpDensityMatrixEigenvalues, TmpDensityMatrixEigenvaluePosition);
      unsigned TmpPos = 0;
      for (; (TmpPos < TmpDensityMatrixEigenvaluePosition) && (DensitySum < 1.0); ++TmpPos)
	{
	  if (TmpDensityMatrixEigenvalues[TmpPos] > 1e-14)
	    {
	      EntanglementEntropy += TmpDensityMatrixEigenvalues[TmpPos] * log(TmpDensityMatrixEigenvalues[TmpPos]);
	      DensitySum += TmpDensityMatrixEigenvalues[TmpPos];
	    }
	}
      double DensitySumError = 0.0;
      for (; TmpPos < TmpDensityMatrixEigenvaluePosition; ++TmpPos)
	DensitySumError += TmpDensityMatrixEigenvalues[TmpPos];
      delete[] TmpDensityMatrixEigenvalues;
      File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << " " << DensitySumError << endl;
    }
  File.close();
}

