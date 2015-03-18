#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"

#include "MathTools/FactorialCoefficient.h"
#include "MathTools/LongRational.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHEMPSCreateStateOperation.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Matrix/SparseComplexMatrix.h"
#include "Matrix/SparseRealMatrix.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;



int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETorusMPSCreateState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager (false, true);

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-fluxquanta", "set the total number of flux quanta and deduce the number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption ('y', "ky-momentum", "constraint on the total momentum modulo the maximum momentum", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "precalculation-blocksize", " indicates the size of the block (i.e. number of B matrices) for precalculations", 1);
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the MPS state into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the MPS state into a text file");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusMPSCreateState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  int TotalKy = 0;
  bool TwistedTorusFlag = false;
  char* OutputFileName = Manager.GetString("bin-output");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  if ((OutputTxtFileName == 0) && (OutputFileName == 0))
    {
      cout << "error, an output file (binary or text) has to be provided" << endl;
      return 0;
    }

  if (Manager.GetInteger("nbr-fluxquanta") <= 0)
    {
      return -1;
    }
  if (Manager.GetInteger("nbr-fluxquanta") > 0)
    {
      NbrFluxQuanta = Manager.GetInteger("nbr-fluxquanta");
    }
  NbrParticles = NbrFluxQuanta / 3;
  TotalKy = Manager.GetInteger("ky-momentum");
  
  double AspectRatio = Manager.GetDouble("aspect-ratio");

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  ParticleOnTorus* Space = 0;
  if (Manager.GetBoolean("boson") == true)
    {
      Space = new BosonOnTorusShort(NbrParticles, NbrFluxQuanta, TotalKy);
    }
  else
    {
#ifdef __64_BITS__
      if (NbrFluxQuanta <= 62)
#else
	if (NbrFluxQuanta <= 30)
#endif
	  {
	    Space = new FermionOnTorus(NbrParticles, NbrFluxQuanta, TotalKy);
	  }
	else
	  {
#ifdef __128_BIT_LONGLONG__
	    if (NbrFluxQuanta <= 126)
#else
	      if (NbrFluxQuanta <= 62)
#endif
		{
		  Space = 0;//new FermionOnSpherePTruncatedLong(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
		}
	      else
		{
#ifdef __128_BIT_LONGLONG__
		  cout << "cannot generate an Hilbert space when nbr-flux > 126" << endl;
#else
		  cout << "cannot generate an Hilbert space when nbr-flux > 62" << endl;
#endif
		  return 0;
		}
	  }
      
    }
  
  cout << "Hilbert space dimension : " << Space->GetLargeHilbertSpaceDimension() << endl;


  SparseRealMatrix* SparseBMatrices = MPSMatrix->GetMatrices();
  cout << "B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;
  SparseRealMatrix StringMatrix = MPSMatrix->GetTorusStringMatrix(NbrParticles);


  int NbrMPSSumIndices = SparseBMatrices[0].GetNbrRow();
  int* MPSSumIndices =  new int[NbrMPSSumIndices];
  for  (int i = 0; i < NbrMPSSumIndices; ++i)
    MPSSumIndices[i] = i;

  RealVector State ;
  ComplexVector ComplexState ;
  if (TwistedTorusFlag == true)
    ComplexState = ComplexVector(Space->GetHilbertSpaceDimension(), true);
  else
    State = RealVector(Space->GetHilbertSpaceDimension(), true);

  if (TwistedTorusFlag == true)
    {
//       FQHEMPSCreateStateOperation Operation(Space, SparseBMatrices, SparseQuasiholeBMatrices, NbrQuasiholes, &ComplexState, MPSRowIndex, MPSColumnIndex,
// 					    Manager.GetInteger("precalculation-blocksize"));
//       Operation.ApplyOperation(Architecture.GetArchitecture());
    }
  else
    {
      Space->CreateStateFromMPSDescription(SparseBMatrices, StringMatrix, State, MPSSumIndices, NbrMPSSumIndices, 1l, 0l, Space->GetLargeHilbertSpaceDimension());
//       FQHEMPSCreateStateOperation Operation(Space, SparseBMatrices, &State, MPSRowIndex, MPSColumnIndex,
// 					    Manager.GetInteger("precalculation-blocksize"));
//       Operation.ApplyOperation(Architecture.GetArchitecture());
    }

  if (Architecture.GetArchitecture()->CanWriteOnDisk() == true)
    {
      if (TwistedTorusFlag == true)
	{
	  ComplexState /= ComplexState.Norm();
	}
      else
	{
	  State /= State.Norm();
	}
      if (TwistedTorusFlag == true)
	{
	  if (OutputTxtFileName != 0)
	    {
	      ofstream File;
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);	
	      for (long i = 0; i < Space->GetLargeHilbertSpaceDimension(); ++i)
		{
		  ComplexState.PrintComponent(File, i) << " ";
		  Space->PrintState(File, i) << endl;
		}
	      File.close();
	    }
	  if (OutputFileName != 0)
	    {
	      ComplexState.WriteVector(OutputFileName);
	    }
	}
      else
	{
	  if (OutputTxtFileName != 0)
	    {
	      ofstream File;
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);	
	      for (long i = 0; i < Space->GetLargeHilbertSpaceDimension(); ++i)
		{
		  State.PrintComponent(File, i) << " ";
		  Space->PrintState(File, i) << endl;
		}
	      File.close();
	    }
	  if (OutputFileName != 0)
	    {
	      State.WriteVector(OutputFileName);
	    }
	}
    }
  
  return 0;
}

