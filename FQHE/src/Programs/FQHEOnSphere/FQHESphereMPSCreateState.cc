#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnSpherePTruncated.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/BosonOnDiskWithSU2Spin.h"

#include "MathTools/ClebschGordanCoefficients.h"
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
  OptionManager Manager ("FQHESphereMPSCreateState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager;

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new BooleanOption  ('\n', "full-basis", "express the final vector in the full Haldane basis for the given root partition");
  (*SystemGroup) += new BooleanOption  ('\n', "use-padding", "root partitions use the extra zero padding");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "precalculation-blocksize", " indicates the size of the block (i.e. number of B matrices) for precalculations", 1);
  (*OutputGroup) += new BooleanOption ('n', "normalize-sphere", "express the MPS in the normalized sphere basis");
  (*OutputGroup) += new BooleanOption ('\n', "normalize-jack", "use the Jack normalization, forcing the first component to be 1");
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the MPS state into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the MPS state into a text file");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMPSCreateState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  int TotalLz = 0;
  int* ReferenceState = 0;
  char* OutputFileName = Manager.GetString("bin-output");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  if ((OutputTxtFileName == 0) && (OutputFileName == 0))
    {
      cout << "error, an output file (binary or text) has to be provided" << endl;
      return 0;
    }

  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
    return -1;

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");
  double AspectRatio = Manager.GetDouble("aspect-ratio");
  double Kappa = 0.0;
  if (CylinderFlag)
    {
       Kappa = (2.0 * M_PI)/sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) * AspectRatio);
       cout<<"Cylinder geometry, kappa= "<<Kappa<<endl;
    }

  int NbrQuasiholes = 0;
  Complex* QuasiholePositions = 0;
  if (Manager.GetString("with-quasiholes") != 0)
    {
      MultiColumnASCIIFile InputQuasiholePosition;
      if (InputQuasiholePosition.Parse(Manager.GetString("with-quasiholes")) == false)
	{
	  InputQuasiholePosition.DumpErrors(cout) << endl;
	  return -1;
	}
      QuasiholePositions = InputQuasiholePosition.GetAsComplexArray(0);
      NbrQuasiholes = InputQuasiholePosition.GetNbrLines();
   }

  ParticleOnSphere* Space = 0;
  if (Manager.GetBoolean("boson") == true)
    {
      Space = new BosonOnSpherePTruncated(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
      return 0;
    }
  else
    {
#ifdef __64_BITS__
	  if (NbrFluxQuanta <= 62)
#else
	  if (NbrFluxQuanta <= 30)
#endif
	    {
	      Space = new FermionOnSpherePTruncated(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
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


  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;
  if (Manager.GetBoolean("use-padding") == true)
    {
      if (Manager.GetBoolean("k-2") == true)
	{
	  if ((Manager.GetInteger("r-index") & 1) == 0)
	    MPSRowIndex = Manager.GetInteger("p-truncation") + (Manager.GetInteger("r-index") / 2);
	  else
	    MPSRowIndex = 2 * Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index") - 1;
	}
      else
	{
	  if (Manager.GetBoolean("rr-3") == true)
	    {
	      MPSRowIndex = 3 * (Manager.GetInteger("p-truncation") + 1);
	    }
	  else
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + ((Manager.GetInteger("laughlin-index") - 1) / 2);
	    }
	}
      MPSColumnIndex = MPSRowIndex;
    }
  else
    {
      if (Manager.GetBoolean("k-2") == true)
	{
	  if ((Manager.GetInteger("r-index") & 1) == 0)
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index");
	      MPSColumnIndex = Manager.GetInteger("p-truncation");
	    }
	  else
	    {
	      MPSRowIndex = 2 * (Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index"));
	      MPSColumnIndex = 2 * Manager.GetInteger("p-truncation");
	    }
	}
      else
	{
	  if (Manager.GetBoolean("rr-3") == true)
	    {
	      if (Manager.GetBoolean("quasihole-sector") == false)
		{
		  MPSRowIndex = 3 * (Manager.GetInteger("p-truncation") + 2);
		  MPSColumnIndex = 3 * Manager.GetInteger("p-truncation");
		}
	      else
		{
		  MPSRowIndex = 9 * Manager.GetInteger("p-truncation") + 9;
		  MPSColumnIndex = 9 * Manager.GetInteger("p-truncation") + 6;
		}		 
	    }
	  else
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + (Manager.GetInteger("laughlin-index") - 1);
	      MPSColumnIndex = Manager.GetInteger("p-truncation");
	    }
	}
    }

  SparseRealMatrix* SparseBMatrices = MPSMatrix->GetMatrices();

  cout << "B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;

  SparseComplexMatrix* SparseQuasiholeBMatrices = 0;
  if (NbrQuasiholes > 0)
    SparseQuasiholeBMatrices = MPSMatrix->GetQuasiholeMatrices(NbrQuasiholes, QuasiholePositions);
  RealVector State ;
  ComplexVector ComplexState ;
  if (NbrQuasiholes > 0)
    ComplexState = ComplexVector(Space->GetHilbertSpaceDimension(), true);
  else
    State = RealVector(Space->GetHilbertSpaceDimension(), true);


  if (NbrQuasiholes > 0)
    {
      FQHEMPSCreateStateOperation Operation(Space, SparseBMatrices, SparseQuasiholeBMatrices, NbrQuasiholes, &ComplexState, MPSRowIndex, MPSColumnIndex,
					    Manager.GetInteger("precalculation-blocksize"));
      Operation.ApplyOperation(Architecture.GetArchitecture());
    }
  else
    {
      FQHEMPSCreateStateOperation Operation(Space, SparseBMatrices, &State, MPSRowIndex, MPSColumnIndex,
					    Manager.GetInteger("precalculation-blocksize"));
      Operation.ApplyOperation(Architecture.GetArchitecture());
    }

  if (Architecture.GetArchitecture()->CanWriteOnDisk() == true)
    {
      if (Manager.GetBoolean("normalize-sphere"))
	{
	  if (NbrQuasiholes == 0)
	    {
	      Space->ConvertFromUnnormalizedMonomial(State);
	    }
	}
      else
	{
	  if (CylinderFlag == true)
	    {
	      if (NbrQuasiholes == 0)
		{
		  State /= State.Norm();
		}
	      else
		{
		  ComplexState /= ComplexState.Norm();
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("normalize-jack") == true)
		{
		  if (NbrQuasiholes == 0)
		    {
		      State /= State[0];
		    }
		  else
		    {
		      ComplexState /= ComplexState[0];
		    }
		}
	    
	    }
	}
      if (Manager.GetBoolean("full-basis") == true)  
	{
	  FermionOnSphereHaldaneBasis* SpaceHaldane = 0;
	  if (Manager.GetBoolean("boson") == true)
	    {
	      cout << "bosons are not yet implemented" << endl;
	      return 0;
	    }
	  else
	    {
#ifdef __64_BITS__
	      if (NbrFluxQuanta <= 62)
#else
		if (NbrFluxQuanta <= 30)
#endif
		  {
		    SpaceHaldane = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);
		  }
		else
		  {
#ifdef __128_BIT_LONGLONG__
		    if (NbrFluxQuanta <= 126)
#else
		      if (NbrFluxQuanta <= 62)
#endif
			{
			  SpaceHaldane = 0;//new FermionOnSpherePTruncatedLong(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
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
	  
	  RealVector NewState;
	  NewState = ((FermionOnSpherePTruncated*) Space)->ConvertToHaldaneBasis(State, *SpaceHaldane);
	  
	  if (OutputTxtFileName != 0)
	    {
	      ofstream File;
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);	
	      for (long i = 0; i < SpaceHaldane->GetLargeHilbertSpaceDimension(); ++i)
		{
		  NewState.PrintComponent(File, i) << " ";
		  Space->PrintStateMonomial(File, i) << endl;
		}
	      File.close();
	    }
	  if (OutputFileName != 0)
	    {
	      NewState.WriteVector(OutputFileName);
	    }
	}
      else 
	{ 
	  if (NbrQuasiholes > 0)
	    {
	      if (OutputTxtFileName != 0)
		{
		  ofstream File;
		  File.open(OutputTxtFileName, ios::binary | ios::out);
		  File.precision(14);	
		  for (long i = 0; i < Space->GetLargeHilbertSpaceDimension(); ++i)
		    {
		      ComplexState.PrintComponent(File, i) << " ";
		      Space->PrintStateMonomial(File, i) << endl;
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
		      Space->PrintStateMonomial(File, i) << endl;
		    }
		  File.close();
		}
	      if (OutputFileName != 0)
		{
		  State.WriteVector(OutputFileName);
		}
	    }
	}
    }
  
  return 0;
}

