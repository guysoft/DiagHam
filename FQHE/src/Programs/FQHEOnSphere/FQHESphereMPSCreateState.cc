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

#include "Tools/FQHEMPS/FQHEMPSLaughlinMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3Matrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHEMPSCreateStateOperation.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

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
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

//   (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
//   (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 20);
//   (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "p-truncation", "truncation level", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-index", "index of the Laughlin state to generate", 3);
  (*SystemGroup) += new BooleanOption  ('\n', "k-2", "consider the (k=2,r) series of clustered states");
 (*SystemGroup) += new BooleanOption  ('\n', "rr-3", "consider the k= 3 Read-Rezayi state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "r-index", "r index of the (k,r) clustered state", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "full-basis", "express the final vector in the full Haldane basis for the given root partition");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "precalculation-blocksize", " indicates the size of the block (i.e. number of B matrices) for precalculations", 1);
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the MPS state into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the MPS state into a text file");
  (*OutputGroup) += new BooleanOption ('n', "normalize-sphere", "express the MPS in the normalized sphere basis");
  (*OutputGroup) += new BooleanOption ('c', "normalize-cylinder", "express the MPS in the normalized cylinder basis");
  (*OutputGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the cylinder", 1);

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

  int LaughlinIndex = Manager.GetInteger("laughlin-index");
  int RIndex = Manager.GetInteger("r-index");
  int NbrBMatrices = 2;
  int MatrixElement = 0;
  if (Manager.GetBoolean("boson"))
    {
      NbrBMatrices = NbrParticles + 1;
    }
  AbstractFQHEMPSMatrix* MPSMatrix = 0; 
  if (Manager.GetBoolean("k-2") == true)
    {
      MPSMatrix = new FQHEMPSClustered2RMatrix(Manager.GetInteger("r-index"), 2, Manager.GetInteger("p-truncation"), NbrBMatrices,
					       CylinderFlag, Kappa);
      if ((Manager.GetInteger("r-index") & 1) == 0)
	MatrixElement = Manager.GetInteger("p-truncation") + (Manager.GetInteger("r-index") / 2);
      else
	MatrixElement = 2 * Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index") - 1;
    }
  else
    {
      if (Manager.GetBoolean("rr-3") == true)
	{
	  MPSMatrix = new FQHEMPSReadRezayi3Matrix(2, Manager.GetInteger("p-truncation"), NbrBMatrices,
						   CylinderFlag, Kappa);
	  MatrixElement = 3 * (Manager.GetInteger("p-truncation") + 1);
	}
      else
	{
	  MPSMatrix = new FQHEMPSLaughlinMatrix(LaughlinIndex, Manager.GetInteger("p-truncation"), NbrBMatrices,
						CylinderFlag, Kappa);
	  MatrixElement = Manager.GetInteger("p-truncation") + ((LaughlinIndex - 1) / 2);
	}
    }

  SparseRealMatrix* SparseBMatrices = MPSMatrix->GetMatrices();

  cout << "B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;

  RealVector State (Space->GetHilbertSpaceDimension(), true);


  FQHEMPSCreateStateOperation Operation(Space, SparseBMatrices, &State, MatrixElement,
					Manager.GetInteger("precalculation-blocksize"));
  Operation.ApplyOperation(Architecture.GetArchitecture());

  if (Architecture.GetArchitecture()->CanWriteOnDisk() == true)
    {
      if (Manager.GetBoolean("normalize-sphere"))
	Space->ConvertFromUnnormalizedMonomial(State);
      if (CylinderFlag)
	State /= State.Norm();
      
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
  
  return 0;
}

