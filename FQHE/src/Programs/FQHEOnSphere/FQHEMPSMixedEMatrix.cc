#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/FermionOnSphereMPSWrapper.h"
#include "HilbertSpace/FermionOnCylinderMPSWrapper.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"
#include "Hamiltonian/TensorProductSparseMatrixSelectedBlockHamiltonian.h"

#include "LanczosAlgorithm/BasicArnoldiAlgorithm.h"
#include "LanczosAlgorithm/BasicArnoldiAlgorithmWithDiskStorage.h"

#include "Matrix/SparseRealMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/FQHEMPSEMatrixMainTask.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

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
  cout.precision(14); 
  
  OptionManager Manager ("FQHEMPSEMatrix" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager(true);

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  OptionGroup* ArnoldiGroup  = new OptionGroup ("Arnoldi options");
  Architecture.AddOptionGroup(&Manager);
  Manager += ArnoldiGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new BooleanOption  ('\n', "diagonal-block", "consider only the block diagonal in P, CFT sector and Q");
  (*SystemGroup) += new BooleanOption  ('\n', "right-eigenstates", "compute the right eigenstates");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "theta", "angle between the components");
  (*SystemGroup) += new BooleanOption  ('\n', "left-eigenstates", "compute the left eigenstates");
  (*SystemGroup) += new BooleanOption  ('\n', "fixed-parity", "compute the left eigenstates");

  
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "ematrix-memory", "amount of memory that can used for precalculations of the E matrix (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
  (*OutputGroup) += new BooleanOption  ('\n', "show-ematrix", "show the transfer matrix");
  (*OutputGroup) += new BooleanOption  ('\n', "ematrix-dimonly", "only compute the dimension of the transfer matrix");
  (*ArnoldiGroup) += new SingleIntegerOption  ('\n', "full-diag", 
					       "maximum Hilbert space dimension for which full diagonalization is applied", 1000);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "disk", "enable disk storage for the Arnoldi algorithm", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "resume", "resume from disk datas", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Arnoldi iteration", false); 
  (*ArnoldiGroup) += new BooleanOption  ('\n', "power-method", "use the power method instead of the Arnoldi algorithm. A single eigenvalue is computed (the one with the largest real part)", false); 
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "arnoldi-memory", "amount of memory when using the Arnoldi algorithm (in Mb)", 500); 
  (*ArnoldiGroup) += new BooleanOption  ('\n', "implicitly-restarted", "use the implicitly restarted Arnoldi algorithm", false); 
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "nbr-excited", "number of eigenvalues to compute above the groundstate", 0);
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "nbr-eigenstates", "number of eigenvalues to compute (if set to zero, this number will be deduced from the state and nbr-excited)", 0);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "sort-real", "sort the eigenvalues with respect to their real part instead of their norm", false); 
  (*ArnoldiGroup) += new SingleDoubleOption ('\n', "arnoldi-precision", "define Arnoldi precision for eigenvalues (0 if automatically defined by the program)", 0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEMPSMixedEMatrix -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0; 
  int NbrFluxQuanta = 1;
  int TotalLz = 0;
  double Theta =  Manager.GetDouble("theta");
  bool ParityFlag = Manager.GetBoolean("fixed-parity");
  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");

  int LandauLevel = 0;

  AbstractFQHEMPSMatrix* MPSLeftMatrix = MPSMatrixManager.GetLeftMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  AbstractFQHEMPSMatrix* MPSRightMatrix = MPSMatrixManager.GetRightMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 

  
  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;

  SparseRealMatrix* TmpSparseBMatrices = MPSLeftMatrix->GetMatrices();
  SparseRealMatrix* TmpSparseRightBMatrices = MPSRightMatrix->GetMatrices();
  int NbrBMatrices = 0;

  for (int i = 0; i < MPSLeftMatrix->GetNbrMatrices(); ++i)
    {
      if (SearchInUnsortedArray<unsigned long>(MPSLeftMatrix->GetPhysicalIndices()[i], MPSRightMatrix->GetPhysicalIndices(), MPSRightMatrix->GetNbrMatrices()) >= 0)
	++NbrBMatrices;
    }
  SparseRealMatrix* SparseBMatrices = new SparseRealMatrix[NbrBMatrices];
  SparseRealMatrix* SparseRightBMatrices = new SparseRealMatrix[NbrBMatrices];
  NbrBMatrices = 0;
  for (int i = 0; i < MPSLeftMatrix->GetNbrMatrices(); ++i)
    {
      int TmpIndex = SearchInUnsortedArray<unsigned long>(MPSLeftMatrix->GetPhysicalIndices()[i], MPSRightMatrix->GetPhysicalIndices(), MPSRightMatrix->GetNbrMatrices());
      if (TmpIndex >= 0)
	{
	  SparseBMatrices[NbrBMatrices] = TmpSparseBMatrices[i];
	  SparseRightBMatrices[NbrBMatrices] = TmpSparseRightBMatrices[TmpIndex];
	  ++NbrBMatrices;
	}
    }

  cout << "Left B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;
  cout << "Right B matrix size = " << SparseRightBMatrices[0].GetNbrRow() << "x" << SparseRightBMatrices[0].GetNbrColumn() << endl;

  char* StateName = 0;
  if (strcmp(MPSLeftMatrix->GetName(), MPSRightMatrix->GetName()) == 0)
    {
      StateName = new char [strlen(MPSLeftMatrix->GetName()) + 1];
      strcpy(StateName, MPSLeftMatrix->GetName());
    }
  else
    {
      StateName = new char [strlen(MPSLeftMatrix->GetName()) + strlen(MPSRightMatrix->GetName()) + 2];
      sprintf (StateName, "%s_%s", MPSLeftMatrix->GetName(), MPSRightMatrix->GetName());
    }

  int NbrEigenstates = MPSLeftMatrix->GetTransferMatrixLargestEigenvalueDegeneracy()*2;
  NbrEigenstates += Manager.GetInteger("nbr-excited");
  if (Manager.GetInteger("nbr-eigenstates") > 0)
    {
      NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
    }
  
  double EnergyShift = 0.0;
  if (Manager.GetBoolean("power-method") == true)
    {
      EnergyShift = 10.0;
      NbrEigenstates = 1;
    }

  char* PrefixOutputFileName = 0;
  char* OutputFileName = 0;
  char * AngleString = new char [50];
  sprintf(AngleString,"theta_%f",Theta);
  char * ParityString = new char [50];
  if (ParityFlag)
    {
      sprintf(ParityString,"%s","fixedparity_");
    }
  else
    sprintf(ParityString,"%s","");


  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
      if (GetExtensionFromFileName(OutputFileName) != 0)
	{
	  int TmpLength = strlen(OutputFileName) - strlen(GetExtensionFromFileName(OutputFileName));
	  PrefixOutputFileName = new char[TmpLength + 1];
	  strncpy(PrefixOutputFileName, OutputFileName, TmpLength);
	  PrefixOutputFileName[TmpLength] = '\0';
	}
    }
  else
    {
      PrefixOutputFileName  = new char [1024];
      if (CylinderFlag == true)
	{
	  if (Manager.GetBoolean("diagonal-block"))
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_cylinder_%s_%s%s_perimeter_%f_plevel_%ld_maxocc_%ld", StateName,ParityString,AngleString,  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"), 
			  Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_cylinder_%s_%s%s_perimeter_%f_plevel_%ld", StateName,ParityString,AngleString,  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  sprintf(PrefixOutputFileName, "ematrix_cylinder_%s_%s%s_perimeter_%f_plevel_%ld_maxocc_%ld", StateName,ParityString,AngleString,
			  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"), 
			  Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_cylinder_%s_%s%s_perimeter_%f_plevel_%ld", StateName,ParityString,AngleString,
			  MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("diagonal-block"))
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_%s_%s%s_plevel_%ld_maxocc_%ld", StateName,ParityString,AngleString,
			  Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_diagblock_%s_%s%s_plevel_%ld", StateName,ParityString,AngleString,
			  Manager.GetInteger("p-truncation"));
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  sprintf(PrefixOutputFileName, "ematrix_%s_%s%s_plevel_%ld_maxocc_%ld", StateName,ParityString,AngleString,
			  Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"));
		}
	      else
		{
		  sprintf(PrefixOutputFileName, "ematrix_%s_%s%s_plevel_%ld", StateName,ParityString,AngleString,
			  Manager.GetInteger("p-truncation"));
		}
	    }
	}
      OutputFileName = new char[strlen(PrefixOutputFileName) + 64];
      sprintf(OutputFileName, "%s.dat", PrefixOutputFileName);
   }


  double Error = 1e-13;
  cout <<"NbrBMatrices = "<<NbrBMatrices<<endl;
  SparseRealMatrix MixedEMatrix;
  SparseRealMatrix*  RightMatrices = new SparseRealMatrix[NbrBMatrices];
  double * Coefficients= new double[NbrBMatrices];

  for (int i =0;i <NbrBMatrices;i++)
    {
      Coefficients[i] = 1;
      RightMatrices[i] = SparseRealMatrix(1,1);
    }

  RightMatrices[0].SetMatrixElement(0,0,sin(Theta));
  RightMatrices[1].SetMatrixElement(0,0,cos(Theta));

  int NbrOrbitals = MPSLeftMatrix->GetNbrOrbitals();
  int NbrStatesPerOrbital = MPSLeftMatrix->GetMaximumOccupation() + 1;
  int NbrStatesPerBlock =  1;
  int NbrRMatrices = 2;
  for (int i = 0; i < NbrOrbitals ; i++)
    NbrStatesPerBlock *= NbrStatesPerOrbital;


  cout << "NbrOrbitals = " << NbrOrbitals << " NbrStatesPerOrbital = " <<NbrStatesPerOrbital<<" NbrStatesPerBlock =" <<NbrStatesPerBlock<<endl; 
  unsigned long * ArrayPhysicalIndice = MPSLeftMatrix->GetPhysicalIndices();
  

  SparseRealMatrix* FusedRMatrices = new SparseRealMatrix [NbrStatesPerBlock];
  
  int NbrUn =0 ;
  int TmpI;
  for(int i =0 ; i < NbrStatesPerBlock; i++)
    { 
      int NbrUn =0 ;
      TmpI = i;
      int Index = SearchInUnsortedArray( (unsigned long)( TmpI %  NbrStatesPerOrbital) , ArrayPhysicalIndice,  NbrRMatrices);
      if (Index <0)
	{
	  FusedRMatrices[i] = SparseRealMatrix(RightMatrices[0].GetNbrRow(),RightMatrices[0].GetNbrColumn());
	}
      else
	{
	  FusedRMatrices[i].Copy(RightMatrices[Index]);
	}
      NbrUn+=TmpI %  NbrStatesPerOrbital;
      TmpI /= NbrStatesPerOrbital;
      for(int p = 1; p < NbrOrbitals ; p++)
	{
	  int Index = SearchInUnsortedArray( (unsigned long)( TmpI %   NbrStatesPerOrbital) , ArrayPhysicalIndice,  NbrRMatrices);
	  if (Index <0)
	    {
	      FusedRMatrices[i].ClearMatrix ();
	    }
	  else
	    {
	      FusedRMatrices[i].Multiply(RightMatrices[Index]);
	    }
          NbrUn+=TmpI %  NbrStatesPerOrbital;
	  TmpI /= NbrStatesPerOrbital;
	}
      if (( NbrUn %2 ==0)&&(ParityFlag))
	Coefficients[i] =0;	
    }

  TensorProductSparseMatrixHamiltonian  * ETransposeHamiltonian =0;
  ETransposeHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, SparseBMatrices, FusedRMatrices, Coefficients, Architecture.GetArchitecture());
  Architecture.GetArchitecture()->SetDimension(SparseBMatrices[0].GetNbrRow());
  
  FQHEMPSEMatrixMainTask TaskLeft(&Manager, ETransposeHamiltonian, NbrEigenstates, false, true, 1e-10, EnergyShift, OutputFileName);
  MainTaskOperation TaskOperationLeft (&TaskLeft);
  TaskOperationLeft.ApplyOperation(Architecture.GetArchitecture());

  return 0;
}

