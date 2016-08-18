#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpin.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/MultipleVectorHamiltonianMultiplyOperation.h"

#include "Hamiltonian/ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHEOnCylinderFileTools.h"

#include "Vector/RealVector.h"


#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity" , "0.01");
  OptionGroup* SystemGroup  = new OptionGroup("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup("precalculation options");
  OptionGroup* OutputGroup  = new OptionGroup ("ouput options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('\n', "input-state", "name of the file containing the state to convert");
  (*SystemGroup) += new SingleStringOption('\n', "degenerate-states", "name of the file containing a list of states (override input-state)");
  (*SystemGroup) += new SingleStringOption('\n', "use-hilbert", "similar to --degenerate-states using a single line text file (starting with Basis=) instead of a single column text file");
  (*SystemGroup) += new SingleStringOption ('\n', "directory", "use a specific directory for the input data instead of the current one");
  (*SystemGroup) += new SingleStringOption ('\n', "occupation-matrices", "use precomputed occupation matrices to evaluate the integrated charge");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "flux-insertion", "include a flux insertion (in flux quantum unit) along the cylinder axis", 0.0);
  (*OutputGroup) += new BooleanOption ('\n', "realspace-density", "plot the density in real space");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-points", "number of points along the cylinder axis", 400);  
  (*OutputGroup) += new SingleDoubleOption  ('\n', "offset", "additional length along the cylinder axis on each side of the [-Lx/2,Lx/2] region where the density should be computed", 5.0);
  (*OutputGroup) += new BooleanOption ('\n', "disable-binary", "do not export the orbital occupation matrices in binary format");
  (*SystemGroup) += new SingleStringOption ('o', "output", "store charge imbalance eigenvalues in a file");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereQuasiholesWithSpinTimeReversalSymmetryDensity -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int LzMax = 0;
  int TotalSz = 0;
  int TotalLz = 0;
  bool Statistics = true;  
  int KValue = 1;
  int RValue = 2;
  bool UseCylinderFlag = false;
  double Ratio = 0.0;
  double Perimeter = 0.0;
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int NbrInputStates = 0;
  char** InputStateNames = 0;
  char* OutputName = 0;
  ofstream FileChargeImbalance;
  RealSymmetricMatrix** OneBodyMatrixElements = new RealSymmetricMatrix*[2];
      
  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0) 
      && (Manager.GetString("use-hilbert") == 0))
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type FQHESphereQuasiholesWithSpinTimeReversalSymmetryConvertStates -h" << endl;
      return -1;
    }
    
  if (Manager.GetString("output"))
    {
      OutputName = new char[strlen(Manager.GetString("output") + 1)];
      strcpy (OutputName, Manager.GetString("output"));
      FileChargeImbalance.open(OutputName, ios::binary | ios::out); 
      FileChargeImbalance.precision(14); 
      FileChargeImbalance << "# i (Q_L-Q_R)_orb (Q_L-Q_R)_real" << endl;
    }
  if (Manager.GetString("input-state") != 0)
    {
      if (FQHEOnCylinderWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, LzMax, TotalLz, TotalSz, Statistics, Ratio, Perimeter) == false)
	{
	  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, LzMax, TotalLz, TotalSz, Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state")  << endl;
	return -1;
	    }
	}
      else
	{
	  UseCylinderFlag = true;
	}
      NbrInputStates = 1;
      InputStateNames = new char*[NbrInputStates];
      InputStateNames[0] = new char [strlen(Manager.GetString("input-state")) + 1];
      strcpy (InputStateNames[0], Manager.GetString("input-state"));
    }
  else
    {
      if (Manager.GetString("use-hilbert") == 0)
	{
	  MultiColumnASCIIFile DegenerateFile;
	  if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	    {
	      DegenerateFile.DumpErrors(cout);
	      return -1;
	    }
	  NbrInputStates = DegenerateFile.GetNbrLines();
	  InputStateNames = new char*[NbrInputStates];
	  for (int i = 0; i < NbrInputStates; ++i)
	    {
	      InputStateNames[i] = new char [strlen(DegenerateFile(0, i)) + 1];
	      strcpy (InputStateNames[i], DegenerateFile(0, i));
	    }
	}
      else
	{
	  ConfigurationParser DegenerateFile;
	  if (DegenerateFile.Parse(Manager.GetString("use-hilbert")) == false)
	    {
	      DegenerateFile.DumpErrors(cout);
	      return -1;
	    }
	  if (DegenerateFile.GetAsStringArray("Basis", ' ', InputStateNames, NbrInputStates) == false)
	    {
	      return 0;
	    }
	}
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  if (FQHEOnCylinderWithSpinFindSystemInfoFromVectorFileName(InputStateNames[i], NbrParticles, LzMax, TotalLz, TotalSz, Statistics, Ratio, Perimeter) == false)
	    {
	      if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(InputStateNames[i], NbrParticles, LzMax, TotalLz, TotalSz, Statistics) == false)
		{
		  cout << "error while retrieving system parameters from file name " << InputStateNames[i] << endl;
		  return -1;
		}
	    }
	  else
	    {
	      UseCylinderFlag = true;
	    }
	}
    }
  RealDiagonalMatrix TmpImbalanceEigenvalues(NbrInputStates);     
  if (Manager.GetString("occupation-matrices") == 0)
    {

      char* FilePrefix = new char[512];
      
      if (UseCylinderFlag == true)
	{
	  if (Perimeter > 0.0)	
	    {
	      if (Statistics == true)
		{
		  sprintf (FilePrefix, "fermions_cylinder_perimeter_%.6f", Perimeter);
		}
	      else
		{
		  sprintf (FilePrefix, "bosons_cylinder_perimeter_%.6f", Perimeter);
		}
	    }
	  else
	    {
	      if (Statistics == true)
		{
		  sprintf (FilePrefix, "fermions_cylinder_ratio_%.6f", Ratio);
		}
	      else
		{
		  sprintf (FilePrefix, "bosons_cylinder_ratio_%.6f", Ratio);
		}      
	    }
	}
      else
	{
	  if (Statistics == true)
	    {
	      sprintf (FilePrefix, "fermions");
	    }
	  else
	    {
	      sprintf (FilePrefix, "bosons");
	    }
	}
      
      RealVector* InputStates = new RealVector[NbrInputStates];
      if (Manager.GetString("input-state") != 0)
	{
	  if (InputStates[0].ReadVector (Manager.GetString("input-state")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	      return -1;      
	    }
	}
      else
	{
	  if (InputStates[0].ReadVector (InputStateNames[0]) == false)
	    {
	      cout << "can't open vector file " << InputStateNames[0]  << endl;
	      return -1;      
	    }	  
	  for (int i = 1; i < NbrInputStates; ++i)
	    {
	      if (InputStates[i].ReadVector (InputStateNames[i]) == false)
		{
		  cout << "can't open vector file " << InputStateNames[i] << endl;
		  return -1;      
		}	  
	      if (InputStates[0].GetVectorDimension() != InputStates[i].GetVectorDimension())
		{
		  cout << "error, " << InputStateNames[0] << " and " <<  InputStateNames[i] << " don't have the same  dimension (" 
		       << InputStates[0].GetVectorDimension() << " and " << InputStates[i].GetVectorDimension()<< ")" << endl;
		  return -1;
		}
	    }
	}
      
      RealVector* TmpStates = new RealVector[NbrInputStates];
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  TmpStates[i] = RealVector(InputStates[0].GetVectorDimension());
	}      
      
      cout << KValue << " " << RValue << " " << TotalLz << " " << LzMax << " " << NbrParticles << " " << TotalSz << " " << Manager.GetString("directory") << " " << FilePrefix << endl;
      QuasiholeOnSphereWithSpinAndPairing* InputSpace;
      if (NbrParticles == 0)
	{
	  InputSpace = new QuasiholeOnSphereWithSpinAndPairing (KValue, RValue, TotalLz, LzMax, TotalSz, Manager.GetString("directory"), FilePrefix, true);
	}
      else
	{
	  InputSpace = new QuasiholeOnSphereWithSpin (KValue, RValue, TotalLz, LzMax, NbrParticles, TotalSz, Manager.GetString("directory"), FilePrefix);
	}
      Architecture.GetArchitecture()->SetDimension(InputSpace->GetHilbertSpaceDimension());
      
      
      double* OneBodyPotentialPairing = 0;
      double* OneBodyPotentialDummy = new double [LzMax + 1];
      double* OneBodyPotentials = new double [LzMax + 1];
      for (int i = 0; i <= LzMax; ++i)
	{
	  OneBodyPotentialDummy[i] = 0.0;
	  OneBodyPotentials[i] = 0.0;
	}
      
      for (int LayerIndex = 0; LayerIndex <= 1; ++LayerIndex)
	{
	  OneBodyMatrixElements[LayerIndex] = new RealSymmetricMatrix[LzMax + 1];
	  for (int MomentumIndex = 0; MomentumIndex <= LzMax; ++MomentumIndex)
	    {      
	      for (int i = 0; i <= LzMax; ++i)
		{
		  OneBodyPotentials[i] = 0.0;
		}
	      OneBodyPotentials[MomentumIndex] = 1.0;
	      AbstractHamiltonian* Hamiltonian = 0;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      if (LayerIndex == 0)
		{
		  Hamiltonian = new ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing (InputSpace, LzMax, OneBodyPotentials, 
														 OneBodyPotentialPairing,
														 OneBodyPotentialPairing,
														 0.0, 0.0,
														 Architecture.GetArchitecture(), 
														 Memory, false, 0);
		}
	      else
		{
		  Hamiltonian = new ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing (InputSpace, LzMax, OneBodyPotentialPairing, 
														 OneBodyPotentials,
														 OneBodyPotentialPairing,
														 0.0, 0.0,
														 Architecture.GetArchitecture(), 
														 Memory, false, 0);
		}
	      
	      MultipleVectorHamiltonianMultiplyOperation TmpOperation(Hamiltonian, InputStates, TmpStates, NbrInputStates);
	      TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	      
	      RealSymmetricMatrix TmpMatrix (NbrInputStates, true);
	      for (int i = 0; i < NbrInputStates; ++i)
		{
		  for (int j = i; j < NbrInputStates; ++j)
		    {
		      TmpMatrix.SetMatrixElement(i, j, TmpStates[i] * InputStates[j]);
		    }
		}
	      OneBodyMatrixElements[LayerIndex][MomentumIndex] = TmpMatrix;
	      delete Hamiltonian;
	    }
	}
      delete InputSpace;
      
      char* DensityMatrixTextFileName = ReplaceExtensionToFileName(InputStateNames[0], "vec", "density.dat");
      ofstream File;
      File.open(DensityMatrixTextFileName, ios::binary | ios::out); 
      File.precision(14); 
      File << "# layer_index m i j <i|a^+_m a_m|j>" << endl;
      if (Manager.GetBoolean("disable-binary") == true)
	{
	  for (int LayerIndex = 0; LayerIndex <= 1; ++LayerIndex)
	    {
	      for (int MomentumIndex = 0; MomentumIndex <= LzMax; ++MomentumIndex)
		{      
		  for (int i = 0; i < NbrInputStates; ++i)
		    {
		      for (int j = i; j < NbrInputStates; ++j)
			{	  
			  double Tmp;
			  OneBodyMatrixElements[LayerIndex][MomentumIndex].GetMatrixElement(i, j, Tmp);
			  File << LayerIndex << " " << MomentumIndex << " " << i << " " << j << " " << Tmp << endl;
			}
		    }
		}
	    }
	}
      else
	{
	  char* DensityMatrixListMatrixFileName = ReplaceExtensionToFileName(InputStateNames[0], "vec", "density.mat.list");
	  ofstream File2;
	  File2.open(DensityMatrixListMatrixFileName, ios::binary | ios::out); 
	  for (int LayerIndex = 0; LayerIndex <= 1; ++LayerIndex)
	    {
	      for (int MomentumIndex = 0; MomentumIndex <= LzMax; ++MomentumIndex)
		{      
		  for (int i = 0; i < NbrInputStates; ++i)
		    {
		      for (int j = i; j < NbrInputStates; ++j)
			{	  
			  double Tmp;
			  OneBodyMatrixElements[LayerIndex][MomentumIndex].GetMatrixElement(i, j, Tmp);
			  File << LayerIndex << " " << MomentumIndex << " " << i << " " << j << " " << Tmp << endl;
			}
		    }
		  char* TmpExtension = new char[256];
		  if (LayerIndex == 0)
		    {
		      sprintf (TmpExtension, "density_up_m_%d.mat", MomentumIndex);
		    }
		  else
		    {
		      sprintf (TmpExtension, "density_down_m_%d.mat", MomentumIndex);
		    }
		  char* TmpDensityMatrixFileName = ReplaceExtensionToFileName(InputStateNames[0], "vec", TmpExtension);
		  if (OneBodyMatrixElements[LayerIndex][MomentumIndex].WriteMatrix(TmpDensityMatrixFileName) == false)
		    {
		      cout << "can't write " << TmpDensityMatrixFileName << endl;
		      return -1;
		    }
		  File2 << TmpDensityMatrixFileName << endl;
		  delete[] TmpExtension;
		  delete[] TmpDensityMatrixFileName;
		}
	    }
	  File2.close();
	}
      File.close();
    }
  else
    {
      // use precomputed occupation matrices

      RealMatrix InputVectors;
      if (Manager.GetString("input-state") != 0)
	{
	  RealVector* TmpVectors = new RealVector[1];
	  if (TmpVectors[0].ReadVector (Manager.GetString("input-state")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	      return -1;      
	    }
	  InputVectors = RealMatrix(TmpVectors, 1);
	  NbrInputStates = 1;
	}
      else
	{
	  MultiColumnASCIIFile DegenerateFile;
	  if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	    {
	      DegenerateFile.DumpErrors(cout);
	      return -1;
	    }
	  NbrInputStates = DegenerateFile.GetNbrLines();
	  RealVector* TmpVectors = new RealVector[NbrInputStates];
	  for (int i = 0; i < NbrInputStates; ++i)
	    {
	      if (TmpVectors[i].ReadVector(DegenerateFile(0, i)) == false)
		{
		  cout << "can't open vector file " <<  DegenerateFile(0, i)<< endl;
		  return -1;      
		}
	    }
	  InputVectors = RealMatrix(TmpVectors, NbrInputStates);
	}  

      MultiColumnASCIIFile OccupationMatrixListFile;
      if (OccupationMatrixListFile.Parse(Manager.GetString("occupation-matrices")) == false)
	{
	  OccupationMatrixListFile.DumpErrors(cout);
	  return -1;
	}
      LzMax = (OccupationMatrixListFile.GetNbrLines() / 2) - 1;
      for (int LayerIndex = 0; LayerIndex <= 1; ++LayerIndex)
	{
	  OneBodyMatrixElements[LayerIndex] = new RealSymmetricMatrix[LzMax + 1];
	  for (int MomentumIndex = 0; MomentumIndex <= LzMax; ++MomentumIndex)
	    {      
	      RealSymmetricMatrix TmpMatrix;
	      if (TmpMatrix.ReadMatrix(OccupationMatrixListFile(0, MomentumIndex + (LayerIndex * (LzMax + 1)))) == false)
		{
		  return -1;
		}
	      RealSymmetricMatrix* TmpMatrix2 = (RealSymmetricMatrix*) TmpMatrix.Conjugate(InputVectors);
	      OneBodyMatrixElements[LayerIndex][MomentumIndex] = *TmpMatrix2;
	      delete TmpMatrix2;
	    }
	}
    }

  cout << "Integrated charge in the upper layer / lower layer:" << endl;
  RealSymmetricMatrix TotalChargeImbalance(NbrInputStates, true);
  for (int MomentumIndex = 0; MomentumIndex <= LzMax; ++MomentumIndex)
    {      
      if (MomentumIndex < (LzMax + 1) / 2)
	{
	  TotalChargeImbalance += OneBodyMatrixElements[0][MomentumIndex];
	  TotalChargeImbalance += OneBodyMatrixElements[1][MomentumIndex];
	}
      else
	{
	  TotalChargeImbalance -= OneBodyMatrixElements[0][MomentumIndex];
	  TotalChargeImbalance -= OneBodyMatrixElements[1][MomentumIndex];
	}
      
      RealSymmetricMatrix TotalChargeUpLayer (NbrInputStates, true);
      RealSymmetricMatrix TotalCharge (NbrInputStates, true);
      for (int MomentumIndex2 = 0; MomentumIndex2 <= MomentumIndex; ++MomentumIndex2)
	{      
	  TotalChargeUpLayer += OneBodyMatrixElements[0][MomentumIndex2];
	  TotalCharge += OneBodyMatrixElements[0][MomentumIndex2];
	}
      RealDiagonalMatrix TmpUpEigenvalues(NbrInputStates);
      TotalChargeUpLayer.LapackDiagonalize(TmpUpEigenvalues);
      RealSymmetricMatrix TotalChargeDownLayer (NbrInputStates, true);
      for (int MomentumIndex2 = 0; MomentumIndex2 <= MomentumIndex; ++MomentumIndex2)
	{      
	  TotalChargeDownLayer += OneBodyMatrixElements[1][LzMax - MomentumIndex2];
	  TotalCharge += OneBodyMatrixElements[1][LzMax - MomentumIndex2];
	}
      RealDiagonalMatrix TmpDownEigenvalues(NbrInputStates);
      TotalChargeDownLayer.LapackDiagonalize(TmpDownEigenvalues);
      RealDiagonalMatrix TmpEigenvalues(NbrInputStates);
      TotalCharge.LapackDiagonalize(TmpEigenvalues);
      cout << MomentumIndex;
      for (int i = 0; i < NbrInputStates; ++i)
	cout << " " << TmpUpEigenvalues[i] << " " << TmpDownEigenvalues[i] << " " << TmpEigenvalues[i];
      cout << endl;
    }
  
   TotalChargeImbalance.LapackDiagonalize(TmpImbalanceEigenvalues);
   for (int i = 0; i < NbrInputStates; ++i)
	cout << "Charge imbalance between left and right (orbital cut) " << (TmpImbalanceEigenvalues[i]) << endl;
   
  if (Manager.GetBoolean("realspace-density") == true)
    {
      char* OutputFileName = 0;
      OutputFileName = ReplaceExtensionToFileName(InputStateNames[0], "vec", "rho.dat");
      ofstream File;
      File.precision(14);
      File.open(OutputFileName, ios::binary | ios::out);

      if (Perimeter == 0.0)
	{
	  Perimeter = sqrt(2.0 * M_PI * (LzMax + 1) * Ratio);
	}
      else
	{
	  Ratio = (Perimeter * Perimeter) / (2.0 * M_PI * (LzMax + 1));
	}
      double CylinderLength = Perimeter / Ratio;
      int NbrPoints = Manager.GetInteger("nbr-points");

      RealSymmetricMatrix TotalChargeImbalance(NbrInputStates, true);
      RealSymmetricMatrix TmpMatrixUpLayer (NbrInputStates, true);
      RealSymmetricMatrix TmpMatrixUpLayerIntegratedCharge (NbrInputStates, true);
      RealSymmetricMatrix TmpMatrixDownLayer (NbrInputStates, true);
      RealSymmetricMatrix TmpMatrixDownLayerIntegratedCharge (NbrInputStates, true);
      RealSymmetricMatrix TmpMatrixTotal (NbrInputStates, true);
      RealSymmetricMatrix TmpMatrixTotalIntegratedCharge (NbrInputStates, true);
      RealDiagonalMatrix TmpImbalanceEigenvaluesRealSpace(NbrInputStates);
      RealDiagonalMatrix TmpEigenvaluesUpLayer(NbrInputStates);
      RealDiagonalMatrix TmpEigenvaluesUpLayerIntegratedCharge(NbrInputStates);
      RealDiagonalMatrix TmpEigenvaluesDownLayer(NbrInputStates);
      RealDiagonalMatrix TmpEigenvaluesDownLayerIntegratedCharge(NbrInputStates);
      RealDiagonalMatrix TmpEigenvaluesTotal(NbrInputStates);
      RealDiagonalMatrix TmpEigenvaluesTotalIntegratedCharge(NbrInputStates);
      double Offset = Manager.GetDouble("offset");
      double XPosition = -(0.5 * CylinderLength) - Offset;
      double XStep = -2.0 * XPosition / ((double) (NbrPoints + 1));
      double TmpPrefactor1 =  1.0 / sqrt(M_PI);
      double TmpPrefactor2 =  2.0 * M_PI / Perimeter;
      double TmpShiftUp = (0.5 * ((double) LzMax) - Manager.GetDouble("flux-insertion")) * TmpPrefactor2;
      double TmpShiftDown = (0.5 * ((double) LzMax) - Manager.GetDouble("flux-insertion")) * TmpPrefactor2;
      for (int TmpX = 0; TmpX <= NbrPoints; ++TmpX)
	{ 
	  TmpMatrixUpLayer.ClearMatrix();
	  TmpMatrixUpLayerIntegratedCharge.ClearMatrix();
	  TmpMatrixDownLayer.ClearMatrix();
	  TmpMatrixDownLayerIntegratedCharge.ClearMatrix();
	  TmpMatrixTotal.ClearMatrix();
	  TmpMatrixTotalIntegratedCharge.ClearMatrix();
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      double TmpFactor = TmpPrefactor1 * exp(-(XPosition + TmpShiftUp - (TmpPrefactor2 * ((double) i))) * (XPosition + TmpShiftUp - (TmpPrefactor2 * ((double) i))));
	      TmpMatrixUpLayer.MultiplyAndAdd(TmpFactor, OneBodyMatrixElements[0][i]);
	      TmpMatrixTotal.MultiplyAndAdd(TmpFactor, OneBodyMatrixElements[0][i]);
	      TmpFactor = 0.5 * (1.0 + erf(XPosition + TmpShiftUp - (TmpPrefactor2 * ((double) i))));
	      TmpMatrixUpLayerIntegratedCharge.MultiplyAndAdd(TmpFactor, OneBodyMatrixElements[0][i]);
	      TmpMatrixTotalIntegratedCharge.MultiplyAndAdd(TmpFactor, OneBodyMatrixElements[0][i]);
	      TmpFactor = TmpPrefactor1 * exp(-(XPosition - TmpShiftDown + (TmpPrefactor2 * ((double) i))) * (XPosition - TmpShiftDown + (TmpPrefactor2 * ((double) i))));
	      TmpMatrixDownLayer.MultiplyAndAdd(TmpFactor, OneBodyMatrixElements[1][LzMax - i]);
	      TmpMatrixTotal.MultiplyAndAdd(TmpFactor, OneBodyMatrixElements[1][LzMax - i]);
	      TmpFactor = 0.5 * (1.0 + erf(XPosition - TmpShiftDown + (TmpPrefactor2 * ((double) i))));
	      TmpMatrixDownLayerIntegratedCharge.MultiplyAndAdd(TmpFactor, OneBodyMatrixElements[1][LzMax - i]);
	      TmpMatrixTotalIntegratedCharge.MultiplyAndAdd(TmpFactor, OneBodyMatrixElements[1][LzMax - i]);
	      
	    }
	  TmpMatrixUpLayer.LapackDiagonalize(TmpEigenvaluesUpLayer);
	  TmpMatrixUpLayerIntegratedCharge.LapackDiagonalize(TmpEigenvaluesUpLayerIntegratedCharge);
	  TmpMatrixDownLayer.LapackDiagonalize(TmpEigenvaluesDownLayer);
	  TmpMatrixDownLayerIntegratedCharge.LapackDiagonalize(TmpEigenvaluesDownLayerIntegratedCharge);
	  TmpMatrixTotal.LapackDiagonalize(TmpEigenvaluesTotal);
	  TmpMatrixTotalIntegratedCharge.LapackDiagonalize(TmpEigenvaluesTotalIntegratedCharge);
	  File << XPosition;
	  
// 	  if (TmpX <= NbrPoints / 2)
// 	    TotalChargeImbalance += TmpMatrixTotal;
// 	  else
// 	    TotalChargeImbalance -= TmpMatrixTotal;
	  for (int i = 0; i < NbrInputStates; ++i)
	    {   
	      File << " " << TmpEigenvaluesUpLayer[i] << " " << TmpEigenvaluesUpLayerIntegratedCharge[i]
		   << " " << TmpEigenvaluesDownLayer[i] << " " << TmpEigenvaluesDownLayerIntegratedCharge[i]
		   << " " << TmpEigenvaluesTotal[i] << " " << TmpEigenvaluesTotalIntegratedCharge[i];
		   
	      
	    }      
	  File << endl;
	  XPosition += XStep;
	}
//       TotalChargeImbalance.LapackDiagonalize(TmpImbalanceEigenvaluesRealSpace);
//       for (int i = 0; i < NbrInputStates; ++i)
// 	{
// 	  cout << "Charge imbalance between left and right (real space cut) " << (TmpImbalanceEigenvaluesRealSpace[i]) * XStep << endl;
// 	  if (OutputName != 0)
// 	    FileChargeImbalance << i << " " << TmpImbalanceEigenvalues[i] << " " << TmpImbalanceEigenvaluesRealSpace[i]* XStep << endl;
// 	}
      
      TotalChargeImbalance.ClearMatrix();
      for (int i = 0; i <= LzMax; ++i)
	{
	  double TmpFactor = -erf(+ TmpShiftUp - (TmpPrefactor2 * ((double) i)));
	  TotalChargeImbalance.MultiplyAndAdd(TmpFactor, OneBodyMatrixElements[0][i]);
	  TmpFactor = -erf(- TmpShiftDown + (TmpPrefactor2 * ((double) i)));
	  TotalChargeImbalance.MultiplyAndAdd(TmpFactor, OneBodyMatrixElements[1][LzMax - i]);	  
	}
      TotalChargeImbalance.LapackDiagonalize(TmpImbalanceEigenvaluesRealSpace);
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  cout << "Charge imbalance between left and right (real space cut) " << TmpImbalanceEigenvaluesRealSpace[i] << endl;
 	  if (OutputName != 0)
 	    FileChargeImbalance << i << " " << TmpImbalanceEigenvalues[i] << " " << TmpImbalanceEigenvaluesRealSpace[i] << endl;
	}
      File.close();
      

      if (OutputFileName != 0)
	FileChargeImbalance.close();
    }
  return 0;

}
