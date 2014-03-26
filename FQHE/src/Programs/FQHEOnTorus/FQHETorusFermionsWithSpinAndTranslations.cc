#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslations.h"
#include "Hamiltonian/ParticleOnTorusCoulombHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/Options.h"

#include "MainTask/FQHEOnTorusMainTask.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);
    
  // some running options and help
  OptionManager Manager ("FQHETorusFermionsWithSpinAndTranslations" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "max-momentum", "maximum momentum for a single particle", 18);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-spin", "total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "x-momentum", "constraint on the total momentum in the x direction (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "constraint on the total momentum in the y direction (negative if none)", -1);
  (*SystemGroup) += new SingleDoubleOption   ('r', "ratio", 
					      "ratio between lengths along the x and y directions (-1 if has to be taken equal to nbr-particles/4)", 1.0);
  (*SystemGroup) += new SingleDoubleOption   ('d', "layerSeparation", 
					      "for bilayer simulations: layer separation in magnetic lengths", 0.0);
  (*SystemGroup) += new  BooleanOption  ('\n', "redundantYMomenta", "Calculate all subspaces up to YMomentum = MaxMomentum-1", false);

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusFermionsWithSpinAndTranslations -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int TotalSpin = Manager.GetInteger("total-spin");
  int NbrIterLanczos = Manager.GetInteger("nbr-iter");
  int NbrFermions = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  int XMomentum = Manager.GetInteger("x-momentum");
  int YMomentum = Manager.GetInteger("y-momentum");
  double XRatio = NbrFermions / 4.0;
  if (Manager.GetDouble("ratio") > 0)
    {
       XRatio = Manager.GetDouble("ratio");
    }
  double LayerSeparation = Manager.GetDouble("layerSeparation");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  int L = 0;
  double GroundStateEnergy = 0.0;

  char* OutputName = new char [512];
  if (LayerSeparation==0.0)
    sprintf (OutputName, "fermions_torus_su2_coulomb_n_%d_2s_%d_sz_%d_ratio_%f.dat", NbrFermions, MaxMomentum, TotalSpin, XRatio);
  else
    sprintf (OutputName, "fermions_torus_d_%f_coulomb_n_%d_2s_%d_sz_%d_ratio_%f.dat", LayerSeparation, 
NbrFermions,MaxMomentum, TotalSpin, XRatio);
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);


  int MomentumModulo = FindGCD(NbrFermions, MaxMomentum);
  int XMaxMomentum = (MomentumModulo - 1);
  if (XMomentum < 0)
    XMomentum = 0;
  else
    XMaxMomentum = XMomentum;
  int YMaxMomentum;
  if (Manager.GetBoolean("redundantYMomenta"))
    YMaxMomentum = (MaxMomentum - 1);
  else
    YMaxMomentum = (MomentumModulo - 1);
  if (YMomentum < 0)
    YMomentum = 0;
  else
    YMaxMomentum = YMomentum; 

  bool FirstRun=true;
  for (; XMomentum <= XMaxMomentum; ++XMomentum)
    for (int YMomentum2 = YMomentum; YMomentum2 <= YMaxMomentum; ++YMomentum2)
      {     
	cout << "----------------------------------------------------------------" << endl;
	cout << " Ratio = " << XRatio << endl;
	FermionOnTorusWithSpinAndMagneticTranslations* TotalSpace = 0 ;
	TotalSpace = new FermionOnTorusWithSpinAndMagneticTranslations (NbrFermions, TotalSpin, MaxMomentum, XMomentum, YMomentum2);	
	cout << " Total Hilbert space dimension = " << TotalSpace->GetHilbertSpaceDimension() << endl;
	cout << "momentum = (" << XMomentum << "," << YMomentum2 << ")" << endl;
	Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());	
	AbstractQHEHamiltonian* Hamiltonian = new ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian
	  (TotalSpace, NbrFermions, MaxMomentum, XMomentum, XRatio, LayerSeparation, 
	   Architecture.GetArchitecture(), Memory);
	char* EigenvectorName = 0;
	if (Manager.GetBoolean("eigenstate"))	
	  {
	    EigenvectorName = new char [512];
	    char* TmpName = RemoveExtensionFromFileName(OutputName, ".dat");
	    if (Manager.GetString("eigenstate-file") == 0)
              sprintf (EigenvectorName, "%s_kx_%d_ky_%d", TmpName, XMomentum, YMomentum);
	    else
              sprintf (EigenvectorName, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), XMomentum, YMomentum);
	    delete [] TmpName;
	  }
	double Shift = -10.0;
	Hamiltonian->ShiftHamiltonian(Shift);      
	FQHEOnTorusMainTask Task (&Manager, TotalSpace, &Lanczos, Hamiltonian, YMomentum2, Shift, OutputName, FirstRun, EigenvectorName);
	Task.SetKxValue(XMomentum);
	MainTaskOperation TaskOperation (&Task);
	TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	if (EigenvectorName != 0)
	  {
	    delete[] EigenvectorName;
	  }
	if (FirstRun == true)
	  FirstRun = false;
	delete Hamiltonian;
	delete TotalSpace;
      }
  File.close();
  delete[] OutputName;
  return 0;
}
