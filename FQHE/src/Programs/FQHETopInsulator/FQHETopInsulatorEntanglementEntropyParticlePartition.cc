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

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

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
  OptionManager Manager ("FQHETopInsulatorEntanglementEntropyParticlePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
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
  FermionOnSquareLatticeMomentumSpace** Spaces = 0;
  ComplexVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  int* TotalKx = 0;
  int* TotalKy = 0;
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  bool Statistics = true;

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
      TotalKy = new int[1];
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
       TotalKy = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalKx[i] = 0;
      TotalKy[i] = 0;
      double Mass = 0.0;
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i],
							      NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], Mass, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
    }


  GroundStates = new ComplexVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[i] << endl;
	return -1;      
      }

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  N    Kx    Ky    lambda";
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

  int MaxSubsystemNbrParticles = (NbrParticles >> 1) + (NbrParticles & 1);
  if (Manager.GetInteger("max-na") > 0)
    MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
  int SubsystemNbrParticles = Manager.GetInteger("min-na");

  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
    }
  File.close();
  return 0;
}
