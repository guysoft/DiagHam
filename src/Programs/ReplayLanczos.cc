#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithGroundStateDiskStorage.h"
#include "Vector/RealVector.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureManager.h"

#include "Options/Options.h"

#include "GeneralTools/Endian.h"


#include <iostream>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("ReplayLanczos" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += LanczosGroup;
  ArchitectureManager Architecture;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigenvalue", "number of eigenvalues to be calculated", 1);
  (*LanczosGroup) += new BooleanOption  ('e', "eigenstate", "compute the ground state", false);
  (*LanczosGroup) += new  SingleStringOption ('o', "ground-filename", "name of the file where the ground state has to be stored (in default ground.vec.XX, with XX=nbr iteration)", 0);
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-iter", "set a new number of lanczos iteration (0 if the one of the lanczos.dat has to be kept)", 0);
  (*LanczosGroup) += new BooleanOption ('r', "reorthogonalized", "indicate whether a reorthogonalized Lanczos algorithm was used");
  (*LanczosGroup) += new BooleanOption ('c', "complex-lanczos", "indicate whether a complex Lanczos algorithm was used");
  (*LanczosGroup) += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type ReplayFastLanczos -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  bool EigenstateFlag = Manager.GetBoolean("eigenstate");
  bool BlockLanczosFlag = ((BooleanOption*) Manager["block-lanczos"])->GetBoolean();
  int NbrEigenvalue = Manager.GetInteger("nbr-eigenvalue");

  char *OutputName;

  if (Manager.GetString("ground-filename")!=NULL)
    {
      OutputName = new char[strlen(Manager.GetString("ground-filename"))+16];
    }
  else
    {
      ifstream File;
      int LanczosIndex;
      File.open("lanczos.dat", ios::binary | ios::in);
      ReadLittleEndian(File, LanczosIndex);
      File.close();
      OutputName = new char[30];
      sprintf(OutputName, "ground.%d", LanczosIndex);
    }
  
  if (!Manager.GetBoolean("complex-lanczos"))
    {
      if (BlockLanczosFlag == false)
	{
	  if (Manager.GetBoolean("reorthogonalized"))
	    {
	      FullReorthogonalizedLanczosAlgorithmWithDiskStorage Lanczos(Architecture.GetArchitecture(), NbrEigenvalue, 0,0);
	      Lanczos.ResumeLanczosAlgorithm();
	      if (EigenstateFlag)
		{
		  char *TmpVectorName = new char[strlen(OutputName)+10];
		  RealVector* Eigenvectors = (RealVector*) Lanczos.GetEigenstates(NbrEigenvalue);
		  for (int i = 0; i < NbrEigenvalue; ++i)
		    {
		      sprintf (TmpVectorName, "%s.%d.vec", OutputName, i);
		      Eigenvectors[i].WriteVector(TmpVectorName);
		    }
		  delete [] TmpVectorName;		  
		}
	    }
	  else
	    {
	      BasicLanczosAlgorithmWithGroundStateDiskStorage Lanczos(Architecture.GetArchitecture(), 0,0);
	      Lanczos.ResumeLanczosAlgorithm();
	      if (EigenstateFlag)
		{
		  char *TmpVectorName = new char[strlen(OutputName)+10];
		  RealVector* Eigenvectors = (RealVector*) Lanczos.GetEigenstates(NbrEigenvalue);
		  for (int i = 0; i < NbrEigenvalue; ++i)
		    {
		      sprintf (TmpVectorName, "%s.%d.vec", OutputName, i);
		      Eigenvectors[i].WriteVector(TmpVectorName);
		    }
		  delete [] TmpVectorName;		  
		}
	    }
	}
      else
	{
	  cout << "Block Lanczos not defined yet"<<endl;
	  exit(-1);
	}
    }
  else // have complex Lanczos data
    {
      if (BlockLanczosFlag == false)
	{
	  FullReorthogonalizedLanczosAlgorithmWithDiskStorage Lanczos(Architecture.GetArchitecture(), NbrEigenvalue, 0,0);
	}
      else
	{
	  cout << "Block Lanczos not defined yet"<<endl;
	  exit(-1);
	}
    }
  delete [] OutputName;
  return 0;  
}
    
