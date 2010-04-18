#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETorusEntanglementSpectrum" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "density-matrix", "file containing the reduced density matrix");  
  (*SystemGroup) += new SingleIntegerOption ('n', "nbr-particles", "number of particles in the A part", 0l);
  (*SystemGroup) += new SingleIntegerOption ('l', "nbr-orbitals", "number of orbitals in the A part", 0l);
  (*SystemGroup) += new SingleIntegerOption ('s', "sz-value", "twice the Sz value of A part (SU(2) mode only)", 0l);
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin (override autodetection from the reduced density matrix file name)");
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output name for the entanglement spectrum (default name replace density-matrix full.ent extension with la_x_na_y.entspec)");
  (*SystemGroup) += new SingleDoubleOption  ('e', "eigenvalue-error", "lowest acceptable reduced density matrix eignvalue", 1e-14);  
  (*SystemGroup) += new BooleanOption ('\n', "show-minmaxkya", "show minimum an maximum Ky value that can be reached");
  (*SystemGroup) += new BooleanOption ('\n', "show-counting", "show degeneracy counting for each Ky value");
  (*SystemGroup) += new BooleanOption ('\n', "particle-entanglement", "compute particle entanglement spectrum");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusEntanglementSpectrum -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrOrbitalsInPartition = Manager.GetInteger("nbr-orbitals");
  int NbrParticlesInPartition = Manager.GetInteger("nbr-particles");
  int TotalSzInPartition = Manager.GetInteger("sz-value");
  double Error = Manager.GetDouble("eigenvalue-error");
  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  int TotalKy = 0;
  int TotalSz = 0;
  bool Statistics = true;
  
  if (Manager.GetString("density-matrix") == 0)
    {
      cout << "a reduced density matrix has to be provided, see man page for option syntax or type FQHETorusEntanglementSpectrum -h" << endl;
      return -1;
    }

  bool SU2SpinFlag = false;
  if (strcasestr(Manager.GetString("density-matrix"), "_su2_") != 0)
    {
      SU2SpinFlag = true;
    }
  else
    {
      SU2SpinFlag = Manager.GetBoolean("su2-spin");      
    }

  if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("density-matrix"), NbrParticles, NbrFluxQuanta, TotalKy, Statistics) == false)
    {
      cout << "can't retrieve system informations from the reduced density matrix file name" << endl;
      return -1;
    } 

  MultiColumnASCIIFile DensityMatrix;
  if (DensityMatrix.Parse(Manager.GetString("density-matrix")) == false)
    {
      DensityMatrix.DumpErrors(cout);
      return -1;
    }

  int MinKya = 1 << 30;
  int MaxKya = -MinKya; 
  int* KyaValueArray = 0;


  if (Manager.GetBoolean("particle-entanglement") == false)
    {
      if (DensityMatrix.GetNbrColumns() < 4)
	{
	  cout << "wrong number of columns in " << Manager.GetString("density-matrix") << endl;
	  return -1;
	}
      
      int* LaValues = DensityMatrix.GetAsIntegerArray(0);
      int* NaValues = DensityMatrix.GetAsIntegerArray(1);
      int* KyValues = DensityMatrix.GetAsIntegerArray(2);
      long Index = 0l;
      long MaxIndex = DensityMatrix.GetNbrLines();
      while ((Index < MaxIndex) && (LaValues[Index] != NbrOrbitalsInPartition))
	++Index;
      
      if (Index < MaxIndex)
	{
	  if (SU2SpinFlag == false)
	    {
	      if (Index < MaxIndex)
		{
		  double* Coefficients = DensityMatrix.GetAsDoubleArray(3);
		  char* OutputFileName = Manager.GetString("output");
		  if (OutputFileName == 0)
		    {
		      char* TmpExtension = new char[256];
		      sprintf(TmpExtension, "la_%d_na_%d.entspec", NbrOrbitalsInPartition, NbrParticlesInPartition);
		      if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent", TmpExtension);
			}
		      else
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent.bz2", TmpExtension);
			}
		    }
		  ofstream File;
		  File.open(OutputFileName, ios::out);
		  File.precision(14);
		  File << "# la na ky shifted_ky lambda -log(lambda)" << endl;
		  if (NbrParticlesInPartition == 0)
		    {
		      int TmpIndex = Index;
		      while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition))
			{
			  double Tmp = Coefficients[Index];
			  if (Tmp > Error)
			    {
			      int TmpKya = (- (KyValues[Index] + ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NaValues[Index])));
			      File << NbrOrbitalsInPartition << " " << NaValues[Index] << " " << KyValues[Index] << " " <<  (0.5 * TmpKya) << " " << Tmp << " " << (-log(Tmp)) << endl;
			      if (TmpKya < MinKya)
				MinKya = TmpKya;
			      if (TmpKya > MaxKya)
				MaxKya = TmpKya;
			    }
			  ++Index;
			}
		      KyaValueArray = new int[(MaxKya - MinKya + 1) >> 1];
		      for (int i = MinKya; i <= MaxKya; i += 2)
			KyaValueArray[(i - MinKya) >> 1] = 0; 
		      Index = TmpIndex;
		      while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition))
			{
			  if (Coefficients[Index] > Error)
			    {
			      int TmpKya = (- (KyValues[Index] + ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NaValues[Index])));
			      KyaValueArray[(TmpKya - MinKya) >> 1]++; 
			    }
			  ++Index;
			}
		    }
		  else
		    {
		      while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
			++Index;
		      if (Index < MaxIndex)
			{
			  int TmpIndex = Index;
			  int Shift = ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NbrParticlesInPartition);
			  while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition) && (NaValues[Index] == NbrParticlesInPartition))
			    {
			      double Tmp = Coefficients[Index];
			      if (Tmp > Error)
				{
				  int TmpKya = (-(KyValues[Index] + Shift));
				  File << NbrOrbitalsInPartition << " " << NbrParticlesInPartition << " " << KyValues[Index] << " " << (0.5 * TmpKya) << " " << Tmp << " " << (-log(Tmp)) << endl;
				  if (TmpKya < MinKya)
				    MinKya = TmpKya;
				  if (TmpKya > MaxKya)
				    MaxKya = TmpKya;
				}
			      ++Index;
			    }	      
			  KyaValueArray = new int[(MaxKya - MinKya + 1) >> 1];
			  for (int i = MinKya; i <= MaxKya; i += 2)
			    KyaValueArray[(i - MinKya) >> 1] = 0; 
			  Index = TmpIndex;
			  while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition) && (NaValues[Index] == NbrParticlesInPartition))
			    {
			      if (Coefficients[Index] > Error)
				{
				  int TmpKya = (-(KyValues[Index] + Shift));
				  KyaValueArray[(TmpKya - MinKya) >> 1]++; 
				}
			      ++Index;
			    }
			}
		      else
			{
			  cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
			  return -1;
			}
		    }
		  File.close();
		}
	    }
	  else
	    {
	      int* SzaValues = DensityMatrix.GetAsIntegerArray(2);
	      while ((Index < MaxIndex) && (SzaValues[Index] != TotalSzInPartition))
		++Index;
	      if (Index < MaxIndex)
		{
		  double* Coefficients = DensityMatrix.GetAsDoubleArray(4);
		  char* OutputFileName = Manager.GetString("output");
		  if (OutputFileName == 0)
		    {
		      char* TmpExtension = new char[256];
		      sprintf(TmpExtension, "la_%d_na_%d_sza_%d_.entspec", NbrOrbitalsInPartition, NbrParticlesInPartition, TotalSzInPartition);
		      if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent", TmpExtension);
			}
		      else
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent.bz2", TmpExtension);
			}
		    }
		  ofstream File;
		  File.open(OutputFileName, ios::out);
		  File.precision(14);
		  File << "# la sza na ky shifted_ky lambda -log(lambda)" << endl;
		  if (NbrParticlesInPartition == 0)
		    {
		      while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition))
			{
			  double Tmp = Coefficients[Index];
			  if (Tmp > Error)
			    {
			      int TmpKya = (-(KyValues[Index] + ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NaValues[Index])));
			      File << NbrOrbitalsInPartition << " " << TotalSzInPartition << " " << NaValues[Index] << " " << KyValues[Index] << " " <<  (0.5 * TmpKya) << " " << Tmp << " " << (-log(Tmp)) << endl;
			      if (TmpKya < MinKya)
				MinKya = TmpKya;
			      if (TmpKya > MaxKya)
				MaxKya = TmpKya;
			    }
			  ++Index;
			}
		    }
		  else
		    {
		      while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
			++Index;
		      if (Index < MaxIndex)
			{
			  double Shift = ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NbrParticlesInPartition);
			  while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition) && (NaValues[Index] == NbrParticlesInPartition))
			    {
			      double Tmp = Coefficients[Index];
			      if (Tmp > Error)
				{
				  int TmpKya = (-(KyValues[Index] + Shift));
				  File << NbrOrbitalsInPartition << " " << TotalSzInPartition << " " << NbrParticlesInPartition << " " << KyValues[Index] << " " << (0.5 * TmpKya) << " " << Tmp << " " << (-log(Tmp)) << endl;
				  if (TmpKya < MinKya)
				    MinKya = TmpKya;
				  if (TmpKya > MaxKya)
				    MaxKya = TmpKya;
				}
			      ++Index;
			    }	      
			}
		      else
			{
			  cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
			  return -1;
			}
		    }
		  File.close();
		}
	      else
		{
		  cout << "error, no entanglement spectrum can be computed from current data (invalid Sz value)" << endl;
		  return -1;
		}
	    }
	}
      else
	{
	  cout << "error, no entanglement spectrum can be computed from current data (invalid number of orbitals)" << endl;
	  return -1;
	}
      
    }
  else
    {
      if (DensityMatrix.GetNbrColumns() < 3)
	{
	  cout << "wrong number of columns in " << Manager.GetString("density-matrix") << endl;
	  return -1;
	}
      
      int* NaValues = DensityMatrix.GetAsIntegerArray(0);
      int* KyValues = DensityMatrix.GetAsIntegerArray(1);
      long Index = 0l;
      long MaxIndex = DensityMatrix.GetNbrLines();
      while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
	++Index;

      if (Index < MaxIndex)
	{
	  double* Coefficients = DensityMatrix.GetAsDoubleArray(2);
	  char* OutputFileName = Manager.GetString("output");
	  if (OutputFileName == 0)
	    {
	      char* TmpExtension = new char[256];
	      sprintf(TmpExtension, "na_%d.parentspec", NbrParticlesInPartition);
	      if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
		{
		  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.parent", TmpExtension);
		}
	      else
		{
		  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.parent.bz2", TmpExtension);
		}
	    }
	  ofstream File;
	  File.open(OutputFileName, ios::out);
	  File.precision(14);
	  File << "# na ky lambda -log(lambda)";
	  File << endl;
	  int TmpIndex = Index;
	  while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
	    {
	      double Tmp = Coefficients[Index];
	      if (Tmp > Error)
		{
		  int TmpKya = KyValues[Index];
		  File << NbrParticlesInPartition << " " << TmpKya << " " << Tmp << " " << (-log(Tmp));
		  File << endl;
		  if (TmpKya < MinKya)
		    MinKya = TmpKya;
		  if (TmpKya > MaxKya)
		    MaxKya = TmpKya;
		}
	      ++Index;
	    }
	  KyaValueArray = new int[(MaxKya - MinKya + 1)];
	  for (int i = MinKya; i <= MaxKya; ++i)
	    KyaValueArray[(i - MinKya)] = 0; 
	  Index = TmpIndex;
	  while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
	    {
	      if (Coefficients[Index] > Error)
		{
		  KyaValueArray[(KyValues[Index] - MinKya)]++; 
		}
	      ++Index;
	    }
	  File.close();	      
	}
      else
	{
	  cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
	  return -1;
	}
    }

  if (Manager.GetBoolean("show-minmaxkya"))
    {
      cout << "min Kya = " << MinKya << endl;
      cout << "max Kya = " << MaxKya << endl;
    }
  
  if ((Manager.GetBoolean("show-counting")) && (KyaValueArray != 0))
    {
      cout << "degeneracy counting : " << endl;
      for (int i = MinKya; i <= MaxKya; ++i)
	{
	  cout << i << " " << KyaValueArray[(i - MinKya)] << endl; 
	}
    }

  return 0;
}

