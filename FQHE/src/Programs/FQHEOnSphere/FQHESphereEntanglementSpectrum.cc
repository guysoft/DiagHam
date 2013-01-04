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

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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
  OptionManager Manager ("FQHESphereEntanglementSpectrum" , "0.01");
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
  (*SystemGroup) += new BooleanOption ('\n', "show-minmaxlza", "show minimum an maximum Lz value that can be reached");
  (*SystemGroup) += new BooleanOption ('\n', "show-counting", "show degeneracy counting for each Lz value");
  (*SystemGroup) += new BooleanOption ('\n', "particle-entanglement", "compute particle entanglement spectrum");
  (*SystemGroup) += new BooleanOption ('\n', "ls-sorted", "for the particle entanglement spectrum with su2-spin on the sphere, sort the spectrum with respect to L and S");
  (*SystemGroup) += new SingleDoubleOption ('\n', "degeneracy-error", "error below which two reduced density matrix eigenvalues are assumed to be degenerated", 1e-12);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereEntanglementSpectrum -h" << endl;
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
  int TotalLz = 0;
  int TotalSz = 0;
  bool Statistics = true;
  
  if (Manager.GetString("density-matrix") == 0)
    {
      cout << "a reduced density matrix has to be provided, see man page for option syntax or type FQHESphereEntanglementSpectrum -h" << endl;
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

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("density-matrix"), NbrParticles, NbrFluxQuanta, TotalLz, Statistics) == false)
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

  int MinLza = 1 << 30;
  int MaxLza = -MinLza; 
  int* LzaValueArray = 0;


  if (Manager.GetBoolean("particle-entanglement") == false)
    {
      if (DensityMatrix.GetNbrColumns() < 4)
	{
	  cout << "wrong number of columns in " << Manager.GetString("density-matrix") << endl;
	  return -1;
	}
      
      int* LaValues = DensityMatrix.GetAsIntegerArray(0);
      int* NaValues = DensityMatrix.GetAsIntegerArray(1);
      int* LzValues = DensityMatrix.GetAsIntegerArray(2);
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
		  File << "# la na lz shifted_lz lambda -log(lambda)" << endl;
		  if (NbrParticlesInPartition == 0)
		    {
		      int TmpIndex = Index;
		      while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition))
			{
			  double Tmp = Coefficients[Index];
			  if (Tmp > Error)
			    {
			      int TmpLza = (- (LzValues[Index] + ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NaValues[Index])));
			      File << NbrOrbitalsInPartition << " " << NaValues[Index] << " " << LzValues[Index] << " " <<  (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp)) << endl;
			      if (TmpLza < MinLza)
				MinLza = TmpLza;
			      if (TmpLza > MaxLza)
				MaxLza = TmpLza;
			    }
			  ++Index;
			}
		      LzaValueArray = new int[((MaxLza - MinLza) >> 1) + 1];
		      for (int i = MinLza; i <= MaxLza; i += 2)
			LzaValueArray[(i - MinLza) >> 1] = 0; 
		      Index = TmpIndex;
		      while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition))
			{
			  if (Coefficients[Index] > Error)
			    {
			      int TmpLza = (- (LzValues[Index] + ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NaValues[Index])));
			      LzaValueArray[(TmpLza - MinLza) >> 1]++; 
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
				  int TmpLza = (-(LzValues[Index] + Shift));
				  File << NbrOrbitalsInPartition << " " << NbrParticlesInPartition << " " << LzValues[Index] << " " << (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp)) << endl;
				  if (TmpLza < MinLza)
				    MinLza = TmpLza;
				  if (TmpLza > MaxLza)
				    MaxLza = TmpLza;
				}
			      ++Index;
			    }	      
			  LzaValueArray = new int[((MaxLza - MinLza) >> 1) + 1];
			  for (int i = MinLza; i <= MaxLza; i += 2)
			    LzaValueArray[(i - MinLza) >> 1] = 0; 
			  Index = TmpIndex;
			  while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition) && (NaValues[Index] == NbrParticlesInPartition))
			    {
			      if (Coefficients[Index] > Error)
				{
				  int TmpLza = (-(LzValues[Index] + Shift));
				  LzaValueArray[(TmpLza - MinLza) >> 1]++; 
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
	      if (((TotalSzInPartition ^ NbrParticlesInPartition) & 1) == 0)
		{
		  int* SzaValues = DensityMatrix.GetAsIntegerArray(3);
		  double* Coefficients = DensityMatrix.GetAsDoubleArray(4);
		  char* OutputFileName = Manager.GetString("output");
		  if (OutputFileName == 0)
		    {
		      char* TmpExtension = new char[256];
		      sprintf(TmpExtension, "la_%d_na_%d_sza_%d.entspec", NbrOrbitalsInPartition, NbrParticlesInPartition, TotalSzInPartition);
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
		  File << "# la sza na lz shifted_lz lambda -log(lambda)" << endl;
		  
		  while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
		    ++Index;
		  if (Index < MaxIndex)
		    {
		      while ((Index < MaxIndex) && (SzaValues[Index] != TotalSzInPartition))
			++Index;
		      if (Index < MaxIndex)
			{
			  double Shift = ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NbrParticlesInPartition);
			  while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition) && (NaValues[Index] == NbrParticlesInPartition) && (SzaValues[Index] == TotalSzInPartition))
			    {
			      double Tmp = Coefficients[Index];
			      if (Tmp > Error)
				{
				  int TmpLza = (-(LzValues[Index] + Shift));
				  File << LaValues[Index] << " " << SzaValues[Index] << " " << NaValues[Index] << " " << LzValues[Index] << " " <<  (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp)) << endl;
				  if (TmpLza < MinLza)
				    MinLza = TmpLza;
				  if (TmpLza > MaxLza)
				    MaxLza = TmpLza;
				}
			      ++Index;
			    }	      
			}
		      else
			{
			  cout << "error, no entanglement spectrum can be computed from current data (invalid Sz value)" << endl;
			  return -1;
			}
		    }
		  else
		    {
		      cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
		      return -1;
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

      if (SU2SpinFlag == false)
	{
	  int* LzValues = DensityMatrix.GetAsIntegerArray(1);
	  double* LValues = 0;
	  double* L2Values = 0;
	  if (DensityMatrix.GetNbrColumns() > 3)
	    {
	      L2Values = DensityMatrix.GetAsDoubleArray(3);
	      LValues = DensityMatrix.GetAsDoubleArray(4);
	    }
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
	      File << "# na lz lambda -log(lambda)";
	      if (LValues != 0)
		{
		  File << " L^2 L";
		}
	      File << endl;
	      int TmpIndex = Index;
	      while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
		{
		  double Tmp = Coefficients[Index];
		  if (Tmp > Error)
		    {
		      int TmpLza = LzValues[Index];
		      File << NbrParticlesInPartition << " " << (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp));
		      if (LValues != 0)
			{
			  File << " " << L2Values[Index] << " " << LValues[Index];
			}
		      File << endl;
		      if (TmpLza < MinLza)
			MinLza = TmpLza;
		      if (TmpLza > MaxLza)
			MaxLza = TmpLza;
		    }
		  ++Index;
		}
	      LzaValueArray = new int[((MaxLza - MinLza ) >> 1) + 1];
	      for (int i = MinLza; i <= MaxLza; i += 2)
		LzaValueArray[(i - MinLza) >> 1] = 0; 
	      Index = TmpIndex;
	      while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
		{
		  if (Coefficients[Index] > Error)
		    {
		      LzaValueArray[(LzValues[Index] - MinLza) >> 1]++; 
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
      else
	{
	  int* LzValues = DensityMatrix.GetAsIntegerArray(4);
	  int* SzValues = DensityMatrix.GetAsIntegerArray(1);
	  int** SzaLzaValueArray = 0;
	  long Index = 0l;
	  long MaxIndex = DensityMatrix.GetNbrLines();
	  while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
	    ++Index;
	  
	  if (Index < MaxIndex)
	    {
	      double* Coefficients = DensityMatrix.GetAsDoubleArray(5);
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
	      if (Manager.GetBoolean("ls-sorted") == false)
		{
		  File << "# na sz lz lambda -log(lambda)";
		  File << endl;
		  int TmpIndex = Index;
		  while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
		    {
		      double Tmp = Coefficients[Index];
		      if (Tmp > Error)
			{
			  int TmpLza = LzValues[Index];
			  File << NbrParticlesInPartition << " " << (0.5 * SzValues[Index]) << " " << (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp));
			  File << endl;
			  if (TmpLza < MinLza)
			    MinLza = TmpLza;
			  if (TmpLza > MaxLza)
			    MaxLza = TmpLza;
			}
		      ++Index;
		    }
		  SzaLzaValueArray = new int* [NbrParticlesInPartition + 1];
		  for (int TmpSza = 0; TmpSza <= NbrParticlesInPartition; ++TmpSza)
		    {
		      SzaLzaValueArray[TmpSza] = new int[((MaxLza - MinLza ) >> 1) + 1];
		      for (int i = MinLza; i <= MaxLza; i += 2)
			SzaLzaValueArray[TmpSza][(i - MinLza) >> 1] = 0; 
		    }
		  Index = TmpIndex;
		  while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
		    {
		      if (Coefficients[Index] > Error)
			{
			  SzaLzaValueArray[(SzValues[Index] + NbrParticlesInPartition) >> 1][(LzValues[Index] - MinLza) >> 1]++; 
			}
		      ++Index;
		    }
		  if (Manager.GetBoolean("show-minmaxlza"))
		    {
		      cout << "min Lza = " << MinLza << endl;
		      cout << "max Lza = " << MaxLza << endl;
		    }
		  
		  if ((Manager.GetBoolean("show-counting")) && (SzaLzaValueArray != 0))
		    {
		      cout << "degeneracy counting (SzA Lza NbrStates) : " << endl;
		      for (int j = 0; j <= NbrParticlesInPartition; ++j)
			{
			  for (int i = MinLza; i <= MaxLza; i += 2)
			    {
			      cout << (j - (NbrParticlesInPartition * 0.5)) << " " << (0.5 * i) << " " << SzaLzaValueArray[j][(i - MinLza) >> 1] << endl; 
			    }
			}
		    }	      
		}
	      else
		{
		  File << "# na s l lambda -log(lambda)";
		  File << endl;
		  int MinSza = 1 << 30;
		  int MaxSza = -MinSza; 
		  int TmpIndex = Index;
		  int LargestSzSector = SzValues[Index];
		  int LargestLzSector = LzValues[Index];
		  while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
		    {
		      if (abs(LzValues[Index]) < abs(LargestLzSector))
			LargestLzSector = LzValues[Index];
		      if (abs(SzValues[Index]) < abs(LargestSzSector))
			LargestSzSector = SzValues[Index];
		      double Tmp = Coefficients[Index];
		      if (Tmp > Error)
			{
			  int TmpLza = LzValues[Index];
			  if ((TmpLza < MinLza) && (TmpLza >= 0))
			    MinLza = TmpLza;
			  if ((TmpLza > MaxLza) && (TmpLza >= 0))
			    MaxLza = TmpLza;
			  int TmpSza = SzValues[Index];
			  if ((TmpSza < MinSza) && (TmpSza >= 0))
			    MinSza = TmpSza;
			  if ((TmpSza > MaxSza) && (TmpSza >= 0))
			    MaxSza = TmpSza;
			}
		      ++Index;
		    }
		  MaxIndex = Index;
		  Index = TmpIndex;
		  while ((Index < MaxIndex) && ((LzValues[Index] != LargestLzSector) || (SzValues[Index] != LargestSzSector)))
		    {
		      ++Index;
		    }
		  int LargestSectorStartingIndex = Index;
		  while ((Index < MaxIndex) && (LzValues[Index] == LargestLzSector) && (SzValues[Index] == LargestSzSector))
		    {
		      ++Index;
		    }
		  int LargestSectorEndingIndex = Index;
		  int LargestSectorNbrIndices = LargestSectorEndingIndex - LargestSectorStartingIndex;
		  int* LSector = new int [LargestSectorNbrIndices];
		  int* SSector = new int [LargestSectorNbrIndices];
		  double* ErrorSector = new double [LargestSectorNbrIndices];
		  LargestLzSector = abs(LargestLzSector);
		  LargestSzSector = abs(LargestSzSector);
		  double DegeneracyError = Manager.GetDouble("degeneracy-error");
		  for (int i = 0; i < LargestSectorNbrIndices; ++i)
		    {
		      double Tmp = Coefficients[i + LargestSectorStartingIndex];
		      if (Tmp > Error)
			{
			  int TmpLSector = LargestLzSector;
			  int TmpSSector = LargestSzSector;
			  Index = TmpIndex;
			  while (Index < MaxIndex)
			    {
			      double Tmp2 = Coefficients[Index];
			      
			      if ((Tmp2 > Error)  && (fabs(Tmp2 - Tmp) < DegeneracyError))
				{
				  if (abs(LzValues[Index]) > TmpLSector)
				    TmpLSector = abs(LzValues[Index]);
				  if (abs(SzValues[Index]) > TmpSSector)
				    TmpSSector = abs(SzValues[Index]);			      
				}
			      ++Index;
			    }
			  LSector[i] = TmpLSector;
			  SSector[i] = TmpSSector;		      
			}
		    }
		  SzaLzaValueArray = new int* [NbrParticlesInPartition + 1];
		  for (int TmpSza = 0; TmpSza <= NbrParticlesInPartition; ++TmpSza)
		    {
		      SzaLzaValueArray[TmpSza] = new int[((MaxLza - MinLza ) >> 1) + 1];
		      for (int i = MinLza; i <= MaxLza; i += 2)
			SzaLzaValueArray[TmpSza][(i - MinLza) >> 1] = 0; 
		    }
		  for (int i = 0; i < LargestSectorNbrIndices; ++i)
		    {
		      double Tmp = Coefficients[i + LargestSectorStartingIndex];
		      if (Tmp > Error)
			{
			  SzaLzaValueArray[SSector[i]][(LSector[i] - MinLza) >> 1]++;
			  File << NbrParticlesInPartition << " " << (0.5 * SSector[i]) << " " << (0.5 * LSector[i]) << " " << Tmp << " " << (-log(Tmp));
			  File << endl;
			}
		    }
		}
	      if (Manager.GetBoolean("show-minmaxlza"))
		{
		  cout << "min Lza = " << MinLza << endl;
		  cout << "max Lza = " << MaxLza << endl;
		}
	      
	      if ((Manager.GetBoolean("show-counting")) && (SzaLzaValueArray != 0))
		{
		  cout << "degeneracy counting (Sa La NbrStates) : " << endl;
		  for (int j = (NbrParticlesInPartition & 1); j <= NbrParticlesInPartition; j += 2)
		    {
		      for (int i = MinLza; i <= MaxLza; i += 2)
			{
			  if (SzaLzaValueArray[j][(i - MinLza) >> 1] > 0)
			    cout << (j * 0.5) << " " << (0.5 * i) << " " << SzaLzaValueArray[j][(i - MinLza) >> 1] << endl; 
			}
		    }
		}
	      File.close();	      

	      return 0;	      
	    }
	  else
	    {
	      cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
	      return -1;
	    }
	}
    }

  if (Manager.GetBoolean("show-minmaxlza"))
    {
      cout << "min Lza = " << MinLza << endl;
      cout << "max Lza = " << MaxLza << endl;
    }
  
  if ((Manager.GetBoolean("show-counting")) && (LzaValueArray != 0))
    {
      cout << "degeneracy counting : " << endl;
      for (int i = MinLza; i <= MaxLza; i += 2)
	{
	  cout << (0.5 * i) << " " << LzaValueArray[(i - MinLza) >> 1] << endl; 
	}
    }

  return 0;
}

