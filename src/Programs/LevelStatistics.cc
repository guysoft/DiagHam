#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Options/Options.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("LevelStatistics" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('s', "spectrum", "name of the file that contains the spectrum");
  (*SystemGroup) += new SingleIntegerOption ('c', "energy-column", "index of the column that contains the energies (0 being the first column)", 1);
  (*SystemGroup) += new BooleanOption('\n', "z2-symmetry", "assume that the spectrum is invariant under the transformation Q<->-Q and that only the Q>=0 are available in the spectrum");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type LevelStatistics -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  MultiColumnASCIIFile SpectrumFile;
  if (SpectrumFile.Parse(Manager.GetString("spectrum")) == false)
    {
      SpectrumFile.DumpErrors(cout);
      return -1;
    }
  double** Spectrum;
  int* SpectrumSize;
  double* SpectrumWeight;
  int NbrSectors = 1;
  if (Manager.GetInteger("energy-column") == 0)
    {
      int TmpSize = SpectrumFile.GetNbrLines();
      SpectrumSize = new int [1];
      Spectrum = new double*[1];
      SpectrumWeight = new double[1];
      Spectrum[0] = SpectrumFile.GetAsDoubleArray(0);
      SpectrumSize[0] = TmpSize;
      SpectrumWeight[0] = 1.0;
    }
  else
    {
      int TmpSize = SpectrumFile.GetNbrLines();
      int NbrQuantumNumber = Manager.GetInteger("energy-column");
      int** QuantumNumbers =  new int*[NbrQuantumNumber];
      for (int i = 0; i < NbrQuantumNumber; ++i)
	QuantumNumbers[i] = SpectrumFile.GetAsIntegerArray(i);
      for (int i = 1; i < TmpSize; ++i)
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][i - 1] != QuantumNumbers[j][i])
		{
		  j = NbrQuantumNumber;
		  ++NbrSectors;
		}
	    }
	}
      Spectrum = new double*[NbrSectors];
      SpectrumSize = new int[NbrSectors];
      SpectrumWeight = new double[NbrSectors];
      NbrSectors = 0;
      int CurrentIndex = 0;
      SpectrumWeight[0] = 1.0;
      if (Manager.GetBoolean("z2-symmetry"))
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][0] > 0)
		{
		  SpectrumWeight[0] *= 2.0;
		}
	    }
	}
      for (int i = 1; i < TmpSize; ++i)
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][i - 1] != QuantumNumbers[j][i])
		{
		  Spectrum[NbrSectors] = new double [i - CurrentIndex];
		  SpectrumSize[NbrSectors] = i - CurrentIndex;
		  CurrentIndex = i;
		  j = NbrQuantumNumber;
		  ++NbrSectors;
		  SpectrumWeight[NbrSectors] = 1.0;
		  if (Manager.GetBoolean("z2-symmetry"))
		    {
		      for (int k = 0; k < NbrQuantumNumber; ++k)
			{
			  if (QuantumNumbers[k][i] > 0)
			    {
			      SpectrumWeight[NbrSectors] *= 2.0;
			    }
			}
		    }
		}
	    }
	}
      Spectrum[NbrSectors] = new double [TmpSize - CurrentIndex];
      SpectrumSize[NbrSectors]= TmpSize - CurrentIndex;
      ++NbrSectors;            
      double* TmpSpectrum = SpectrumFile.GetAsDoubleArray(Manager.GetInteger("energy-column"));
      Spectrum[0][0] = TmpSpectrum[0];
      NbrSectors = 0;
      CurrentIndex = 0;
      for (int i = 1; i < TmpSize; ++i)
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][i - 1] != QuantumNumbers[j][i])
		{
		  CurrentIndex = i;
		  j = NbrQuantumNumber;
		  ++NbrSectors;
		}
	    }
	  Spectrum[NbrSectors][i - CurrentIndex] = TmpSpectrum[i];
	}
      ++NbrSectors;            
      for (int i = 0; i < NbrQuantumNumber; ++i)
	delete[] QuantumNumbers[i];
      delete[] QuantumNumbers;
    }

  double Min = 0.0;
  double Max = 0.0;
  double TmpInfDiff;
  double TmpSupDiff;
  for (int i = 0; i < NbrSectors; ++i)
    {     
      if (SpectrumSize[i] > 2)
	{
	  int Lim = SpectrumSize[i] - 1;
	  for (int j = 1; j < Lim; ++j)
	    {
	      TmpInfDiff = Spectrum[i][j] - Spectrum[i][j - 1];
	      TmpSupDiff = Spectrum[i][j + 1] - Spectrum[i][j];
	      if (TmpInfDiff > TmpSupDiff)
		{
		  Min += TmpSupDiff * SpectrumWeight[i];
		  Max += TmpInfDiff * SpectrumWeight[i];
		}
	      else
		{
		  Max += TmpSupDiff * SpectrumWeight[i];
		  Min += TmpInfDiff * SpectrumWeight[i];
		}
	    }
	}
//       for (int j = 0; j < SpectrumSize[i]; ++j)
// 	cout << i << " " << j << " " << Spectrum[i][j] << endl;
    }
  cout << (Min / Max) << endl;
  return 0;
}
