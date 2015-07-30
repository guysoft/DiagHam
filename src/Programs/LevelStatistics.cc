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


// parse the information contained in a single spectrum file
//
// spectrumFileName = specrtum file name
// manager = pointer to the option manager
// spectrum = reference on the two dimensional array where the spectrum will be stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minAverageSpacing = reference on the minimum level spacing 
// maxAverageSpacing = reference on the maximum level spacing 
// averageSpacing = reference on the average level spacing 
// nbrSpacings = reference on the number of level spacings
// return value = true if no error occured
bool LevelStatisticsParseSpectrumFile(char* spectrumFileName, OptionManager* manager, double**& spectrum, int*& spectrumSize, 
				      int*& spectrumWeight, int& nbrSectors, double& minAverageSpacing, double& maxAverageSpacing, double& averageSpacing,
				      long& nbrSpacings);


// perform level statistics on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minAverageSpacing = minimum level spacing 
// maxAverageSpacing = maximum level spacing 
// averageSpacing =average level spacing 
// nbrSpacings = number of level spacings
// nbrBins = number of bins for th level statistics
// binSize = level spacing range for each bin
// nbrSpacingPerBin = array that contains the number of level spacing per bin
// nbrRejectedSpacings = reference on the number of rejected level spacings (i.e. that cannot be stored in any bin)
// nbrAcceptedSpacings = reference on the number of accepted level spacings (i.e. that can be stored in a bin)
void LevelStatisticsPerformLevelStatistics(double** spectrum, int* spectrumSize, int* spectrumWeight, int nbrSectors, 
					   double minAverageSpacing, double maxAverageSpacing, double averageSpacing, long nbrSpacings,
					   int nbrBins, double binSize, long* nbrSpacingPerBin, long& nbrRejectedSpacings, long& nbrAcceptedSpacings);



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
  (*SystemGroup) += new SingleStringOption ('\n', "multiple-spectra", "column formated ASCII file that lists all the spectra to analyze");
  (*SystemGroup) += new SingleIntegerOption ('c', "energy-column", "index of the column that contains the energies (0 being the first column)", 1);
  (*SystemGroup) += new SingleIntegerOption ('d', "degeneracy-column", "index of the optional column that contains the energies degeneracy not appearing explicitly in the spectrum (must be larger than energy-column)", 0);
  (*SystemGroup) += new BooleanOption('\n', "z2-symmetry", "assume that the spectrum is invariant under the transformation Q<->-Q and that only the Q>=0 are available in the spectrum");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the output file where the level statistics will be stored (if none, try to deduce it from either --spectrum or --multiple-spectra replacing the dat extension with levelstat)");
  (*OutputGroup) += new SingleDoubleOption ('b', "bin-spacing", "spacing range for each bin", 0.1);
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

  if ((Manager.GetString("multiple-spectra") == 0) && (Manager.GetString("spectrum") == 0))
    {
      cout << "error, at least one spectrum should be provided" << endl;
      return 0;
    }

  char* OutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      if (Manager.GetString("multiple-spectra") == 0)
	{  
	  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectrum"), "dat", "levelstat");
	  if (OutputFileName == 0)
	    {
	      cout << "error, can't guess output file name from " << Manager.GetString("spectrum") << endl;
	      return 0;
	    }
	}
      else
	{
	  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("multiple-spectra"), "dat", "levelstat");
	  if (OutputFileName == 0)
	    {
	      cout << "error, can't guess output file name from " << Manager.GetString("multiple-spectra") << endl;
	      return 0;
	    }
	}
    }
  else
    {
      OutputFileName = new char [strlen (Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }


  double AverageSpacing = 0.0;
  double MinAverageSpacing = 1.0e300;
  double MaxAverageSpacing = 0.0;
  long NbrSpacings = 0l;
  int* SpectrumSize;
  double** Spectrum;
  int* SpectrumWeight;
  int NbrSectors;

  if (Manager.GetString("multiple-spectra") == 0)
    {
      if (LevelStatisticsParseSpectrumFile(Manager.GetString("spectrum"), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					   MinAverageSpacing, MaxAverageSpacing, AverageSpacing, NbrSpacings) == false)
	{
	  return 0;
	}
      for (int i = 0; i < NbrSectors; ++i)
	{    
	  if (SpectrumSize[i] > 0)
	    delete[] Spectrum[i];      
	}
      delete[] Spectrum;
      delete[] SpectrumSize;
      delete[] SpectrumWeight;
    }
  else
    {
      MultiColumnASCIIFile SpectraFile;
      if (SpectraFile.Parse(Manager.GetString("multiple-spectra")) == false)
	{
	  SpectraFile.DumpErrors(cout);
	  return false;
	}
      for (int i = 0; i < SpectraFile.GetNbrLines(); ++i)
	{
	  cout << "checking file " << SpectraFile(0, i) << endl;
	  if (LevelStatisticsParseSpectrumFile(SpectraFile(0, i), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					       MinAverageSpacing, MaxAverageSpacing, AverageSpacing, NbrSpacings) == false)
	    {
	      return 0;
	    }
	  for (int i = 0; i < NbrSectors; ++i)
	    {    
	      if (SpectrumSize[i] > 0)
		delete[] Spectrum[i];      
	    }
	  delete[] Spectrum;
	  delete[] SpectrumSize;
	  delete[] SpectrumWeight;	  
	}
    }
  
  AverageSpacing /= (double) NbrSpacings;
  MinAverageSpacing /= AverageSpacing;
  MaxAverageSpacing /= AverageSpacing;

  double BinSize = Manager.GetDouble("bin-spacing");
  int NbrBins = ((int) (MaxAverageSpacing / BinSize)) + 1;
  long* NbrSpacingPerBin = new long [NbrBins];
  for (int i = 0 ; i < NbrBins; ++i)
    NbrSpacingPerBin[i] = 0l;
  long NbrRejectedSpacings = 0l;
  long NbrAcceptedSpacings = 0l;


  if (Manager.GetString("multiple-spectra") == 0)
    {
      double DummyAverageSpacing = 0.0;
      double DummyMinAverageSpacing = 1.0e300;
      double DummyMaxAverageSpacing = 0.0;
      long DummyNbrSpacings = 0l; 
      if (LevelStatisticsParseSpectrumFile(Manager.GetString("spectrum"), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					   DummyMinAverageSpacing, DummyMaxAverageSpacing, DummyAverageSpacing, DummyNbrSpacings) == false)
	{
	  return 0;
	}

      LevelStatisticsPerformLevelStatistics(Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					    MinAverageSpacing, MaxAverageSpacing, AverageSpacing, NbrSpacings,
					    NbrBins, BinSize, NbrSpacingPerBin, NbrRejectedSpacings, NbrAcceptedSpacings);
      for (int i = 0; i < NbrSectors; ++i)
	{    
	  if (SpectrumSize[i] > 0)
	    delete[] Spectrum[i];      
	}
      delete[] Spectrum;
      delete[] SpectrumSize;
      delete[] SpectrumWeight;	  
    }
  else
    {
      MultiColumnASCIIFile SpectraFile;
      if (SpectraFile.Parse(Manager.GetString("multiple-spectra")) == false)
	{
	  SpectraFile.DumpErrors(cout);
	  return false;
	}
      for (int i = 0; i < SpectraFile.GetNbrLines(); ++i)
	{
	  double DummyAverageSpacing = 0.0;
	  double DummyMinAverageSpacing = 1.0e300;
	  double DummyMaxAverageSpacing = 0.0;
	  long DummyNbrSpacings = 0l; 
	  cout << "processing file " << SpectraFile(0, i) << endl;
	  if (LevelStatisticsParseSpectrumFile(SpectraFile(0, i), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					       DummyMinAverageSpacing, DummyMaxAverageSpacing, DummyAverageSpacing, DummyNbrSpacings) == false)
	    {
	      return 0;
	    }
	  
	  LevelStatisticsPerformLevelStatistics(Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
						MinAverageSpacing, MaxAverageSpacing, AverageSpacing, NbrSpacings,
						NbrBins, BinSize, NbrSpacingPerBin, NbrRejectedSpacings, NbrAcceptedSpacings);
	  for (int i = 0; i < NbrSectors; ++i)
	    {    
	      if (SpectrumSize[i] > 0)
		delete[] Spectrum[i];      
	    }
	  delete[] Spectrum;
	  delete[] SpectrumSize;
	  delete[] SpectrumWeight;	  
	}
    }

  double Sum = 0.0;
  for (int i = 0 ; i < NbrBins; ++i)
    {
      double PSpacing = (((double) NbrSpacingPerBin[i]) / ((double) NbrSpacings));
      Sum += PSpacing;
    }

  long MaxNbrSpacingPerBin = 0l;
  long MinNbrSpacingPerBin = NbrSpacings;
  for (int i = 0 ; i < NbrBins; ++i)
    {
      if (NbrSpacingPerBin[i] > MaxNbrSpacingPerBin)
	{
	  MaxNbrSpacingPerBin = NbrSpacingPerBin[i];
	}
      if (NbrSpacingPerBin[i] < MinNbrSpacingPerBin)
	{
	  MinNbrSpacingPerBin = NbrSpacingPerBin[i];
	}
    }

  ofstream File;
  File.open(OutputFileName, ios::binary | ios::out);
  File.precision(14);
  File << "# Min spacing " << MinAverageSpacing << endl;
  File << "# Max spacing " << MaxAverageSpacing << endl;
  File << "# Average spacing = " << AverageSpacing << endl;
  File << "# Nbr points = " << NbrSpacings << endl;
  File << "# Nbr rejected points = " << NbrRejectedSpacings << endl;
  File << "# int ds P(s) = " << Sum << endl;
  File << "# Min nbr spacings per bin " << MinNbrSpacingPerBin << endl;
  File << "# Max nbr spacings per bin " << MaxNbrSpacingPerBin << endl;
  File << "# Max spacing " << MaxAverageSpacing << endl;
  File << "# s P(s)" << endl;
  for (int i = 0 ; i < NbrBins; ++i)
    {
      double PSpacing = (((double) NbrSpacingPerBin[i]) / ((double) NbrSpacings));
      File << (BinSize * (((double) i) + 0.5)) << " " << (PSpacing / BinSize) << endl;
    }  
  File.close();

  cout << "Min spacing " << MinAverageSpacing << endl;
  cout << "Max spacing " << MaxAverageSpacing << endl;
  cout << "Average spacing = " << AverageSpacing << endl;
  cout << "Nbr points = " << NbrSpacings << endl;
  cout << "Nbr rejected points = " << NbrRejectedSpacings << endl;
  cout << "int ds P(s) = " << Sum << endl;
  cout << "Min nbr spacings per bin " << MinNbrSpacingPerBin << endl;
  cout << "Max nbr spacings per bin " << MaxNbrSpacingPerBin << endl;
  delete[] NbrSpacingPerBin;

  return 0;
}

// parse the information contained in a single spectrum file
//
// spectrumFileName = specrtum file name
// manager = pointer to the option manager
// spectrum = reference on the two dimensional array where the spectrum will be stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minAverageSpacing = reference on the minimum level spacing 
// maxAverageSpacing = reference on the maximum level spacing 
// averageSpacing = reference on the average level spacing 
// nbrSpacings = reference on the number of level spacings
// return value = true if no error occured

bool LevelStatisticsParseSpectrumFile(char* spectrumFileName, OptionManager* manager, double**& spectrum, int*& spectrumSize, 
				      int*& spectrumWeight, int& nbrSectors, double& minAverageSpacing, double& maxAverageSpacing, double& averageSpacing,
				      long& nbrSpacings)
{
  MultiColumnASCIIFile SpectrumFile;
  if (SpectrumFile.Parse(spectrumFileName) == false)
    {
      SpectrumFile.DumpErrors(cout);
      return false;
    }
  nbrSectors = 1;
  if (manager->GetInteger("energy-column") == 0)
    {
      int TmpSize = SpectrumFile.GetNbrLines();
      spectrumSize = new int [1];
      spectrum = new double*[1];
      spectrumWeight = new int[1];
      spectrum[0] = SpectrumFile.GetAsDoubleArray(0);
      spectrumSize[0] = TmpSize;
      spectrumWeight[0] = 1;
    }
  else
    {
      int TmpSize = SpectrumFile.GetNbrLines();
      int NbrQuantumNumber = manager->GetInteger("energy-column");
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
		  ++nbrSectors;
		}
	    }
	}
      spectrum = new double*[nbrSectors];
      spectrumSize = new int[nbrSectors];
      spectrumWeight = new int[nbrSectors];
      int* TmpDegeneracy = 0;
      if (manager->GetInteger("degeneracy-column") >= 0)
	{
	  TmpDegeneracy = SpectrumFile.GetAsIntegerArray(manager->GetInteger("degeneracy-column"));
	}
      nbrSectors = 0;
      int CurrentIndex = 0;
      spectrumWeight[0] = 1;
      if (manager->GetBoolean("z2-symmetry"))
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][0] > 0)
		{
		  spectrumWeight[0] *= 2;
		}
	    }
	}
      else
	{
	  if (TmpDegeneracy != 0)
	    {
	      spectrumWeight[0] = TmpDegeneracy[0];
	    }
	}
      for (int i = 1; i < TmpSize; ++i)
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][i - 1] != QuantumNumbers[j][i])
		{
		  spectrum[nbrSectors] = new double [i - CurrentIndex];
		  spectrumSize[nbrSectors] = i - CurrentIndex;
		  CurrentIndex = i;
		  j = NbrQuantumNumber;
		  ++nbrSectors;
		  spectrumWeight[nbrSectors] = 1;
		  if (manager->GetBoolean("z2-symmetry"))
		    {
		      for (int k = 0; k < NbrQuantumNumber; ++k)
			{
			  if (QuantumNumbers[k][i] > 0)
			    {
			      spectrumWeight[nbrSectors] *= 2;
			    }
			}
		    }
		  else
		    {
		      if (TmpDegeneracy != 0)
			{
			  spectrumWeight[nbrSectors] = TmpDegeneracy[i];
			}
		    }
		}
	    }
	}
      spectrum[nbrSectors] = new double [TmpSize - CurrentIndex];
      spectrumSize[nbrSectors]= TmpSize - CurrentIndex;
      ++nbrSectors;            
      double* TmpSpectrum = SpectrumFile.GetAsDoubleArray(manager->GetInteger("energy-column"));
      spectrum[0][0] = TmpSpectrum[0];
      nbrSectors = 0;
      CurrentIndex = 0;
      for (int i = 1; i < TmpSize; ++i)
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][i - 1] != QuantumNumbers[j][i])
		{
		  CurrentIndex = i;
		  j = NbrQuantumNumber;
		  ++nbrSectors;
		}
	    }
	  spectrum[nbrSectors][i - CurrentIndex] = TmpSpectrum[i];
	}
      ++nbrSectors;            
      for (int i = 0; i < NbrQuantumNumber; ++i)
	delete[] QuantumNumbers[i];
      delete[] QuantumNumbers;
    }

//   double Min = 0.0;
//   double Max = 0.0;
//   double TmpInfDiff;
//   double TmpSupDiff;
//   for (int i = 0; i < nbrSectors; ++i)
//     {     
//       if (spectrumSize[i] > 2)
// 	{
// 	  int Lim = spectrumSize[i] - 1;
// 	  for (int j = 1; j < Lim; ++j)
// 	    {
// 	      TmpInfDiff = spectrum[i][j] - spectrum[i][j - 1];
// 	      TmpSupDiff = spectrum[i][j + 1] - spectrum[i][j];
// 	      if (TmpInfDiff > TmpSupDiff)
// 		{
// 		  Min += TmpSupDiff * ((double) spectrumWeight[i]);
// 		  Max += TmpInfDiff * ((double) spectrumWeight[i]);
// 		}
// 	      else
// 		{
// 		  Max += TmpSupDiff * ((double) spectrumWeight[i]);
// 		  Min += TmpInfDiff * ((double) spectrumWeight[i]);
// 		}
// 	    }
// 	}
//     }

  for (int i = 0; i < nbrSectors; ++i)
    {     
      if (spectrumSize[i] > 1)
	{
	  nbrSpacings += spectrumSize[i] - 1;
	  int Lim = spectrumSize[i];	  
	  for (int j = 1; j < Lim; ++j)
	    {
	      double TmpDiff = spectrum[i][j] - spectrum[i][j - 1];
	      averageSpacing += TmpDiff;
	      if (TmpDiff > maxAverageSpacing)
		{
		  maxAverageSpacing = TmpDiff;
		}
	      if (TmpDiff < minAverageSpacing)
		{
		  minAverageSpacing = TmpDiff;
		}
	    }
	}
    }
  return true;
}

// perform level statistics on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minAverageSpacing = minimum level spacing 
// maxAverageSpacing = maximum level spacing 
// averageSpacing =average level spacing 
// nbrSpacings = number of level spacings
// nbrBins = number of bins for th level statistics
// binSize = level spacing range for each bin
// nbrSpacingPerBin = array that contains the number of level spacing per bin
// nbrRejectedSpacings = reference on the number of rejected level spacings (i.e. that cannot be stored in any bin)
// nbrAcceptedSpacings = reference on the number of accepted level spacings (i.e. that can be stored in a bin)

void LevelStatisticsPerformLevelStatistics(double** spectrum, int* spectrumSize, int* spectrumWeight, int nbrSectors, 
					   double minAverageSpacing, double maxAverageSpacing, double averageSpacing, long nbrSpacings,
					   int nbrBins, double binSize, long* nbrSpacingPerBin, long& nbrRejectedSpacings, long& nbrAcceptedSpacings)
{
  for (int i = 0; i < nbrSectors; ++i)
    {     
      if (spectrumSize[i] > 1)
	{
	  int Lim = spectrumSize[i];	  
	  for (int j = 1; j < Lim; ++j)
	    {
	      double TmpDiff = (spectrum[i][j] - spectrum[i][j - 1]) / averageSpacing;
	      int TmpIndex = int (TmpDiff / binSize);
	      if (TmpIndex < nbrBins)
		{
		  nbrSpacingPerBin[TmpIndex]++;
		  ++nbrAcceptedSpacings;
		}
	      else
		{
		  ++ nbrRejectedSpacings;
		}
	    }
	}
    }
}
