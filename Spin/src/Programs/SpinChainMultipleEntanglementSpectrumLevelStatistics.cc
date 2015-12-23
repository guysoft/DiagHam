#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/Endian.h"

#include "MathTools/NumericalAnalysis/Constant1DRealFunction.h"
#include "MathTools/NumericalAnalysis/Tabulated1DRealFunction.h"
#include "MathTools/NumericalAnalysis/Linear1DRealFunction.h"

#include "MathTools/NumericalAnalysis/Constant1DRealFunction.h"
#include "MathTools/NumericalAnalysis/Tabulated1DRealFunction.h"
#include "MathTools/NumericalAnalysis/Linear1DRealFunction.h"

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
using std::ifstream;
using std::ios;


// parse the information contained in a single spectrum file
//
// spectrumFileName = specrtum file name
// manager = pointer to the option manager
// spectrum = reference on the two dimensional array where the spectrum will be stored
// spectrumSize = number of sattes per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = reference on the minimum energy
// maxEnergy = reference on the maximum energy
// nbrStates = reference on the number of states
// return value = true if no error occured
bool DensityOfStatesParseSpectrumFile(char* spectrumFileName, OptionManager* manager, double**& spectrum, int*& spectrumSize, 
				      int*& spectrumWeight, int& nbrSectors, double& minEnergy, double& maxEnergy, long& nbrStates, Abstract1DRealFunction* densityOfStates);


// perform the density of states on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states
// nbrBins = number of bins for the density of states
// binSize = density of states range for each bin
// nbrStatePerBin = array that contains the number of level spacing per bin
// nbrRejectedStates = reference on the number of states (i.e. that cannot be stored in any bin)
// nbrAcceptedStates = reference on the number of states (i.e. that can be stored in a bin)
void DensityOfStatesPerformDensityOfStates(double** spectrum, int* spectrumSize, int* spectrumWeight, int nbrSectors, 
					   double minEnergy, double maxEnergy, long nbrStates,
					   int nbrBins, double binSize, long* nbrStatePerBin, long& nbrRejectedStates, long& nbrAcceptedStates);


// perform the density of states on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states
// nbrBins = number of bins for the density of states
// binSize = density of states range for each bin
// nbrStatePerBin = array that contains the number of level spacing per bin
// nbrRejectedStates = reference on the number of states (i.e. that cannot be stored in any bin)
// nbrAcceptedStates = reference on the number of states (i.e. that can be stored in a bin)
void DensityOfStatesParseSpectrumFile(ifstream& inputFile, double*& spectrum, int& spectrumSize,
				      double& minEnergy, double& maxEnergy, long& nbrStates,
				      int nbrBins, double binSize, long* nbrStatePerBin, long& nbrRejectedStates, long& nbrAcceptedStates);

// perform the density of states on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states
void DensityOfStatesParseSpectrumFile(ifstream& inputFile, double*& spectrum, int& spectrumSize,
				      double& minEnergy, double& maxEnergy, long& nbrStates);


// perform the density of states on a parsed spectrum with a cut-off on the entanglement energies
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states
// minEnergyCutOff = mininum entanglement energy cut-off
// minEnergyCutOff = maxinum entanglement energy cut-off
void DensityOfStatesParseSpectrumFile(ifstream& inputFile, double*& spectrum, int& spectrumSize,
				      double& minEnergy, double& maxEnergy, long& nbrStates, double minEnergyCutOff, double maxEnergyCutOff);

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
bool LevelStatisticsParseSpectrumFile(double* spectrum, int spectrumSize,
				      double& minAverageSpacing, double& maxAverageSpacing, double& averageSpacing,
				      long& nbrSpacings, Abstract1DRealFunction* densityOfStates);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("SpinChainMultipleEntanglementSpectrumLevelStatistics" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('s', "spectra", "name of the file that contains the binary entanglement spectra");
  (*SystemGroup) += new SingleIntegerOption('\n', "window-min", "set the index of the first entanglement spectrum to consider", 0);
  (*SystemGroup) += new SingleIntegerOption('\n', "window-max", " set the index of the last entanglement spectrum to consider (negative if up to the last available entanglement spectrum)", -1);
  (*SystemGroup) += new SingleDoubleOption('\n', "min-entenergy", "reject all entanglement energies below a given energy", 0.0);
  (*SystemGroup) += new SingleDoubleOption('\n', "max-entenergy", " reject all entanglement energies above a given energy (negative if no cut-off should be applied)", -1.0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the output file where the density of states will be stored (if none, try to deduce it from either --spectra replacing the dat extension with dos/levelstat)");
  (*OutputGroup) += new SingleIntegerOption ('\n', "nbr-bins", "number of bins for the density of states", 100);
  (*OutputGroup) += new SingleDoubleOption ('\n', "bin-spacing", "spacing range for each bin for the level statistics", 0.1);
  (*OutputGroup) += new BooleanOption ('\n', "discard-outputfiles", "do not save any results on disk, just display the summary using the standard output");
  (*OutputGroup) += new SingleIntegerOption ('\n', "extract-singlespectrum", "extract a single entanglement spectrum from --spectra instead of computing level statistics", -1);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainMultipleEntanglementSpectrumLevelStatistics -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if (Manager.GetString("spectra") == 0)
    {
      cout << "error, no entanglement spectrum was provided" << endl;
      return 0;
    }
  
  
  double MaxEnergy = 0.0;
  double MinEnergy = 1e300;
  double MinEnergyCutOff = 0.0;
  double MaxEnergyCutOff = 1e300;
  bool CutOffFlag = false;
  if ((Manager.GetDouble("min-entenergy") > 0.0) || (Manager.GetDouble("max-entenergy") > 0.0))
    {
      CutOffFlag  = true;
      MinEnergyCutOff = Manager.GetDouble("min-entenergy");
      if (Manager.GetDouble("max-entenergy") > 0.0)
	{
	  MaxEnergyCutOff = Manager.GetDouble("max-entenergy");
	}
    }
  long NbrLevels = 0l;
  int* SpectrumSize;
  double** Spectrum;
  int* SpectrumWeight;
  int NbrSectors;
  ifstream File;
  File.open(Manager.GetString("spectra"), ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "cannot open the file " << Manager.GetString("spectra") << endl;
      return 0;
    }

  int NbrEntanglementSpectra = 0;
  ReadLittleEndian(File, NbrEntanglementSpectra);  
  int NbrSzASectors = 0;
  ReadLittleEndian(File, NbrSzASectors);  
  int* SzaSectors = new int [NbrSzASectors];
  ReadBlockLittleEndian(File, SzaSectors, NbrSzASectors);

  int MinEntanglementSpectrumIndex = Manager.GetInteger("window-min");
  int MaxEntanglementSpectrumIndex = Manager.GetInteger("window-max");
  if (MaxEntanglementSpectrumIndex < MinEntanglementSpectrumIndex)
    {
      MaxEntanglementSpectrumIndex = NbrEntanglementSpectra - 1;
    }
  int TotalNbrEntanglementSpectra = (MaxEntanglementSpectrumIndex - MinEntanglementSpectrumIndex + 1);
//  int TotalNbrEntanglementSpectra = NbrSzASectors * (MaxEntanglementSpectrumIndex - MinEntanglementSpectrumIndex + 1);
  if (Manager.GetInteger("extract-singlespectrum") >= 0)
    {
      MaxEntanglementSpectrumIndex = Manager.GetInteger("extract-singlespectrum");
      MinEntanglementSpectrumIndex = MaxEntanglementSpectrumIndex;
    }
  cout << "number of entanglement spectra = " << NbrEntanglementSpectra << endl;
  cout << "will process the entanglement spectra between " << MinEntanglementSpectrumIndex << " and " << MaxEntanglementSpectrumIndex << endl;
  Spectrum = new double*[TotalNbrEntanglementSpectra];
  SpectrumSize = new int[TotalNbrEntanglementSpectra];
  int Index = 0;

  for (int i = 0; i < MinEntanglementSpectrumIndex; ++i)
    {
      for (int j = 0; j < NbrSzASectors; ++j)
	{
	  int TmpSize = 0;
	  ReadLittleEndian(File, TmpSize);
	  File.seekg (TmpSize * sizeof(double), ios::cur);
	}
    }
  if (Manager.GetInteger("extract-singlespectrum") >= 0)
    {
      char* SingleSpectrumExtension = new char[128];
      sprintf (SingleSpectrumExtension, "spec_%d.ent", MaxEntanglementSpectrumIndex);
      char* SingleSpectrumOutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectra"), "ent", SingleSpectrumExtension);
      delete[] SingleSpectrumExtension;
      ofstream File2;
      File2.open(SingleSpectrumOutputFileName, ios::binary | ios::out);
      File2.precision(14);
      for (int j = 0; j < NbrSzASectors; ++j)
	{
	  int TmpSize = 0;
	  ReadLittleEndian(File, TmpSize);
	  double* TmpArray = new double[TmpSize];
	  ReadBlockLittleEndian(File, TmpArray, TmpSize); 
	  for (int k = 0; k <  TmpSize; ++k)
	    File2 << SzaSectors[j] << " " << TmpArray[k] << " " << (-log(TmpArray[k])) << endl;
	  delete[] TmpArray;
	}      
      File2.close();
      File.close();
      return 0;
    }
  if (true)
    {
       for (int i = MinEntanglementSpectrumIndex; i <= MaxEntanglementSpectrumIndex; ++i)
	{
	  for (int j = 0; j < NbrSzASectors; ++j)
	    {
	      if (SzaSectors[j] == 0)
		{
		  if (CutOffFlag == false)
		    DensityOfStatesParseSpectrumFile(File, Spectrum[Index], SpectrumSize[Index], MinEnergy, MaxEnergy, NbrLevels);
		  else
		    DensityOfStatesParseSpectrumFile(File, Spectrum[Index], SpectrumSize[Index], MinEnergy, MaxEnergy, NbrLevels, MinEnergyCutOff, MaxEnergyCutOff);
		  ++Index;
		}
	      else
		{
		  int TmpSize = 0;
		  ReadLittleEndian(File, TmpSize);
		  File.seekg (TmpSize * sizeof(double), ios::cur);
		}
	    }
	}
   }
  else
    {
      for (int i = MinEntanglementSpectrumIndex; i <= MaxEntanglementSpectrumIndex; ++i)
	{
	  for (int j = 0; j < NbrSzASectors; ++j)
	    {
	      DensityOfStatesParseSpectrumFile(File, Spectrum[Index], SpectrumSize[Index], MinEnergy, MaxEnergy, NbrLevels);
	      ++Index;
	    }
	}
    }
  File.close();

  char* DensityOfStateOutputFileName = 0;
  char* LevelStatisticOutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      char* DensityOfStateExtension = new char[128];
      sprintf (DensityOfStateExtension, "range_%d_%d.ent.dos", MinEntanglementSpectrumIndex, MaxEntanglementSpectrumIndex);
      DensityOfStateOutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectra"), "ent", DensityOfStateExtension);
      delete[] DensityOfStateExtension;
      char*  LevelStatisticExtension = new char[128];
      sprintf (LevelStatisticExtension, "range_%d_%d.ent.levelstat", MinEntanglementSpectrumIndex, MaxEntanglementSpectrumIndex);
      LevelStatisticOutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectra"), "ent", LevelStatisticExtension);
      delete[] LevelStatisticExtension;
    }
  else
    {
      LevelStatisticOutputFileName = new char [strlen (Manager.GetString("output-file")) + 1];
      strcpy (LevelStatisticOutputFileName, Manager.GetString("output-file"));
    }

  
  int NbrBins = Manager.GetInteger("nbr-bins");
  MaxEnergy +=  0.001 * (MaxEnergy - MinEnergy) / ((double) NbrBins);
  double BinSize = (MaxEnergy - MinEnergy) / ((double) NbrBins);
  long* NbrStatePerBin = new long [NbrBins];
  for (int i = 0 ; i < NbrBins; ++i)
    NbrStatePerBin[i] = 0l;
  long NbrRejectedStates = 0l;
  long NbrAcceptedStates = 0l;

  for (int i = 0; i < TotalNbrEntanglementSpectra; ++i)
    {
      for (int j = 0; j < SpectrumSize[i]; ++j)
	{
	  int TmpIndex = int ((Spectrum[i][j] - MinEnergy) / BinSize);
	  if (TmpIndex < NbrBins)
	    {
	      NbrStatePerBin[TmpIndex]++;
	      ++NbrAcceptedStates;
	    }
	  else
	    {
	      ++NbrRejectedStates;
	    }      
	}
    }


  long MaxNbrStatePerBin = 0l;
  long MinNbrStatePerBin = NbrLevels;
  for (int i = 0 ; i < NbrBins; ++i)
    {
      if (NbrStatePerBin[i] > MaxNbrStatePerBin)
	{
	  MaxNbrStatePerBin = NbrStatePerBin[i];
	}
      if (NbrStatePerBin[i] < MinNbrStatePerBin)
	{
	  MinNbrStatePerBin = NbrStatePerBin[i];
	}
    }
  
  double Sum = 0.0;
  double* DensityOfStateEnergies =  new double[NbrBins];
  double* DensityOfStateValues =  new double[NbrBins];
  for (int i = 0 ; i < NbrBins; ++i)
    {
      double DOS = (((double) NbrStatePerBin[i]) / ((double) NbrLevels));
      DensityOfStateEnergies[i] = (MinEnergy + (BinSize * ((double) i)));
      DensityOfStateValues[i] = (DOS / BinSize);
      Sum += DOS;
    }

  if ((DensityOfStateOutputFileName != 0) && (Manager.GetBoolean("discard-outputfiles") == false))
    {
      ofstream File;
      File.open(DensityOfStateOutputFileName, ios::binary | ios::out);
      File.precision(14);
      File << "# Min energy " << MinEnergy << endl;
      File << "# Max energy " << MaxEnergy << endl;
      File << "# Nbr of states = " << NbrLevels << endl;
      File << "# Nbr rejected states = " << NbrRejectedStates << endl;
      File << "# int dE rho(E) = " << Sum << endl;
      File << "# Min nbr of states per bin " << MinNbrStatePerBin << endl;
      File << "# Max nbr of states per bin " << MaxNbrStatePerBin << endl;
      File << "# Bin size " << BinSize << endl;
      File << "# E P(E)" << endl;
      for (int i = 0 ; i < NbrBins; ++i)
	{
	  double DOS = (((double) NbrStatePerBin[i]) / ((double) NbrLevels));
	  File << (MinEnergy + (BinSize * ((double) i))) << " " << (DOS / BinSize) << endl;
	}  
      File.close();
    }

  cout << "Min energy " << MinEnergy << endl;
  cout << "Max energy " << MaxEnergy << endl;
  cout << "Nbr of states = " << NbrLevels << endl;
  cout << "Nbr rejected states = " << NbrRejectedStates << endl;
  cout << "int dE rho(E) = " << Sum << endl;
  cout << "Min nbr of states per bin " << MinNbrStatePerBin << endl;
  cout << "Max nbr of states per bin " << MaxNbrStatePerBin << endl;
  cout << "Bin size " << BinSize << endl;

  Abstract1DRealFunction* DensityOfStates = 0;
  double DensityOfStatesThreshold = 1e-5;
  Tabulated1DRealFunction TmpFunction(DensityOfStateEnergies, DensityOfStateValues, NbrBins);
  DensityOfStates = TmpFunction.GetPrimitive();
  //      DensityOfStates = new Linear1DRealFunction(1.0);

  double AverageSpacing = 0.0;
  double MinAverageSpacing = 1.0e300;
  double MaxAverageSpacing = 0.0;
  long NbrSpacings = 0l;
  for (int i = 0; i < TotalNbrEntanglementSpectra; ++i)
    {
      LevelStatisticsParseSpectrumFile(Spectrum[i], SpectrumSize[i],
				       MinAverageSpacing, MaxAverageSpacing, AverageSpacing, NbrSpacings, DensityOfStates);
    }


  cout << MinAverageSpacing << " " << MaxAverageSpacing << " " << AverageSpacing << endl;
  AverageSpacing /= (double) NbrSpacings;
  MinAverageSpacing /= AverageSpacing;
  MaxAverageSpacing /= AverageSpacing;

  double SpacingBinSize = Manager.GetDouble("bin-spacing");
  int SpacingNbrBins = ((int) (MaxAverageSpacing / SpacingBinSize)) + 1;
  long* NbrSpacingPerBin = new long [SpacingNbrBins];
  for (int i = 0 ; i < SpacingNbrBins; ++i)
    NbrSpacingPerBin[i] = 0l;
  long NbrRejectedSpacings = 0l;
  long NbrAcceptedSpacings = 0l;
  double AverageR = 0.0;
  double VarianceAverageR = 0.0;

  for (int i = 0; i < TotalNbrEntanglementSpectra; ++i)
    {
      if (SpectrumSize[i] > 1)
	{
	  double TmpAverageR = 0.0;
	  double TmpNbrRatios = 0.0;
	  double TmpVarianceAverageR = 0.0;
	  int Lim = SpectrumSize[i] - 1;	  
	  double TmpInfDiff;
	  double TmpSupDiff;
	  for (int j = 1; j < Lim; ++j)
	    {
	      double TmpDiff = (Spectrum[i][j] - Spectrum[i][j - 1]) / AverageSpacing;
	      int TmpIndex = int (TmpDiff / SpacingBinSize);
	      if (TmpIndex < SpacingNbrBins)
		{
		  NbrSpacingPerBin[TmpIndex]++;
		  ++NbrAcceptedSpacings;
		}
	      else
		{
		  ++ NbrRejectedSpacings;
		}
	      double TmpDiff2 = (Spectrum[i][j + 1] - Spectrum[i][j]) / AverageSpacing;
	      if (TmpDiff2 > TmpDiff)
		{
		  if (TmpDiff2 > MACHINE_PRECISION)
		    {
		      double Tmp = TmpDiff / TmpDiff2;
		      TmpAverageR += Tmp;	
		      TmpVarianceAverageR += Tmp * Tmp;	      
		      TmpNbrRatios += 1.0;
		    }
		}
	      else
		{
		  if (TmpDiff > MACHINE_PRECISION)
		    {
		      double Tmp =  TmpDiff2 / TmpDiff;
		      TmpAverageR += Tmp;	
		      TmpVarianceAverageR += Tmp * Tmp;	      
		      TmpNbrRatios += 1.0;
		    }
		}
	    }
	  double TmpDiff = (Spectrum[i][Lim] - Spectrum[i][Lim - 1]) / AverageSpacing;
	  int TmpIndex = int (TmpDiff / SpacingBinSize);
	  if (TmpIndex < SpacingNbrBins)
	    {
	      NbrSpacingPerBin[TmpIndex]++;
	      ++NbrAcceptedSpacings;
	    }
	  else
	    {
	      ++NbrRejectedSpacings;
	    }
	  if (TmpNbrRatios > 0.0)
	    {
	      AverageR += (TmpAverageR / TmpNbrRatios);
	      VarianceAverageR += (TmpAverageR / TmpNbrRatios) * (TmpAverageR / TmpNbrRatios);
	    }
	}
    }

  Sum = 0.0;
  for (int i = 0 ; i < SpacingNbrBins; ++i)
    {
      double PSpacing = (((double) NbrSpacingPerBin[i]) / ((double) NbrSpacings));
      Sum += PSpacing;
    }

  long MaxNbrSpacingPerBin = 0l;
  long MinNbrSpacingPerBin = NbrSpacings;
  for (int i = 0 ; i < SpacingNbrBins; ++i)
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

  if (Manager.GetBoolean("discard-outputfiles") == false)
    {
      ofstream File2;
      File2.open(LevelStatisticOutputFileName, ios::binary | ios::out);
      File2.precision(14);
      File2 << "# Min spacing " << MinAverageSpacing << endl;
      File2 << "# Max spacing " << MaxAverageSpacing << endl;
      File2 << "# Average spacing = " << AverageSpacing << endl;
      File2 << "# Nbr points = " << NbrSpacings << endl;
      File2 << "# Nbr rejected points = " << NbrRejectedSpacings << endl;
      File2 << "# int ds P(s) = " << Sum << endl;
      File2 << "# Min nbr spacings per bin " << MinNbrSpacingPerBin << endl;
      File2 << "# Max nbr spacings per bin " << MaxNbrSpacingPerBin << endl;
      File2 << "# Max spacing " << MaxAverageSpacing << endl;
      File2 << "# s P(s)" << endl;
      File2 << "# <r>=" << (AverageR / TotalNbrEntanglementSpectra) << " " 
	    << sqrt (((VarianceAverageR / TotalNbrEntanglementSpectra) - ((AverageR / TotalNbrEntanglementSpectra) * (AverageR / TotalNbrEntanglementSpectra))) / (TotalNbrEntanglementSpectra - 1)) << endl;
      for (int i = 0 ; i < SpacingNbrBins; ++i)
	{
	  double PSpacing = (((double) NbrSpacingPerBin[i]) / ((double) NbrSpacings));
	  File2 << (SpacingBinSize * (((double) i) + 0.5)) << " " << (PSpacing / SpacingBinSize) << endl;
	}  
      File2.close();
    }

  cout << "Min spacing " << MinAverageSpacing << endl;
  cout << "Max spacing " << MaxAverageSpacing << endl;
  cout << "Average spacing = " << AverageSpacing << endl;
  cout << "Nbr points = " << NbrSpacings << endl;
  cout << "Nbr rejected points = " << NbrRejectedSpacings << endl;
  cout << "int ds P(s) = " << Sum << endl;
  cout << "Min nbr spacings per bin " << MinNbrSpacingPerBin << endl;
  cout << "Max nbr spacings per bin " << MaxNbrSpacingPerBin << endl;
  cout << "<r>=" << (AverageR / TotalNbrEntanglementSpectra) << " " 
       << sqrt (((VarianceAverageR / TotalNbrEntanglementSpectra) - ((AverageR / TotalNbrEntanglementSpectra) * (AverageR / TotalNbrEntanglementSpectra))) / (TotalNbrEntanglementSpectra - 1)) << endl;
  delete[] NbrSpacingPerBin;
  delete[] NbrStatePerBin;

  return 0;
}

// parse the information contained in a single spectrum file
//
// spectrumFileName = specrtum file name
// manager = pointer to the option manager
// spectrum = reference on the two dimensional array where the spectrum will be stored
// spectrumSize = number of sattes per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = reference on the minimum energy
// maxEnergy = reference on the maximum energy
// nbrStates = reference on the number of states
// return value = true if no error occured

bool DensityOfStatesParseSpectrumFile(char* spectrumFileName, OptionManager* manager, double**& spectrum, int*& spectrumSize, 
				      int*& spectrumWeight, int& nbrSectors, double& minEnergy, double& maxEnergy, long& nbrStates, Abstract1DRealFunction* densityOfStates)
{
  MultiColumnASCIIFile SpectrumFile;
  if (SpectrumFile.Parse(spectrumFileName) == false)
    {
      SpectrumFile.DumpErrors(cout);
      return false;
    }
  nbrSectors = 1;
  if ((manager->GetInteger("energy-column") == 0) || (manager->GetBoolean("discard-quantumnumbers") == true))
    {
      spectrumSize = new int [1];
      spectrum = new double*[1];
      spectrumWeight = new int[1]; 
      double* TmpSpectrum = SpectrumFile.GetAsDoubleArray(manager->GetInteger("energy-column"));
      int TmpSize = SpectrumFile.GetNbrLines();
      if (manager->GetBoolean("discard-quantumnumbers") == true)
	{
	  SortArrayUpOrdering<double>(TmpSpectrum, TmpSize);
	}
      if ((manager->GetInteger("window-min") > 0) || (manager->GetInteger("window-max") >= 0))
	{
	  int MaxIndex = manager->GetInteger("window-max");
	  if (MaxIndex < 0)
	    MaxIndex = TmpSize - 1;
	  int MinIndex = manager->GetInteger("window-min");
	  spectrumSize[0] = MaxIndex - MinIndex + 1;
	  spectrumWeight[0] = 1;
	  spectrum[0] = new double[spectrumSize[0]];
	  for (int i = MinIndex; i <= MaxIndex; ++i)
	    spectrum[0][i - MinIndex] = TmpSpectrum[i];
	  delete[] TmpSpectrum;
	}
      else
	{
	  spectrum[0] = TmpSpectrum;
	  spectrumSize[0] = TmpSize;
	  spectrumWeight[0] = 1;
	}
   }
  else
    {
      int TmpSize = SpectrumFile.GetNbrLines();
      int NbrQuantumNumber = manager->GetInteger("energy-column");
      int** QuantumNumbers =  new int*[NbrQuantumNumber];
      double* TmpSpectrum = 0;
      int* TmpDegeneracy = 0;
      if (manager->GetInteger("filter-column") >= 0)
	{
	  int FilterColumn = manager->GetInteger("filter-column");
	  int FilterValue = manager->GetInteger("filter-value");
	  int TmpActualSize = 0;
	  int* TmpFilterArray = SpectrumFile.GetAsIntegerArray(FilterColumn);
	  for (int i = 0; i < TmpSize; ++i)
	    if (TmpFilterArray[i] == FilterValue)
	      ++TmpActualSize;
	  if (TmpActualSize == 0)
	    {
	      cout << "error, applying filter leads to an empty spectrum" << endl;
	      return false;
	    }
	  for (int i = 0; i < NbrQuantumNumber; ++i)
	    {
	      QuantumNumbers[i] = new int [TmpActualSize];	  
	      int* TmpArray = SpectrumFile.GetAsIntegerArray(i);
	      TmpActualSize = 0;
	      for (int j = 0; j < TmpSize; ++j)
		{
		  if (TmpFilterArray[j] == FilterValue)
		    {
		      QuantumNumbers[i][TmpActualSize] = TmpArray[j];
		      ++TmpActualSize;
		    }
		}
	      delete[] TmpArray;
	    }
	  double* TmpSpectrum2 = SpectrumFile.GetAsDoubleArray(manager->GetInteger("energy-column"));
	  TmpSpectrum = new double[TmpActualSize];
	  TmpActualSize = 0;
	  for (int j = 0; j < TmpSize; ++j)
	    {
	      if (TmpFilterArray[j] == FilterValue)
		{
		  TmpSpectrum[TmpActualSize] = TmpSpectrum2[j];
		  ++TmpActualSize;
		}
	    }	  
	  if (manager->GetInteger("degeneracy-column") >= 0)
	    {
	      TmpDegeneracy = new int[TmpActualSize];
	      int* TmpDegeneracy2 = SpectrumFile.GetAsIntegerArray(manager->GetInteger("degeneracy-column"));
	      TmpActualSize = 0;
	      for (int j = 0; j < TmpSize; ++j)
		{
		  if (TmpFilterArray[j] == FilterValue)
		    {
		      TmpDegeneracy[TmpActualSize] = TmpDegeneracy2[j];
		      ++TmpActualSize;
		    }
		}	  
	      delete[] TmpDegeneracy2;
	    }
	  delete[] TmpSpectrum2;
	  delete[] TmpFilterArray;
	  TmpSize = TmpActualSize;
	}     
      else
	{
	  for (int i = 0; i < NbrQuantumNumber; ++i)
	    QuantumNumbers[i] = SpectrumFile.GetAsIntegerArray(i);
	  TmpSpectrum = SpectrumFile.GetAsDoubleArray(manager->GetInteger("energy-column"));
	  if (manager->GetInteger("degeneracy-column") >= 0)
	    {
	      TmpDegeneracy = SpectrumFile.GetAsIntegerArray(manager->GetInteger("degeneracy-column"));
	    }
	}
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
      spectrum[0][0] = (*densityOfStates)(TmpSpectrum[0]);
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
	  spectrum[nbrSectors][i - CurrentIndex] = (*densityOfStates)(TmpSpectrum[i]);
	}
      ++nbrSectors;            
      for (int i = 0; i < NbrQuantumNumber; ++i)
	delete[] QuantumNumbers[i];
      delete[] QuantumNumbers;
    }

  for (int i = 0; i < nbrSectors; ++i)
    {     
      if (spectrumSize[i] > 0)
	{
	  nbrStates += spectrumSize[i];
	  int Lim = spectrumSize[i];	  
	  for (int j = 0; j < Lim; ++j)
	    {
	      double TmpEnergy =  spectrum[i][j];
	      if (TmpEnergy > maxEnergy)
		{
		  maxEnergy = TmpEnergy;
		}
	      if (TmpEnergy < minEnergy)
		{
		  minEnergy = TmpEnergy;
		}
	    }
	}
    }
  return true;
}

// perform the density of states on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states
// nbrBins = number of bins for the density of states
// binSize = density of states range for each bin
// nbrStatePerBin = array that contains the number of level spacing per bin
// nbrRejectedStates = reference on the number of states (i.e. that cannot be stored in any bin)
// nbrAcceptedStates = reference on the number of states (i.e. that can be stored in a bin)

void DensityOfStatesParseSpectrumFile(ifstream& inputFile, double*& spectrum, int& spectrumSize,
				      double& minEnergy, double& maxEnergy, long& nbrStates,
				      int nbrBins, double binSize, long* nbrStatePerBin, long& nbrRejectedStates, long& nbrAcceptedStates)
{
  ReadLittleEndian(inputFile, spectrumSize);
  spectrum = new double[spectrumSize];
  ReadBlockLittleEndian(inputFile, spectrum, spectrumSize);  
  nbrStates += (long) spectrumSize;
  for (int i = 0; i < spectrumSize; ++i)
    {     
      if (spectrum[i] < minEnergy)
	minEnergy = spectrum[i];
      if (spectrum[i] > maxEnergy)
	maxEnergy = spectrum[i];
      int TmpIndex = int (spectrum[i] / binSize);
      if (TmpIndex < nbrBins)
	{
	  nbrStatePerBin[TmpIndex]++;
	  ++nbrAcceptedStates;
	}
      else
	{
	  ++ nbrRejectedStates;
	}      
    }  
}

// perform the density of states on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states

void DensityOfStatesParseSpectrumFile(ifstream& inputFile, double*& spectrum, int& spectrumSize,
				      double& minEnergy, double& maxEnergy, long& nbrStates)
{
  ReadLittleEndian(inputFile, spectrumSize);
  spectrum = new double[spectrumSize];
  ReadBlockLittleEndian(inputFile, spectrum, spectrumSize);  
  nbrStates += (long) spectrumSize;
  for (int i = 0; i < spectrumSize; ++i)
    {     
      spectrum[i] = -log(spectrum[i]);
      if (spectrum[i] < minEnergy)
	minEnergy = spectrum[i];
      if (spectrum[i] > maxEnergy)
	maxEnergy = spectrum[i];
     }  
}

// perform the density of states on a parsed spectrum with a cut-off on the entanglement energies
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states
// minEnergyCutOff = mininum entanglement energy cut-off
// minEnergyCutOff = maxinum entanglement energy cut-off

void DensityOfStatesParseSpectrumFile(ifstream& inputFile, double*& spectrum, int& spectrumSize,
				      double& minEnergy, double& maxEnergy, long& nbrStates, double minEnergyCutOff, double maxEnergyCutOff)
{
  int TmpSpectrumSize = 0;
  ReadLittleEndian(inputFile, TmpSpectrumSize);
  double* TmpSpectrum = new double[TmpSpectrumSize];
  ReadBlockLittleEndian(inputFile, TmpSpectrum, TmpSpectrumSize);  
  spectrumSize = 0;
  for (int i = 0; i < TmpSpectrumSize; ++i)
    {     
      TmpSpectrum[i] = -log(TmpSpectrum[i]);
      if ((TmpSpectrum[i] >= minEnergyCutOff) && (TmpSpectrum[i] <= maxEnergyCutOff))
	++spectrumSize;
    }
  spectrum = new double[spectrumSize];
  nbrStates += (long) spectrumSize;
  spectrumSize = 0;
   for (int i = 0; i < TmpSpectrumSize; ++i)
    {     
      if ((TmpSpectrum[i] >= minEnergyCutOff) && (TmpSpectrum[i] <= maxEnergyCutOff))
	{
	  if (TmpSpectrum[i] < minEnergy)
	    minEnergy = TmpSpectrum[i];
	  if (TmpSpectrum[i] > maxEnergy)
	    maxEnergy = TmpSpectrum[i];
	  spectrum[spectrumSize] = TmpSpectrum[i];
	  ++spectrumSize;
	}  
    }
   delete[] TmpSpectrum;
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

bool LevelStatisticsParseSpectrumFile(double* spectrum, int spectrumSize,
				      double& minAverageSpacing, double& maxAverageSpacing, double& averageSpacing,
				      long& nbrSpacings, Abstract1DRealFunction* densityOfStates)
{
  for (int i = 0; i < spectrumSize; ++i)
    spectrum[i] = (*densityOfStates)(spectrum[i]);
  nbrSpacings += spectrumSize - 1;
  for (int j = 1; j < spectrumSize; ++j)
    {
      double TmpDiff = spectrum[j] - spectrum[j - 1];
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
  return true;
}

