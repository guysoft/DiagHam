#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ArrayTools.h"

#include "Tools/FQHEFiles/FQHEOnCylinderFileTools.h"

#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpin.h"

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


// extract the spectrum from a file
//
// spectrumFile = spectrum file name
// nbrLzSectors = reference on the number of momentum sectors
// minLzSector = reference on the minimal momentum
// nbrEnergie =  reference on the array where the number of levels per momentum sector
// energies =  reference on the array where the emergies will be stored
// return value = true if no error occured
bool FQHECylinderQuasiholesWithSpinTimeReversalSymmetryExtractSpectrum(char* spectrumFile, int& nbrLzSectors, int& minLzSector,
								       int*& nbrEnergie, double**& energies);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis" , "0.01");
  OptionGroup* SystemGroup  = new OptionGroup("system options");
  OptionGroup* OutputGroup  = new OptionGroup ("ouput options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('s', "singlelayer-spectra", "name of the file containing the list of single layer spectra");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-states", "maximum number of states to generate per momentum sector", 100);  
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the two layer system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-ky", "twice the total momentum of the two layer system", 0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "charging-energy", "factor in front of the charging energy (i.e 1/(2C))", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "average-nbrparticles", "average number of particles", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "fix-nbrparticles", "fix the number of particles for the two layer system (no restriction if negative)", -1);  
  (*SystemGroup) += new SingleDoubleOption ('\n', "degeneracy-error", "difference below which two energies are considered to be degenerate", 1.0e-12);
  (*SystemGroup) += new SingleStringOption ('\n', "directory", "use a specific directory for the input data instead of the current one (only useful when building the eigenstates in the full quasihole basis)");
  
  (*OutputGroup) += new BooleanOption ('\n', "build-eigenstates", "build the eigenstates in the full quasihole basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("singlelayer-spectra") == 0)
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis -h" << endl;
      return -1;
    }

  int TwoLayerSzSector = Manager.GetInteger("total-sz");
  int TwoLayerLzSector = Manager.GetInteger("total-ky");
  int FixTotalNbrParticles = Manager.GetInteger("fix-nbrparticles");

  int MaxNbrParticlesPerLayer = -1;
  int NbrFluxQuanta = 0;
  bool Statistics = true;
  double Ratio = 0.0;
  double Perimeter = 0.0;
  int KValue = 1;
  int RValue = 2;

  MultiColumnASCIIFile SingleLayerSpectrumFile;
  if (SingleLayerSpectrumFile.Parse(Manager.GetString("singlelayer-spectra")) == false)
    {
      SingleLayerSpectrumFile.DumpErrors(cout);
      return -1;
    }
  for (int i = 0; i < SingleLayerSpectrumFile.GetNbrLines(); ++i)
    {
      int TmpNbrParticles = 0;
      if (FQHEOnCylinderFindSystemInfoFromFileName(SingleLayerSpectrumFile(0, i), TmpNbrParticles, NbrFluxQuanta,
						   Statistics, Ratio, Perimeter) == false)
	{
	  cout << "can't extract system information from file name " << SingleLayerSpectrumFile(0, i) << endl;
	  return -1;
	}
      if (TmpNbrParticles > MaxNbrParticlesPerLayer)
	{
	  MaxNbrParticlesPerLayer = TmpNbrParticles;
	}
    }
  char** SpectrumFileNames = new char* [MaxNbrParticlesPerLayer + 1];
  char* TmpOutputFileName;
  int* NbrLzSectors = new int [MaxNbrParticlesPerLayer + 1];
  int* MinLzValues = new int [MaxNbrParticlesPerLayer + 1];
  int** NbrEnergies = new int* [MaxNbrParticlesPerLayer + 1];
  double*** Energies = new double** [MaxNbrParticlesPerLayer + 1];
  for (int i = 0; i <= MaxNbrParticlesPerLayer; ++i)
    {
      NbrLzSectors[i] = 0;
    }
  
  for (int i = 0; i < SingleLayerSpectrumFile.GetNbrLines(); ++i)
    {
      int TmpNbrParticles = 0;
      if (FQHEOnCylinderFindSystemInfoFromFileName(SingleLayerSpectrumFile(0, i), TmpNbrParticles, NbrFluxQuanta,
						   Statistics, Ratio, Perimeter) == false)
	{
	  cout << "can't extract system information from file name " << SingleLayerSpectrumFile(0, i) << endl;
	  return -1;
	}
      if (FQHECylinderQuasiholesWithSpinTimeReversalSymmetryExtractSpectrum(SingleLayerSpectrumFile(0, i), NbrLzSectors[TmpNbrParticles], MinLzValues[TmpNbrParticles],
									    NbrEnergies[TmpNbrParticles], Energies[TmpNbrParticles]) == false)
	{
	  return -1;
	}
      SpectrumFileNames[TmpNbrParticles] = new char[strlen(SingleLayerSpectrumFile(0, i)) + 1];
      strcpy (SpectrumFileNames[TmpNbrParticles], SingleLayerSpectrumFile(0, i));
      if (TmpNbrParticles == 0)
	{
	  char* TmpString = new char[128];
	  sprintf(TmpString, "sz_%d_lz_%d", TwoLayerSzSector, TwoLayerLzSector);
	  TmpOutputFileName = ReplaceString(SingleLayerSpectrumFile(0, i), "sz_0_lz", TmpString);
	  delete[] TmpString;
	}
    }

  int TotalNbrLevels = 0;

  for (int UpLayerNbrParticles = 0; UpLayerNbrParticles <= MaxNbrParticlesPerLayer; ++UpLayerNbrParticles)
    {
      if (NbrLzSectors[UpLayerNbrParticles] != 0)
	{
	  int DownLayerNbrParticles = UpLayerNbrParticles - TwoLayerSzSector;
	  if ((DownLayerNbrParticles >= 0) && (DownLayerNbrParticles <= MaxNbrParticlesPerLayer) && (NbrLzSectors[DownLayerNbrParticles] != 0) && ((FixTotalNbrParticles < 0) || (FixTotalNbrParticles == (DownLayerNbrParticles + UpLayerNbrParticles))))
	    {
	      for (int UpLayerLzSector = 0; UpLayerLzSector < NbrLzSectors[UpLayerNbrParticles]; ++UpLayerLzSector)
		{
		  int DownLayerLzSector = ((((UpLayerLzSector * 2) + MinLzValues[UpLayerNbrParticles]) - TwoLayerLzSector) -  MinLzValues[DownLayerNbrParticles]) / 2;
		  TotalNbrLevels += NbrEnergies[UpLayerNbrParticles][UpLayerLzSector] *  NbrEnergies[DownLayerNbrParticles][DownLayerLzSector];
		}
	    }
	} 
    }
  int* TwoLayerIndices = new int [TotalNbrLevels];
  int* TwoLayerNbrParticles = new int [TotalNbrLevels];
  int* TwoLayerUpIndices = new int [TotalNbrLevels];
  int* TwoLayerUpLzValues = new int [TotalNbrLevels];
  int* TwoLayerDownIndices = new int [TotalNbrLevels];
  double* TwoLayerEnergies = new double [TotalNbrLevels];
  TotalNbrLevels = 0;
  for (int UpLayerNbrParticles = 0; UpLayerNbrParticles <= MaxNbrParticlesPerLayer; ++UpLayerNbrParticles)
    {
      if (NbrLzSectors[UpLayerNbrParticles] != 0)
	{
	  int DownLayerNbrParticles = UpLayerNbrParticles - TwoLayerSzSector;
	  if ((DownLayerNbrParticles >= 0) && (DownLayerNbrParticles <= MaxNbrParticlesPerLayer) && (NbrLzSectors[DownLayerNbrParticles] != 0) && ((FixTotalNbrParticles < 0) || (FixTotalNbrParticles == (DownLayerNbrParticles + UpLayerNbrParticles))))
	    {
	      int TmpNbrParticles = UpLayerNbrParticles + DownLayerNbrParticles;
	      double TmpEnergyShift = (Manager.GetDouble("average-nbrparticles") - ((double) TmpNbrParticles));
	      TmpEnergyShift *= TmpEnergyShift;	      
	      TmpEnergyShift *= Manager.GetDouble("charging-energy");
	      for (int UpLayerLzSector = 0; UpLayerLzSector < NbrLzSectors[UpLayerNbrParticles]; ++UpLayerLzSector)
		{
		  int ShiftedUpLayerLzSector = ((UpLayerLzSector * 2) + MinLzValues[UpLayerNbrParticles]);
		  int DownLayerLzSector = ((ShiftedUpLayerLzSector - TwoLayerLzSector) - MinLzValues[DownLayerNbrParticles]) / 2;
		  for (int i = 0; i < NbrEnergies[UpLayerNbrParticles][UpLayerLzSector]; ++i)
		    {
		      for (int j = 0; j < NbrEnergies[DownLayerNbrParticles][DownLayerLzSector]; ++j)
			{
			  TwoLayerEnergies[TotalNbrLevels] = Energies[UpLayerNbrParticles][UpLayerLzSector][i] + Energies[DownLayerNbrParticles][DownLayerLzSector][j] + TmpEnergyShift;
			  TwoLayerNbrParticles[TotalNbrLevels] = TmpNbrParticles;
			  TwoLayerIndices[TotalNbrLevels] = TotalNbrLevels;
			  TwoLayerUpLzValues[TotalNbrLevels] = ShiftedUpLayerLzSector;
			  TwoLayerUpIndices[TotalNbrLevels] = i;
			  TwoLayerDownIndices[TotalNbrLevels] = j;
			  ++TotalNbrLevels;
			}
		    }
		}
	    }
	} 
    }

  cout << "total number of levels in the Sz=" << TwoLayerSzSector << " Ky=" << TwoLayerLzSector << " : " << TotalNbrLevels << endl; 
  SortArrayUpOrdering<int>(TwoLayerEnergies, TwoLayerIndices, TotalNbrLevels);
  int EffectiveSubspaceDimension = Manager.GetInteger("nbr-states");
  double Error = Manager.GetDouble("degeneracy-error");
  if (EffectiveSubspaceDimension == TotalNbrLevels)
    {
      --EffectiveSubspaceDimension;
    }
  while ((EffectiveSubspaceDimension > 0) && (fabs (TwoLayerEnergies[EffectiveSubspaceDimension - 1] - TwoLayerEnergies[EffectiveSubspaceDimension]) < Error))
    {
      --EffectiveSubspaceDimension;
    }
  
  int** LargestUsedEigenstate = new int* [MaxNbrParticlesPerLayer + 1];
  for (int i = 0; i <= MaxNbrParticlesPerLayer; ++i)
    {
      LargestUsedEigenstate[i] = new int [NbrLzSectors[i]];
      for (int j = 0; j < NbrLzSectors[i]; ++j)
	{
	  LargestUsedEigenstate[i][j] = -1;
	}
    }
  char* OutputFileNameExtension = new char[128];
  sprintf (OutputFileNameExtension, "_effective_%d.basis", EffectiveSubspaceDimension);
  char* OutputFileName = ReplaceString(TmpOutputFileName, ".dat", OutputFileNameExtension);
  ofstream File;  
  File.open(OutputFileName, ios::binary | ios::out); 
  File.precision(14); 
  File << "# N Sz Lz N_u N_d Lz_u Lz_d E file_up file_down" << endl;
  char* TmpExtension = new char[128];
  for (int i = 0; i < EffectiveSubspaceDimension; ++i)
    {
      int UpLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] + TwoLayerSzSector) / 2;
      int DownLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] - TwoLayerSzSector) / 2;  
      int DownLayerLzValue = TwoLayerUpLzValues[TwoLayerIndices[i]] - TwoLayerLzSector;
      File << TwoLayerNbrParticles[TwoLayerIndices[i]] << " " << TwoLayerSzSector << " " << TwoLayerLzSector
	   << " " << UpLayerNbrParticles << " " << DownLayerNbrParticles << " "
	   << TwoLayerUpLzValues[TwoLayerIndices[i]] << " " << DownLayerLzValue << " " << TwoLayerEnergies[i];
      sprintf(TmpExtension, "_%d.%d.vec", TwoLayerUpLzValues[TwoLayerIndices[i]], TwoLayerUpIndices[TwoLayerIndices[i]]);
      char* TmpUpLayerVectorFileName = ReplaceExtensionToFileName(SpectrumFileNames[UpLayerNbrParticles], ".dat", TmpExtension);
      File << " " << TmpUpLayerVectorFileName;
      sprintf(TmpExtension, "_%d.%d.vec", DownLayerLzValue, TwoLayerDownIndices[TwoLayerIndices[i]]);
      char* TmpDownLayerVectorFileName = ReplaceExtensionToFileName(SpectrumFileNames[UpLayerNbrParticles], ".dat", TmpExtension);
      File << " " << TmpUpLayerVectorFileName<< " " << TmpDownLayerVectorFileName;
      File << endl;
      if (LargestUsedEigenstate[UpLayerNbrParticles][(TwoLayerUpLzValues[TwoLayerIndices[i]] - MinLzValues[UpLayerNbrParticles]) / 2] < TwoLayerUpIndices[TwoLayerIndices[i]])
	{
	  LargestUsedEigenstate[UpLayerNbrParticles][(TwoLayerUpLzValues[TwoLayerIndices[i]] - MinLzValues[UpLayerNbrParticles]) / 2] = TwoLayerUpIndices[TwoLayerIndices[i]];
	}
      if (LargestUsedEigenstate[DownLayerNbrParticles][(DownLayerLzValue - MinLzValues[DownLayerNbrParticles]) / 2] < TwoLayerDownIndices[TwoLayerIndices[i]])
	{
	  LargestUsedEigenstate[DownLayerNbrParticles][(DownLayerLzValue - MinLzValues[DownLayerNbrParticles]) / 2] = TwoLayerDownIndices[TwoLayerIndices[i]];
	}
      delete[] TmpUpLayerVectorFileName;
      delete[] TmpDownLayerVectorFileName;
    }
  File.close();

  for (int i = 0; i <= MaxNbrParticlesPerLayer; ++i)
    {
      for (int j = 0; j < NbrLzSectors[i]; ++j)
	{
	  if (LargestUsedEigenstate[i][j] >= 0)
	    {
	      cout << "require " << (LargestUsedEigenstate[i][j] + 1) << " eigenstates for N=" << i << " and Ky=" << ((j * 2) + MinLzValues[i]) << endl;
	    }
	}
    }

  if (Manager.GetBoolean("build-eigenstates"))
    {
      char* FilePrefix = new char[512];
      
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
      QuasiholeOnSphereWithSpinAndPairing* InputSpace;
      InputSpace = new QuasiholeOnSphereWithSpinAndPairing (KValue, RValue, TwoLayerLzSector, NbrFluxQuanta, TwoLayerSzSector, 
							    Manager.GetString("directory"), FilePrefix, true);
      for (int i = 0; i < EffectiveSubspaceDimension; ++i)
	{
	  int UpLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] + TwoLayerSzSector) / 2;
	  int DownLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] - TwoLayerSzSector) / 2;  
	  int DownLayerLzValue = TwoLayerUpLzValues[TwoLayerIndices[i]] - TwoLayerLzSector;
	  sprintf(TmpExtension, "_%d.%d.vec", TwoLayerUpLzValues[TwoLayerIndices[i]], TwoLayerUpIndices[TwoLayerIndices[i]]);
	  char* TmpUpLayerVectorFileName = ReplaceExtensionToFileName(SpectrumFileNames[UpLayerNbrParticles], ".dat", TmpExtension);
	  RealVector TmpUpLayerVector;
	  if (TmpUpLayerVector.ReadVector(TmpUpLayerVectorFileName) == false)
	    {
	      cout << "error, can't read " << TmpUpLayerVectorFileName << endl;
	  return -1;
	    }
	  sprintf(TmpExtension, "_%d.%d.vec", DownLayerLzValue, TwoLayerDownIndices[TwoLayerIndices[i]]);
	  char* TmpDownLayerVectorFileName = ReplaceExtensionToFileName(SpectrumFileNames[UpLayerNbrParticles], ".dat", TmpExtension);
	  RealVector TmpDownLayerVector;
	  if (TmpDownLayerVector.ReadVector(TmpDownLayerVectorFileName) == false)
	    {
	      cout << "error, can't read " << TmpUpLayerVectorFileName << endl;
	      return -1;
	    }
	  char* OutputVectorFileNameExtension = new char[128];
	  sprintf (OutputVectorFileNameExtension, "_effective_%d_fixedn_%d_n_0_", EffectiveSubspaceDimension, (UpLayerNbrParticles + DownLayerNbrParticles));
	  char* OutputVectorFileName1 = ReplaceString(TmpOutputFileName, "_n_0_", OutputVectorFileNameExtension);
	  sprintf (OutputVectorFileNameExtension, "%d.vec", i);
	  char* OutputVectorFileName2 = ReplaceExtensionToFileName(OutputVectorFileName1, "dat", OutputVectorFileNameExtension);
	  RealVector TmpTwoLayerVector = InputSpace->BuildFromTwoSingleLayerEigenstates(TmpUpLayerVector, UpLayerNbrParticles, TwoLayerUpLzValues[TwoLayerIndices[i]], 
											TmpDownLayerVector, DownLayerNbrParticles, DownLayerLzValue);
	  if (TmpTwoLayerVector.WriteVector(OutputVectorFileName2) == false)
	    {
	      cout << "error, can't write " << OutputVectorFileName2 << endl;
	      return -1;
	    }	  
	  cout << "generating eignestate " << OutputVectorFileName2 << endl;
	  delete[] OutputVectorFileName1;
	  delete[] OutputVectorFileName2;
	  delete[] TmpUpLayerVectorFileName;
	  delete[] TmpDownLayerVectorFileName;
	}
    }
      
  delete[] TwoLayerEnergies;
  delete[] TwoLayerNbrParticles;
  delete[] TwoLayerIndices;
  delete[] TwoLayerUpIndices;
  delete[] TwoLayerDownIndices;
  return 0;
}

// extract the spectrum from a file
//
// spectrumFile = spectrum file name
// nbrLzSectors = reference on the number of momentum sectors
// minLzSector = reference on the minimal momentum
// nbrEnergies =  reference on the array where the number of levels per momentum sector
// energies =  reference on the array where the emergies will be stored
// return value = true if no error occured

bool FQHECylinderQuasiholesWithSpinTimeReversalSymmetryExtractSpectrum(char* spectrumFile, int& nbrLzSectors, int& minLzSector,
								       int*& nbrEnergies, double**& energies)
{
  MultiColumnASCIIFile SpectrumFile;
  if (SpectrumFile.Parse(spectrumFile) == false)
    {
      SpectrumFile.DumpErrors(cout);
      return false;
    }
  if (SpectrumFile.GetNbrColumns() < 2)
    {
      cout << "too few columns in " << spectrumFile << endl;
      return false;
    }
  int TmpNbrValues = SpectrumFile.GetNbrLines();
  int* TmpLzSectors = SpectrumFile.GetAsIntegerArray(0);
  double* TmpEnergies = SpectrumFile.GetAsDoubleArray(1);
  minLzSector = TmpLzSectors[0];
  int CurrentLzValue = minLzSector;
  nbrLzSectors = 1;
  for (int i = 1; i < TmpNbrValues; ++i)
    {
      if (CurrentLzValue != TmpLzSectors[i])
	{
	  CurrentLzValue = TmpLzSectors[i];
	  ++nbrLzSectors;
	  if (minLzSector > CurrentLzValue)
	    {
	      minLzSector = CurrentLzValue;
	    }	  
	}
    }
  nbrEnergies = new int [nbrLzSectors];
  energies = new double* [nbrLzSectors];
  for (int i = 0; i < nbrLzSectors; ++i)
    {
      nbrEnergies[i] = 0;
    }
  for (int i = 0; i < TmpNbrValues; ++i)
    {
      ++nbrEnergies[(TmpLzSectors[i] - minLzSector) / 2];
    }
  for (int i = 0; i < nbrLzSectors; ++i)
    {
      energies[i] = new double[nbrEnergies[i]];
      nbrEnergies[i] = 0;
    }
  for (int i = 0; i < TmpNbrValues; ++i)
    {
      energies[(TmpLzSectors[i] - minLzSector) / 2][nbrEnergies[(TmpLzSectors[i] - minLzSector) / 2]] = TmpEnergies[i];
      ++nbrEnergies[(TmpLzSectors[i] - minLzSector) / 2];
    }
  delete[] TmpEnergies;
  delete[] TmpLzSectors;
  return true;
}
