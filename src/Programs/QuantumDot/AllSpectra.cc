#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/QuantumDot/Spectra/Spectra.h"
#include "Tools/QuantumDot/Spectra/OverlapSpectra.h"
#include "Tools/QuantumDot/Spectra/AverageSpectra.h"
#include "Tools/QuantumDot/Spectra/DOSSpectra.h"
#include "Tools/QuantumDot/Spectra/PeriodicSpectra.h"
#include "Tools/QuantumDot/Spectra/XYReflexionSymmetricPeriodicSpectra.h"

#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/XYReflexionSymmetricPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXPairYPeriodic3DOneParticle.h"

#include "BitmapPicture/BmpFormat.h"


#include <iostream>
#include <fstream>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;

int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("AllSpectra" , "0.01");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* FileGroup =  new OptionGroup ("File and energy options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += HilbertSpaceGroup;
  Manager += FileGroup;
  Manager += MiscGroup;

  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statex", "number of states in x direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statey", "number of states in y direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statez", "number of states in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);
  (*HilbertSpaceGroup) += new BooleanOption ('\n', "pairX", "pair function in X direction", false);
  (*HilbertSpaceGroup) += new BooleanOption ('\n', "pairY", "pair function in Y direciton", false);

  (*FileGroup) += new SingleStringOption('\n', "file1", "name of the electron state file", "");
  (*FileGroup) += new SingleStringOption('\n', "file2", "name of the hole state file", "");  
  (*FileGroup) += new SingleDoubleOption('\n', "e1", "energy of the first particle", 0.0);
  (*FileGroup) += new SingleDoubleOption('\n', "e2", "energy of the second particle", 0.0);
  (*FileGroup) += new SingleDoubleOption('g', "gap", "energy of the normal gap (no strain)", 1.52);  
  (*FileGroup) += new SingleDoubleOption('X', "sizeX", "size of the sample (in Angstrom unit)", 500);
  (*FileGroup) += new SingleDoubleOption('Y', "sizeY", "size of the sample (in Angstrom unit)", 500);
  (*FileGroup) += new SingleDoubleOption('Z', "sizeZ", "size of the sample (in Angstrom unit)", 500);

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrStateX = ((SingleIntegerOption*) Manager["nbr-statex"])->GetInteger();
  int NbrStateY = ((SingleIntegerOption*) Manager["nbr-statey"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  bool PairX = ((BooleanOption*) Manager["pairX"])->GetBoolean();
  bool PairY = ((BooleanOption*) Manager["pairY"])->GetBoolean();

  char* FileName1 = ((SingleStringOption*) Manager["file1"])->GetString();
  char* FileName2 = ((SingleStringOption*) Manager["file2"])->GetString();
  double Energy1 = ((SingleDoubleOption*) Manager["e1"])->GetDouble();
  double Energy2 = ((SingleDoubleOption*) Manager["e2"])->GetDouble();
  double Gap = ((SingleDoubleOption*) Manager["gap"])->GetDouble();
  double SizeX = ((SingleDoubleOption*) Manager["sizeX"])->GetDouble();
  double SizeY = ((SingleDoubleOption*) Manager["sizeY"])->GetDouble();
  double SizeZ = ((SingleDoubleOption*) Manager["sizeZ"])->GetDouble();

   // define Hilbert space    
  XYReflexionSymmetricPeriodic3DOneParticle GeneralSpace(NbrStateX / 2, NbrStateY / 2, NbrStateZ, LowImpulsionZ);
  XYReflexionSymmetricPeriodic3DOneParticle* Space;
  if (PairX)
    if (PairY)
      Space = new PairXPairYPeriodic3DOneParticle(GeneralSpace);     
    else
      Space = new PairXImpairYPeriodic3DOneParticle(GeneralSpace); 
  else
     if (PairY)
      Space = new ImpairXPairYPeriodic3DOneParticle(GeneralSpace);     
    else
      Space = new ImpairXImpairYPeriodic3DOneParticle(GeneralSpace);

  XYReflexionSymmetricPeriodicSpectra spectra(Space, FileName1);  

  //void GetDerivedOverlap (XYReflexionSymmetricPeriodic3DOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &realOverlap, double &imaginaryOverlap, double &realOverlapX, double &imaginaryOverlapX, double &realOverlapY, double &imaginaryOverlapY);

  double real, imaginary, realX, imaginaryX, realY, imaginaryY;
  spectra.GetDerivedOverlap(Space, FileName2, SizeX, SizeY, SizeZ, real, imaginary, realX, imaginaryX, realY, imaginaryY);
  
  double re1, re2, im1, im2;
  
  re1 = real - 23 * 1.6 * (realX - realY)/ ((Gap + Energy1 + Energy2) * (Gap + Energy1 + Energy2));
  re2 = real + 23 * 1.6 * (realX - realY)/ ((Gap + Energy1 + Energy2) * (Gap + Energy1 + Energy2));
  im1 = imaginary -  23 * 1.6 * (imaginaryX - imaginaryY) / ((Gap + Energy1 + Energy2) * (Gap + Energy1 + Energy2));
  im2 = imaginary +  23 * 1.6 * (imaginaryX - imaginaryY) / ((Gap + Energy1 + Energy2) * (Gap + Energy1 + Energy2));

  double tmp1 = re1 * re1 + im1 * im1;
  double tmp2 = re2 * re2 + im2 * im2;

  cout << "Polarization degree is: " << ((tmp1 - tmp2) / (tmp1 + tmp2)) << endl;

  XYReflexionSymmetricPeriodicSpectra spectra2(Space, FileName2);
  int Number = 50;
  ofstream Electron ("Function_Electron.txt");
  ofstream Hole ("Function_Hole.txt");
  double density = 0.0;
  for (int j = 0; j <= Number; ++j)
    {
      for (int i = 0; i <= Number; ++i)
	{
	  density = spectra.PlanarProbabilityDensity(i * SizeX / Number, SizeX, j * SizeY / Number, SizeY);
	  Electron << density << " ";
	  density = spectra2.PlanarProbabilityDensity(i * SizeX / Number, SizeX, j * SizeY / Number, SizeY);
	  Hole << density << " ";
	}
      Electron << '\n'; Hole << '\n';
    }
  Electron.close(); Hole.close();
  
  /*
  // some running options and help 
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFile('\n', "input", "name of the input file", 0);
  SingleIntegerOption XCell('X', "Xcells", "number of the cells in X direction", 161);
  SingleIntegerOption YCell('Y', "Ycells", "number of the cells in Y direction", 161);
  SingleIntegerOption ZCell('Z', "Zcells", "number of the cells in Z direction", 21);
  SingleIntegerOption Division('D', "division", "number of the sub-divisions in each direction", 5);
  SingleStringOption Output('r', "output", "name of the output file", "default_output.txt");
  BooleanOption PairXOption ('\n', "pairX", "pair function in X direction", false);
  BooleanOption PairYOption ('\n', "pairY", "pair function in Y direciton", false);
  BooleanOption PairX2Option ('\n', "pairX2", "pair function in X direction of other particles", false);
  BooleanOption PairY2Option ('\n', "pairY2", "pair function in Y direciton of other particles", false);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &InputFile;
  OptionList += &XCell;
  OptionList += &YCell;
  OptionList += &ZCell;
  OptionList += &Division;
  OptionList += &Output;
  OptionList += &PairXOption;
  OptionList += &PairYOption;
  OptionList += &PairX2Option;
  OptionList += &PairY2Option;

  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "bad options" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }

  int M, N, H;
  char* FileName;
  int Number;
  
  FileName = InputFile.GetString();
  M = XCell.GetInteger();
  N = YCell.GetInteger();
  H = ZCell.GetInteger();
  Number = Division.GetInteger();
  char * out = Output.GetString();
  bool PairX = PairXOption.GetBoolean();
  bool PairY = PairYOption.GetBoolean();
  bool PairX2 = PairX2Option.GetBoolean();
  bool PairY2 = PairY2Option.GetBoolean();
 
  
  XYReflexionSymmetricPeriodic3DOneParticle GeneralSpace(40, 40, 21, -10);
  XYReflexionSymmetricPeriodic3DOneParticle* Space;
  if (PairX)
    if (PairY)
      Space = new PairXPairYPeriodic3DOneParticle(GeneralSpace);
    else
      Space = new PairXImpairYPeriodic3DOneParticle(GeneralSpace);
  else
     if (PairY)
      Space = new ImpairXPairYPeriodic3DOneParticle(GeneralSpace);
    else
      Space = new ImpairXImpairYPeriodic3DOneParticle(GeneralSpace);

  XYReflexionSymmetricPeriodic3DOneParticle* Space2;
  if (PairX2)
    if (PairY2)
      Space2 = new PairXPairYPeriodic3DOneParticle(GeneralSpace);
    else
      Space2 = new PairXImpairYPeriodic3DOneParticle(GeneralSpace);
  else
     if (PairY2)
      Space2 = new ImpairXPairYPeriodic3DOneParticle(GeneralSpace);
    else
      Space2 = new ImpairXImpairYPeriodic3DOneParticle(GeneralSpace);
  XYReflexionSymmetricPeriodicSpectra spectra(Space, FileName); 
  
  //Periodic3DOneParticle* Space = new Periodic3DOneParticle(M / 2 + 1 , -M / 4, N / 2 + 1, -N / 4, H, -H / 2);
  //PeriodicSpectra spectra(Space, FileName);
   
  double Lx = 5.65, Ly = 5.65, Lz = 5;
  double SizeX = M * Lx, SizeY = N * Ly, SizeZ = H * Lz;
  
  ofstream OutFile(out);    
  double PositionX = (double(M) * Lx / 2.0);
  double RealValue = 0.0, ImaginaryValue = 0.0;
  for (int j = 0; j <= N; ++j)
    {      
      for (int k = 0; k <= H; ++k)
  	{
  	  spectra.WaveFunctionValue(PositionX, SizeX, j * Ly, SizeY, k * Lz, SizeZ, RealValue, ImaginaryValue);
  	  OutFile << (RealValue * RealValue + ImaginaryValue * ImaginaryValue) << " "; 
  	}
      OutFile << '\n';
    }
  OutFile.close(); 
  */
  /*
  double ReX, ImX, ReY, ImY, ReZ, ImZ;
  char** Files = new char* [200];
  ifstream energy("eigenvalues");
  double fundamental;
  energy >> fundamental;
  double tmpE;  
  ofstream polarization ("Polarization.txt");

  
  ofstream PX("PolarizationX.txt");
  ofstream PY("PolarizationY.txt");  
  ofstream PZ("PolarizationZ.txt");
  
  for (int i = 1; i < 100; ++i)
    {
      Files[i] = new char[80];
      AddString(Files[i], "eigenvector.", i, "");
      spectra.GetImpulsion(Space2, Files[i], SizeX, SizeY, SizeZ, ReX, ImX, ReY, ImY, ReZ, ImZ);
      energy >> tmpE;
      polarization << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) / (tmpE - fundamental) << '\t' << ((ReY * ReY) + (ImY * ImY)) / (tmpE - fundamental) << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental) << endl;
      PX << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) / (tmpE - fundamental) << endl;
      PY << tmpE - fundamental << '\t' << ((ReY * ReY) + (ImY * ImY)) / (tmpE - fundamental) << endl;
      PZ << tmpE - fundamental << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental) << endl;
      cout << i << endl;
    }
  
  double EigenZ[] = {0.013066932129562, 0.1236666401825, 0.26229073118036};
  fundamental = -0.13225305464614; 
  for (int p = 0; p < 3; ++p)
    {
      Files[p] = new char[80];
      AddString(Files[p], "eigenvector.", p, "");
      for (int m = 0; m < 20; ++m)
	for (int n = 0; n < 20; ++n)
	  {
	    tmpE = EigenZ[p] + 150.4 * (m * m + n * n) / (0.07 * SizeX * SizeX);
	    if (tmpE < 0.4)
	      {	     
		spectra.GetImpulsionWithContinuum(m, n, 21, -10, Files[p], SizeZ, ReZ, ImZ);
		PZ << tmpE - fundamental << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental) << endl;	    
	      }
	  }      
    }
  
  PX.close(); PY.close(); PZ.close();
  energy.close(); polarization.close();
  */
  /*
  // void PeriodicSpectra::DensityProbability(double x, double SizeX, double y, double SizeY, double z, double SizeZ, double& Real, double& Imaginary)
  double PositionX = (double(M) * Lx / 2.0);
  double RealValue = 0.0, ImaginaryValue = 0.0;  
  for (int j = 0; j <= N; ++j)
    {      
      for (int k = 0; k <= H; ++k)
  	{
  	  spectra.WaveFunctionValue(PositionX, SizeX, j * Ly, SizeY, k * Lz, SizeZ, RealValue, ImaginaryValue);
  	  OutFile << (RealValue * RealValue + ImaginaryValue * ImaginaryValue) << " "; 
  	}
      OutFile << '\n';
    }
  OutFile.close();  
  
  double square = 0.0; double moyenne = 0.0;
  
  for (int i = 1; i < 150; ++i)
    {
      Files[i] = new char[80];
      AddString(Files[i], "eigenvector.", i, "");
      PeriodicSpectra* spectra = new PeriodicSpectra(Space, Files[i]);
      spectra->GetMeanValueZ(square);
      energy >> tmpE;
      mean << tmpE - fundamental << '\t' << moyenne << '\t' << square << '\n';
      delete spectra;
    }
  */

  /*
  Periodic3DOneParticle* Space = new Periodic3DOneParticle(M, M / 2, N, N / 2, H / 2 + 1, H / 4);
  
  PeriodicSpectra spectra(Space, FileName);
  ofstream OutFile(out);  

  // void PeriodicSpectra::DensityProbability(double x, double SizeX, double y, double SizeY, double z, double SizeZ, double& Real, double& Imaginary)

  double Lx = 2.97, Ly = 2.97, Lz = 2.64;
  double SizeX = M * Lx, SizeY = N * Ly, SizeZ = H * Lz;
  double PositionX = (double(M) * Lx / 2.0);
  double RealValue = 0.0, ImaginaryValue = 0.0;

  for (int j = 0; j <= N; ++j)
    {      
      for (int k = 0; k <= H; ++k)
	{
	  spectra.WaveFunctionValue(PositionX, SizeX, j * Ly, SizeY, k * Lz, SizeZ, RealValue, ImaginaryValue);
	  OutFile << (RealValue * RealValue + ImaginaryValue * ImaginaryValue) << " "; 
	}
      OutFile << '\n';
    }
  OutFile.close();

  ofstream File;
  File.open("MeanH.txt", ios::binary | ios::out | ios::app);
  double squareX, squareY, squareZ;
  File << spectra.GetMeanValueX(squareX) << '\t';
  File << squareX << '\t';
  File << spectra.GetMeanValueY(squareY) << '\t';
  File << squareY << '\t';
  File << spectra.GetMeanValueZ(squareZ) << '\t';
  File << squareZ << '\n';
  File.close();
  */
  /*
  ThreeDPotential potential(50, 50, 30, 6, 10);
  potential.ReadDiagram("Diagram/Diagram/0.175/h/Diagram.txt");
  PicRGB background (255, 255, 255);
  PicRGB InN (0, 0, 0);
  PicRGB GaN (231, 231, 231);
  // bool Potential::SaveBmpPicture(int under, int above, int startX, int endX, int startY, int endY, int choice, int sizeX, int sizeY, PicRGB& InN, PicRGB& GaN, PicRGB& background, int NbrX, char* fileName);
  potential.SaveBmpPicture(9, 20, 0, 50, 0, 50, 1, 5, 5, InN, GaN, background, 4, "Diagram/Diagram/0.175/h/Diagram.bmp");
  */
  /*
  char** Files = new char* [1]; int* State = new int[1];
  for (int i = 0; i < 1; ++i)
    {
      State[i] = 100;
      Files[i] = new char[80];
      Files[0] = FileName;
    }
  DOSSpectra DOS(1, Files, State, 4e-3, -0.2, 0.5, 2e-4);
  DOS.WriteSpectra(out);
  */

  /*
  for (int n = 102; n < 110; ++n)
    {  
      char* EE = new char [50]; char* ES = new char [50];
      char* HE = new char [50]; char* HS = new char [50];
      
      AddString(ES, "./", n, "/Etat_E.txt");
      AddString(HS, "./", n, "/Etat_H.txt");
      AddString(EE, "./", n, "/Energy_E.txt");
      AddString(HE, "./", n, "/Energy_H.txt");
      
      OverlapSpectra a(ES, EE, 9, HS, HE, 9, 14400);
      char* Over = new char [50];
      AddString(Over, "./", n, "/Overlap.txt");
      a.WriteSquareOverlap(Over);
      delete[] EE;  delete[] ES;  delete[] HE;  delete[] HS;
      EE = 0; ES = 0; HE = 0; HS = 0;
    }
*/
//  OverlapSpectra(char* ElectronStateFile, char* ElectronEnergyFile ,int ElectronNumber, char* HoleStateFile, char* HoleEnergyFile, int HoleNumber, int NumberState);
//  AverageSpectra(char* StateFile, char* EnergyFile ,int Number, int M, int N, double a, double b)
  //AverageSpectra(char** StateFile, char* EnergyFile ,int Number, int M, int N, int H, double a, double b, double c);
  /*
  char** a = new char*[2]; a[0] = new char[40];  a[1] = new char[40];
  a[0] = "eige.0"; a[1] = "eige.1";
  
  AverageSpectra b (a, "eiva", 2, 40, 40, 25, 2.97, 2.97, 2.64);
  b.WriteXYZ("mean.txt");
  */
/*
  for (int n = 60; n < 110; ++n)
    {
      char* State = new char [50];
      char* Energy = new char [50];
      AddString(State, "./", n, "/Etat_H.txt");
      AddString(Energy, "./", n, "/Energy_H.txt");
      AverageSpectra b(State, Energy, 9, 120, 120, 2.97, 2.97);
      char* File = new char [50];
      // AddString(File, "./", n, "/Mean_E.txt");
      //b.WriteXY(File);
      
      b.WriteVarMean("Dist_H.txt");

      delete[] State; delete[] Energy; delete[] File; 
      State = 0; Energy = 0; File = 0;
    }
*/
/*
  int Nbr = 1;
  char** Files = new char* [Nbr]; int* State = new int[Nbr];
  for (int i = 0; i < Nbr; ++i)
    {
      State[i] = 99;
      Files[i] = new char[80];
      Files[i] = FileName;
    }
  Spectra Absorption (Nbr, Files, State, 4e-3, 0.0, 0.6, 2e-4);
  Absorption.WriteSpectra(out);
*/

//Spectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE);

  return 0;
}
