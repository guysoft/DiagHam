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

#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"

#include "BitmapPicture/BmpFormat.h"

#include "Vector/RealVector.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;

int main(int argc, char** argv)
{

  // some running options and help 
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFile('\n', "input", "name of the input file", 0);
  SingleIntegerOption XCell('X', "Xcells", "number of the cells in X direction", 61);
  SingleIntegerOption YCell('Y', "Ycells", "number of the cells in Y direction", 61);
  SingleIntegerOption ZCell('Z', "Zcells", "number of the cells in Z direction", 60);
  SingleIntegerOption Division('D', "division", "number of the sub-divisions in each direction", 5);
  SingleStringOption Output('r', "output", "name of the output file", "default_output.txt");

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &InputFile;
  OptionList += &XCell;
  OptionList += &YCell;
  OptionList += &ZCell;
  OptionList += &Division;
  OptionList += &Output;

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

  //ofstream File;
  //File.open("MeanE.txt", ios::binary | ios::out | ios::app);
  double squareX, squareY, squareZ;
  cout << spectra.GetMeanValueX(squareX) << '\t';
  cout << squareX << '\t';
  cout << spectra.GetMeanValueY(squareY) << '\t';
  cout << squareY << '\t';
  cout << spectra.GetMeanValueZ(squareZ) << '\t';
  cout << squareZ << '\n';
  //File.close();

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
  char** Files = new char* [50]; int* State = new int[50];
  for (int i = 0; i < 50; ++i)
    {
      State[i] = 9;
      Files[i] = new char[80];
      AddString(Files[i], "./", i + 60, "/Energy_E.txt");
    }
  DOSSpectra DOS(50, Files, State, 1e-3, -0.02, 0.03, 2e-5);
  DOS.WriteSpectra("DOS_E.txt");
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
  char** Files = new char* [50]; int* State = new int[50];
  for (int i = 0; i < 50; ++i)
    {
      State[i] = 81;
      Files[i] = new char[80];
      AddString(Files[i], "./", i + 60, "/Overlap.txt");
    }
  Spectra Absorption (50, Files, State, 1e-3, -0.03, 0.03, 2e-5);
  Absorption.WriteSpectra("Absorption.txt");
*/    
//Spectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE);
/*
  ifstream File("./46/Etat_E.txt");
  ifstream File1("./46/Etat_E1.txt");  
  ofstream Out("./46/Etat_E2.txt");

  RealVector* a = new RealVector[9];
  RealVector* b = new RealVector[9];
  for (int i = 0; i < 9; ++i)
    {
      a[i].ReadVector(File);
      b[i] = RealVector(3600);
      for (int j = 0; j < 3600; ++j)
	{
	  File1 >> b[i][j];
	  Out << a[i][j] << " ";
	}
      Out << '\n';
    }
  cout << a[8][35] << endl;
  cout << b[8][35] << endl;
  File.close(); Out.close();
*/
  return 0;
}
