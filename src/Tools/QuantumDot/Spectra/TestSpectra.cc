#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/QuantumDot/Spectra/Spectra.h"
#include "Tools/QuantumDot/Spectra/DOSSpectra.h"
#include "Tools/QuantumDot/Spectra/OverlapSpectra.h"
#include "Tools/QuantumDot/Spectra/AverageSpectra.h"
#include "Tools/QuantumDot/Spectra/TimeResolvedPLSpectra.h"
#include "Tools/QuantumDot/Spectra/CylinderInMagneticFieldSpectra.h"

#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/XYReflexionSymmetricPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/VerticalPeriodicParticleInMagneticField.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::ifstream;
using std::endl;
using std::ofstream;
using std::ios;

int main(int argc, char** argv)
{
  cout.precision(14);
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFile('\n', "input", "name of the input file", 0);
  SingleIntegerOption NumberRValueOption ('R', "R-states", "number of states in plane", 100);
  SingleIntegerOption NumberZValueOption ('Z', "Z-states", "number of cells in z direction", 100);
  SingleIntegerOption NumberMValueOption ('m', "momentum", "quantum number of kinetic in z direction", 0);
  SingleDoubleOption MagneticFieldOption ('b', "magnetic", "magnetic field in Z direction (in Tesla unit)", 30);
  SingleDoubleOption SizeZOption ('z', "size-z", "size of sample in Z direction (in Angstrom unit)", 118.65);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &InputFile;
  OptionList += &NumberRValueOption;
  OptionList += &NumberZValueOption;
  OptionList += &NumberMValueOption;
  OptionList += &MagneticFieldOption;
  OptionList += &SizeZOption;

  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type ExplicitMatrixExample -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }

  char* FileName = InputFile.GetString();
  int NbrStateR = NumberRValueOption.GetInteger();
  int NbrStateZ = NumberZValueOption.GetInteger();
  int NumberM = NumberMValueOption.GetInteger();
  double MagneticField = MagneticFieldOption.GetDouble();
  double SizeZ = SizeZOption.GetDouble();
  
  VerticalPeriodicParticleInMagneticField* space = new VerticalPeriodicParticleInMagneticField (0, NbrStateR, NbrStateZ, -NbrStateZ / 2);
  VerticalPeriodicParticleInMagneticField* space2 = new VerticalPeriodicParticleInMagneticField (NumberM, NbrStateR, NbrStateZ, -NbrStateZ / 2);
  CylinderInMagneticFieldSpectra* spectra = new CylinderInMagneticFieldSpectra(space, FileName, MagneticField);
  
  double delta = SizeZ / 100.0; double p = 0.0;
  double shift = 0.0;
  for (double z = shift; z <= (SizeZ + shift); z += delta)
    {
      p = spectra->ZProbabilityDensity(z, SizeZ);
      cout << (z - shift) << '\t' << p << '\n';
    }
  
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
      spectra->GetImpulsion(space2, Files[i], SizeZ, ReX, ImX, ReY, ImY, ReZ, ImZ);
      energy >> tmpE;
      polarization << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) / (tmpE - fundamental) << '\t' << ((ReY * ReY) + (ImY * ImY)) / (tmpE - fundamental) << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental) << endl;
      PX << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) / (tmpE - fundamental)  << endl;
      PY << tmpE - fundamental << '\t' << ((ReY * ReY) + (ImY * ImY)) / (tmpE - fundamental)  << endl;
      PZ << tmpE - fundamental << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental)  << endl;
      //cout << i << endl;
    }
  PX.close(); PY.close(); PZ.close();
  energy.close(); polarization.close();
  */

/*
  //DOSSpectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE)
  char** name = new char*[1];
  name[0] = "Sinus50.txt";
  int state[1]; state[0] = 10;
  DOSSpectra spectra(1, name, state, 5e-4, -0.011, 0.18, 1e-4);
  spectra.WriteSpectra("Sinus50");
*/
/*
  // OverlapSpectra(char* ElectronStateFile, char* ElectronEnergyFile ,int ElectronNumber, char* HoleStateFile, char* HoleEnergyFile, int HoleNumber);
  //  OverlapSpectra x (36);
  //  Spectra (5);
  // e.WriteSpectra("yes");
  aa->WriteVector("hehe");
  bb->WriteVector("hehe");
  ifstream te("hehe");
  RealVector c(6), d(3); 
  te >> c >> d;
  c.WriteAsciiVector("he1");
  d.WriteAsciiVector("he2");
  te.close();
  */
/*
  OverlapSpectra x("Etat_E.txt", "Energy_E.txt", 9, "Etat_H.txt", "Energy_H.txt", 9, 14400);
  //cout << "hehe" << "\n";
  //OverlapSpectra x("ES", "EE", 1, "HS", "HE", 1, 4);
  //x.WriteSpectra("Overlap1.txt");
  x.WriteSquareOverlap("Overlap.txt");
*/  

//  AbsorptionSpectra::AbsorptionSpectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE)
/*
  char** f = new char*[3]; f[0] = "1"; f[1] = "2"; f[2] = "3";
  int* s = new int[3]; s[0] = 4; s[1] = 5; s[2] = 5;
  AbsorptionSpectra y(3, f, s, 0.2, -2, 6, 0.1);
  y.WriteSpectra("Absorption.txt");
*/
  //AverageSpectra::AverageSpectra(char* StateFile, char* EnergyFile ,int Number, int M, int N, double a, double b)
  //AverageSpectra z("Etat_E.txt", "Energy_E.txt", 9, 120, 120, 2.97, 2.97);
  //z.WriteXY("Mean_E.txt");
  /*
  char* s = new char[80];
  AddString(s, "./", 4, "eh.txt");
  cout << s << '\n';
  */
  //TimeResolvedPLSpectra(int FileNumber, char** FileName, int* ElectronState, int* HoleState, double Emin, double Emax, double T, int NT);
  /* 
     char** f = new char*[50];
     int* es = new int[50];
     int* hs = new int[50];
     
     for (int i = 0; i < 50; ++i)
     {
     f[i]  = new char[60];
     AddString(f[i], "./", i + 60, "/Overlap.txt");
     es[i] = 9;
     hs[i] = 9;
     }
     TimeResolvedPLSpectra a(50, f, es, hs, -0.02, -0.01, 20, 20000);
     a.WriteSpectra("PLTest.txt");
  */
  return 0;
}
