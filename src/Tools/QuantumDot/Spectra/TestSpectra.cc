#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
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
#include "Tools/QuantumDot/Spectra/CylinderQuantumDotSpectra.h"

#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/XYReflexionSymmetricPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/VerticalPeriodicParticleInMagneticField.h"

#include "Tools/QuantumDot/Potential/QuantumDotThreeDConstantCylinderPotential.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::ifstream;
using std::endl;
using std::ofstream;
using std::ios;

int main(int argc, char** argv)
{
  /*
  cout.precision(14);
  OptionManager Manager ("CylinderQuantumDotInMagneticField" , "0.01");
  OptionGroup* PotentialGroup = new OptionGroup ("potential options");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += PotentialGroup;
  Manager += HilbertSpaceGroup;
  Manager += LanczosGroup;
  Manager += MiscGroup;

  (*PotentialGroup) += new SingleDoubleOption ('\n', "radius", "radius of the supercylinder (in Angstrom unit)", 1000);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "barrier", "number of cells in the well barrier", 10.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "below", "width of the layer below the wetting layer (in Angstrom unit)", 10.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "wetting", "width of the wetting layer (in Angstrom unit)", 5.0);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "nbr-dot", "number of uniformly high layer in the dot", 3);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "base", "base radius in Angstrom unit", 100.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "height", "height of dot in Angstrom unit", 17.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "top", "top radius in Anstrom unit", 74.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "above", "width of the layer above the dot layer (in Angstrom unit)", 70.0);

  (*HilbertSpaceGroup) += new SingleIntegerOption ('R', "R-states", "number of states in plane", 50);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('Z', "Z-states", "number of cells in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('m', "momentum", "quantum number of kinetic in z direction", 0);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('k', "wave", "wave vector of Bloch function in Z direction (in 1/Angstrom unit)", 0.0);

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");
  (*MiscGroup) += new SingleStringOption ('\n', "input", "file input", ""); 

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

  double SuperCylinderRadius = ((SingleDoubleOption*) Manager["radius"])->GetDouble();
  double Barrier = ((SingleDoubleOption*) Manager["barrier"])->GetDouble();
  double Below = ((SingleDoubleOption*) Manager["below"])->GetDouble();
  double WettingWidth = ((SingleDoubleOption*) Manager["wetting"])->GetDouble();
  double BaseRadius = ((SingleDoubleOption*) Manager["base"])->GetDouble();
  double DotHeight = ((SingleDoubleOption*) Manager["height"])->GetDouble();
  int DotNbr = ((SingleIntegerOption*) Manager["nbr-dot"])->GetInteger();
  double TopRadius = ((SingleDoubleOption*) Manager["top"])->GetDouble();
  double Above = ((SingleDoubleOption*) Manager["above"])->GetDouble();

  int NbrStateR = ((SingleIntegerOption*) Manager["R-states"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["Z-states"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  int NumberM = ((SingleIntegerOption*) Manager["momentum"])->GetInteger();
  double WaveVector = ((SingleDoubleOption*) Manager["wave"])->GetDouble();

  char* FileName = ((SingleStringOption*) Manager["input"])->GetString();

  QuantumDotThreeDConstantCylinderPotential* potential = new QuantumDotThreeDConstantCylinderPotential(Below, WettingWidth, DotNbr, DotHeight, BaseRadius, TopRadius, Above, Barrier, SuperCylinderRadius);

  VerticalPeriodicParticleInMagneticField* Space = new VerticalPeriodicParticleInMagneticField(NumberM, NbrStateR, NbrStateZ, LowImpulsionZ); 

  CylinderQuantumDotSpectra* spectra = new CylinderQuantumDotSpectra(Space, FileName, 0.0);
  cout << "The probability to find the particle in the dot is: " << spectra->GetDotProbability(potential) << endl;
  */

  cout.precision(14);
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFile('\n', "input", "name of the input file", 0);
  SingleIntegerOption NumberRValueOption ('R', "R-states", "number of states in plane", 100);
  SingleIntegerOption NumberZValueOption ('Z', "Z-states", "number of states in z direction", 100);
  SingleIntegerOption LowZOption ('\n', "lowz", "lower impulsion in z direction", -10);
  SingleIntegerOption NumberMValueOption ('m', "momentum", "quantum number of kinetic in z direction", 0);
  SingleDoubleOption MagneticFieldOption ('b', "magnetic", "magnetic field in Z direction (in Tesla unit)", 30);
  SingleDoubleOption SizeZOption ('z', "size-z", "size of sample in Z direction (in Angstrom unit)", 118.65);
  SingleDoubleOption SizeROption ('r', "size-r", "size of sample in plane (in Angstrom unit)", 1000);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &InputFile;
  OptionList += &NumberRValueOption;
  OptionList += &NumberZValueOption;
  OptionList += &LowZOption;
  OptionList += &NumberMValueOption;
  OptionList += &MagneticFieldOption;
  OptionList += &SizeZOption;
  OptionList += &SizeROption;

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
  int LowZ = LowZOption.GetInteger();
  int NumberM = NumberMValueOption.GetInteger();
  double MagneticField = MagneticFieldOption.GetDouble();
  double SizeZ = SizeZOption.GetDouble();
  double SizeR = SizeROption.GetDouble();
  /*
  VerticalPeriodicParticleInMagneticField* space = new VerticalPeriodicParticleInMagneticField (0, NbrStateR, NbrStateZ, LowZ);
  VerticalPeriodicParticleInMagneticField* space2 = new VerticalPeriodicParticleInMagneticField (NumberM, NbrStateR, NbrStateZ, LowZ);
  CylinderQuantumDotSpectra* spectra = new CylinderQuantumDotSpectra(space, FileName, MagneticField);

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
  for (int i = 0; i < 4; ++i)
    {
      Files[i] = new char[80];
      AddString(Files[i], "eigenvector.", i, "");
      spectra->GetImpulsion(space2, Files[i], SizeZ, SizeR, ReX, ImX, ReY, ImY, ReZ, ImZ);
      energy >> tmpE;
      polarization << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) / (tmpE - fundamental) << '\t' << ((ReY * ReY) + (ImY * ImY)) / (tmpE - fundamental) << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental) << endl;
      PX << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) / (tmpE - fundamental)  << endl;
      PY << tmpE - fundamental << '\t' << ((ReY * ReY) + (ImY * ImY)) / (tmpE - fundamental)  << endl;
      PZ << tmpE - fundamental << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental)  << endl;
      //cout << i << endl;
    } 
  */

  
  VerticalPeriodicParticleInMagneticField* space = new VerticalPeriodicParticleInMagneticField (0, NbrStateR, NbrStateZ, LowZ);
  VerticalPeriodicParticleInMagneticField* space2 = new VerticalPeriodicParticleInMagneticField (NumberM, NbrStateR, NbrStateZ, LowZ);
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
  char** Files = new char* [NbrStateR];
  int* StateNumber = new int [NbrStateR];
  for (int i = 0; i < NbrStateR; ++i)
    {
      StateNumber[i] = 100;
      Files[i] = new char [80];
      if (i < 9)
	AddString(Files[i], "./0.00", i + 1, "/eigenvalues");
      else
	AddString(Files[i], "./0.0", i + 1, "/eigenvalues");
    }
  DOSSpectra spectra(NbrStateR, Files, StateNumber, 4e-3, -0.20, 0.5, 1e-4);
  spectra.WriteSpectra(FileName);
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
