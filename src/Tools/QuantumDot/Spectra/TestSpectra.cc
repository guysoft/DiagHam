#include "Tools/QuantumDot/Potential/ThreeDPotential.h"

#include "Tools/QuantumDot/Spectra/Spectra.h"
#include "Tools/QuantumDot/Spectra/DOSSpectra.h"
#include "Tools/QuantumDot/Spectra/OverlapSpectra.h"
#include "Tools/QuantumDot/Spectra/AverageSpectra.h"
#include "Tools/QuantumDot/Spectra/TimeResolvedPLSpectra.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::ifstream;
using std::endl;

int main()
{
  //DOSSpectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE)
  char** name = new char*[1];
  name[0] = "Sinus50.txt";
  int state[1]; state[0] = 10;
  DOSSpectra spectra(1, name, state, 5e-4, -0.011, 0.18, 1e-4);
  spectra.WriteSpectra("Sinus50");


  // OverlapSpectra(char* ElectronStateFile, char* ElectronEnergyFile ,int ElectronNumber, char* HoleStateFile, char* HoleEnergyFile, int HoleNumber);
  //  OverlapSpectra x (36);
  //  Spectra (5);
  // e.WriteSpectra("yes");
  /*
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
