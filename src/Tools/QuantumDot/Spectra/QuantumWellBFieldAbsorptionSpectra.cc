////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class for absorption spectra of quantum well            //
//                                in magnetic field                           //
//                                                                            //
//                        last modification : 06/12/2005                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"

#include "Tools/QuantumDot/Spectra/QuantumWellBFieldAbsorptionSpectra.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include <fstream>
#include <math.h>
#include <stdlib.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;


// constructor from a set of energy files. Each peak is assimilated to a Lorentzian function.
//
// fileNumber=  number of files
// files = name of files
// stateNumber = integer array containing number of states in each file
// thetaPolarizationAngle = angle between the z axis and the polarization vector
// phiPolarizationAngle = angle between the x axis and the projection of the polarization vector onto the x-y plane
// zSize = system dimension in the z direction (in Angstrom unit)
// gamma = lorentzian broadening parameter
// beta = beta factor (inverse of the temperature in energy unit 1/kT) 
// eMin = photon minimum energy
// eMax = photon maximum energy
// deltaE = photon energy step

QuantumWellBFieldAbsorptionSpectra::QuantumWellBFieldAbsorptionSpectra(int fileNumber, int nbrInitialStates, char** initialStateSpectrumFiles, char** initialStateEigenstateFiles, 	  
								       int nbrFinalStates, char** finalStateSpectrumFiles, char** finalStateEigenstateFiles,
								       double thetaPolarizationAngle, double phiPolarizationAngle, double zSize, 
								       double gamma, double beta, double eMin, double eMax, double deltaE)
{
  int N = (int) ((eMax - eMin) / deltaE);
  double* Energy = new double [N];
  double* Absorption = new double [N]; 
  double tmp1 = eMin; 
  double tmp2 = 0.0; 
  double g = Gamma * Gamma * 0.25;  

  for (int i = 0; i < N; ++i)
    {
      Energy[i] = tmp1;
      Absorption[i] = 0.0;
      tmp1 += deltaE;
    }

  this->ComputeOscillatorStrengthMatrix(thetaPolarizationAngle, phiPolarizationAngle, zSize);

  for (int i = 0; i < fileNumber; ++i)
    {
      int TmpNbrStates = StateNumber[i];
      double* tmp = new double [TmpNbrStates];
      ifstream file;
      file.open(Files[i],ios::out);
      if (!file.is_open())
        {
	  cout << "Error in open the file: " << Files[i] << "Exit now" << endl;
	  exit(0);
	}
      for (int j = 0; j < TmpNbrStates; ++j)
	{
	  file >> tmp[j];
	}
      file.close();
      
      tmp1 = eMin; 
      for (int j = 0; j < N; ++j)
	{
	  tmp2 = 0.0;
	  for (int k = 0; k < nbrStates; ++k)
	    {
	      tmp3 = 0.0;
	      for (int l = k + 1; l < nbrStates; ++l)
		tmp3 += 1.0 / ((tmp1 - tmp[k]) * (tmp1 - tmp[k]) + g);
	      
	    }
	  Absorption[j] += tmp2;
	  tmp1 += deltaE;
	}
      tmp = 0;
      delete[] tmp;
    }
  double tmp3 = Gamma / (2.0 * M_PI * ((double) FileNumber));
  for (int i = 0; i < N; ++i)
    Absorption[i] *= tmp3;

  this->AxeX = new RealVector(Energy, N);
  this->AxeY = new RealVector(Absorption, N);
  this->PointNumber = N;
}


void QuantumWellBFieldAbsorptionSpectra::AddSort (int nbrInitialStates, char* initialStateSpectrumFileName, char** initialEigenstateFileNames,
						  int nbrFinalStates, char* finalStateSpectrumFileName, char** finalEigenstateFileNames)
{
  int NbrFinalStates = 2 * nbrInitialStates;
  double* InitialEnergies = new double [nbrInitialStates];
  double** FinalEnergies = new double [NbrFinalStates];

  ifstream File1;
  File1.open(initialStateSpectrumFileName, ios::out);
  if (!File1.is_open())
    {
      cout << "Error in open the file: " << initialStateSpectrumFileName << "Exit now" << endl;
      exit(0);
    }
  for (int j = 0; j < nbrInitialStates; ++j)
    {
      File1 >> InitialEnergies[j];
    }
  File1.close();
  ifstream File2;
  File2.open(finalStateSpectrumFileName, ios::out);
  if (!File2.is_open())
    {
      cout << "Error in open the file: " << finalStateSpectrumFileName << "Exit now" << endl;
      exit(0);
    }
  for (int j = 0; j < NbrFinalStates; ++j)
    {
      File2 >> FinalEnergies[j];
    }
  File2.close();
      
  ComplexVector TmpInitialVector(nbrInitialStates);
  ComplexVector TmpInitialVector2(nbrInitialStates);
  ComplexVector TmpFinalVector(NbrFinalStates);
  Complex Tmp;
  tmp = 0.0;
  for (int i = 0; i < NbrFinalStates; ++i)
    {
      TmpFinalVector.ReadVector(finalEigenstateFileNames[i]);   
      TmpInitialVector2.Multiply(this->OscillatorStrength, TmpFinalVector);
      for (int j = 0; j < nbrInitialStates; ++j)
	if (FinalEnergies[i] > InitialEnergies[j])
	  {
	    TmpInitialVector.ReadVector(initialEigenstateFileNames[j]);
	    Tmp = (TmpInitialVector * TmpInitialVector2) * exp (this->Beta * this->FinalEnergies[i]);
	    for (int k = 0; k < this->NbrValues; ++k)
	      {
		double TmpEnergy = this->Energies[k];
		Tmp;

		for 
		tmp += / (((this->FinalEnergies[i] - this->InitialEnergy[j] - TmpEnergy) * (this->FinalEnergies[i] - this->InitialEnergy[j] - TmpEnergy)) 
			  + (this->Gamma * this->Gamma));
	      }
	  tmp *= ;
	}
    }

  delete[] InitialEnergies;
  delete[] FinalEnergies;
}

// compute the oscillator strength matrix 
//
// thetaPolarizationAngle = angle between the z axis and the polarization vector
// phiPolarizationAngle = angle between the x axis and the projection of the polarization vector onto the x-y plane
// zSize = system dimension in the z direction (in Angstrom unit)

void ComputeOscillatorStrengthMatrix(double thetaPolarizationAngle, double phiPolarizationAngle, double zSize)
{
  double CosTheta = cos (thetaPolarizationAngle);
  this->OscillatorStrength = ComplexMatrix (this->NbrInitialStates, this->NbrFinalStates, true);
  int HalfNbrFinalStates = this->NbrFinalStates / 2;

  double ZFactor = cos (thetaPolarizationAngle) * 2.0 * zSize / (M_PI * M_PI);
  double Tmp = 0;
  int Diff = 2 - 1;
  int Sum = 2 + 1;
  if ((Diff & 1) != 0)
    Tmp += 2.0 / ((double) (Diff * Diff));
  if ((Sum & 1) != 0)
    Tmp -= 2.0 / ((double) (Sum * Sum));
  ZFactor *= Tmp;

  for (int m = 0; m < HalfNbrFinalStates; ++m)
    for (int n = 0; n < this->NbrInitialStates; ++n)
      {
	this->OscillatorStrength.SetMatrixElement(n, 2 * m, ZFactor);
      }
}
