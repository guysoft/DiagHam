
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//        class for periodic average spectra with Fourier-Bessel basis        //
//                                                                            //
//                        last modification : 07/05/2004                      //
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


#include "Tools/QuantumDot/Spectra/CylinderQuantumDotSpectra.h"
#include "MathTools/BesselJZeros.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;


// constructor from a Hilbert space and a file
//
// space = Hilbert space describing the particle
// fileName = name of the state file
// bz = magnetic field in Z direction
CylinderQuantumDotSpectra::CylinderQuantumDotSpectra(VerticalPeriodicParticleInMagneticField* space, char* fileName, double bz)
{
  this->NumberM = space->GetQuantumNumberM();
  this->NbrStateR = space->GetNbrStateR();
  this->NbrStateZ = space->GetNbrStateZ();
  this->LowerImpulsionZ = space->GetLowerImpulsionZ();
  this->Bz = bz;

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }

  this->RealCoefficients = new double* [this->NbrStateR];
  this->ImaginaryCoefficients = new double* [this->NbrStateR];
  for (int i = 0; i < this->NbrStateR; ++i)
    {
      this->RealCoefficients[i] = new double [this->NbrStateZ];
      this->ImaginaryCoefficients[i] = new double [this->NbrStateZ];         
      for (int k = 0; k < this->NbrStateZ; ++k)		
	File >> this->RealCoefficients[i][k] >> this->ImaginaryCoefficients[i][k];	    	
    } 
  File.close();
}

// get the value of impulsion operators with another wavefunction <this|p|another>
//
// space = Hilbert space describing the other particle
// fileName = the file to stock the other function
// sizeZ = size of sample in Z direction
// sizeR = size of the super-cylinder in plane
// impulsionX, impulsionY, impulsionZ = reference to the return values

void CylinderQuantumDotSpectra::GetImpulsion(VerticalPeriodicParticleInMagneticField* space, char* fileName, double sizeZ, double sizeR, double &realImpulsionX, double &imaginaryImpulsionX, double &realImpulsionY, double &imaginaryImpulsionY, double &realImpulsionZ, double &imaginaryImpulsionZ)
{
  int numberM = space->GetQuantumNumberM();
  int nbrStateR = space->GetNbrStateR();
  int nbrStateZ = space->GetNbrStateZ();
  int lowerImpulsionZ = space->GetLowerImpulsionZ();

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }
  double** realCoefficients = new double* [nbrStateR];
  double** imaginaryCoefficients = new double* [nbrStateR];
  for (int i = 0; i < nbrStateR; ++i)
    {
      realCoefficients[i] = new double [nbrStateZ];
      imaginaryCoefficients[i] = new double [nbrStateZ];         
      for (int k = 0; k < nbrStateZ; ++k)		
	File >> realCoefficients[i][k] >> imaginaryCoefficients[i][k]; 	
    }
  File.close();  

  int MaxLowerZ = this->LowerImpulsionZ;
  if (this->LowerImpulsionZ < lowerImpulsionZ)
    MaxLowerZ = lowerImpulsionZ;

  int MinUpperR = this->NbrStateR, MinUpperZ = this->LowerImpulsionZ + this->NbrStateZ;
  if (MinUpperR > (nbrStateR))
    MinUpperR = nbrStateR;
  if (MinUpperZ > (lowerImpulsionZ + nbrStateZ))
    MinUpperZ = lowerImpulsionZ + nbrStateZ;

  realImpulsionX = 0.0;
  realImpulsionY = 0.0;
  realImpulsionZ = 0.0;
  imaginaryImpulsionX = 0.0;
  imaginaryImpulsionY = 0.0;
  imaginaryImpulsionZ = 0.0;
  double TmpRe = 0.0, TmpIm = 0.0;
  double* Re1; double* Im1; double* Re2; double* Im2;
  double TmpRe1 = 0.0, TmpIm1 = 0.0;
  // x and y directions
  if (this->NumberM == numberM)
    {
      realImpulsionX = 0.0;
      imaginaryImpulsionX = 0.0;
      realImpulsionY = 0.0;
      imaginaryImpulsionY = 0.0;
    }
  else
    {
      double* InverseZeros1 = new double [MinUpperR];
      double* InverseZeros2 = new double [MinUpperR];
      for (int n = 0; n < MinUpperR; ++n)
	{
	  //InverseZeros1[n] = BesselJZeros[1][n] / BesselJZeros[0][n];
	  //InverseZeros2[n] = BesselJZeros[0][n] / BesselJZeros[1][n];
	}
      TmpRe = 0.0; TmpIm = 0.0;
      for (int n1 = 0; n1 < MinUpperR; ++n1)
	{
	  Re1 = this->RealCoefficients[n1];
	  Im1 = this->ImaginaryCoefficients[n1];
	  for (int n2 = 0; n2 < MinUpperR; ++n2)
	    {
	      Re2 = realCoefficients[n2];
	      Im2 = imaginaryCoefficients[n2];
	      TmpRe1 = 0.0; TmpIm1 = 0.0;
	      for (int p = MaxLowerZ; p < MinUpperZ; ++p) 
		{
		  TmpRe1 += (Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
		  TmpIm1 += (-Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);
		}
	      double fraction = BesselJZeros[0][n1] / BesselJZeros[1][n2];
	      TmpRe += (TmpRe1 / (fraction - 1.0 / fraction));
	      TmpIm += (TmpIm1 / (fraction - 1.0 / fraction));
	    }
	}
      realImpulsionX = TmpRe / sizeR; realImpulsionY = TmpRe / sizeR; 
      imaginaryImpulsionX = TmpIm / sizeR; imaginaryImpulsionY = TmpIm / sizeR; 
    }
  // z direction
  if (this->NumberM != numberM)
    {
      realImpulsionZ = 0.0;
      imaginaryImpulsionZ = 0.0;
    }
  else
    {          
      for (int n = 0; n < MinUpperR; ++n)
	{	      
	  Re1 = this->RealCoefficients[n];
	  Im1 = this->ImaginaryCoefficients[n];
	  Re2 = realCoefficients[n];
	  Im2 = imaginaryCoefficients[n];
	  for (int p = MaxLowerZ; p < MinUpperZ; ++p)
	    {
	      realImpulsionZ += ((Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]) * double(p));
	      imaginaryImpulsionZ += ((Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] - Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]) * double(p));
	    }
	}
      realImpulsionZ = realImpulsionZ * 2.0 * M_PI / sizeZ; imaginaryImpulsionZ = imaginaryImpulsionZ * 2.0 * M_PI / sizeZ;	    
    }  
}
