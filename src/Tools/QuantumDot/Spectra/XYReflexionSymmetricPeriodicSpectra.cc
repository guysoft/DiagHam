////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//       class for periodic average spectra with XY reflexion symmetry        //
//                                                                            //
//                        last modification : 04/04/2004                      //
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


#include "Tools/QuantumDot/Spectra/XYReflexionSymmetricPeriodicSpectra.h"

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

XYReflexionSymmetricPeriodicSpectra::XYReflexionSymmetricPeriodicSpectra(XYReflexionSymmetricPeriodic3DOneParticle* space, char* fileName)
{
  this->NbrStateX = space->GetNbrStateX();
  this->NbrStateY = space->GetNbrStateY();
  this->NbrStateZ = space->GetNbrStateZ();
  this->LowerImpulsionX = space->GetLowerImpulsionX();
  this->LowerImpulsionY = space->GetLowerImpulsionY();
  this->LowerImpulsionZ = space->GetLowerImpulsionZ();

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }

  this->RealCoefficients = new double** [this->NbrStateX];
  this->ImaginaryCoefficients = new double** [this->NbrStateX];
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      this->RealCoefficients[i] = new double* [this->NbrStateY];
      this->ImaginaryCoefficients[i] = new double* [this->NbrStateY];   
      for (int j = 0; j < this->NbrStateY; ++j)
	{
	  this->RealCoefficients[i][j] = new double [this->NbrStateZ];
	  this->ImaginaryCoefficients[i][j] = new double [this->NbrStateZ];
	  for (int k = 0; k < this->NbrStateZ; ++k)
	    File >> this->RealCoefficients[i][j][k] >> this->ImaginaryCoefficients[i][j][k];
	}
    }
  File.close();
}

// get the value of impulsion operators with another wavefunction <this|p|another>
//
// space = Hilbert space describing the other particle
// fileName = the file to stock the other function
// sizeX, sizeY, sizeZ = size of sample in X, Y and Z directions
// impulsionX, impulsionY, impulsionZ = reference to the return values

void XYReflexionSymmetricPeriodicSpectra::GetImpulsion(XYReflexionSymmetricPeriodic3DOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &realImpulsionX, double &imaginaryImpulsionX, double &realImpulsionY, double &imaginaryImpulsionY, double &realImpulsionZ, double &imaginaryImpulsionZ)
{
  int nbrStateX = space->GetNbrStateX();
  int nbrStateY = space->GetNbrStateY();
  int nbrStateZ = space->GetNbrStateZ();
  int lowerImpulsionX = space->GetLowerImpulsionX();
  int lowerImpulsionY = space->GetLowerImpulsionY();
  int lowerImpulsionZ = space->GetLowerImpulsionZ();

  int MaxLowerX = this->LowerImpulsionX, MaxLowerY = this->LowerImpulsionY, MaxLowerZ = this->LowerImpulsionZ;
  if (this->LowerImpulsionX < lowerImpulsionX)
    MaxLowerX = lowerImpulsionX;
  if (this->LowerImpulsionY < lowerImpulsionY)
    MaxLowerY = lowerImpulsionY;
  if (this->LowerImpulsionZ < lowerImpulsionZ)
    MaxLowerZ = lowerImpulsionZ;

  int MinUpperX = this->LowerImpulsionX + this->NbrStateX, MinUpperY = this->LowerImpulsionY + this->NbrStateY, MinUpperZ = this->LowerImpulsionZ + this->NbrStateZ;
  if (MinUpperX > (lowerImpulsionX + nbrStateX))
    MinUpperX = lowerImpulsionX + nbrStateX;
  if (MinUpperY > (lowerImpulsionY + nbrStateY))
    MinUpperY = lowerImpulsionY + nbrStateY;
  if (MinUpperZ > (lowerImpulsionZ + nbrStateZ))
    MinUpperZ = lowerImpulsionZ + nbrStateZ;

  cout << MaxLowerX << '\t' <<  MaxLowerY << '\t' <<  MaxLowerZ << endl;
  cout << MinUpperX << '\t' << MinUpperY << '\t' << MinUpperZ << endl;

  int MinStateX = MinUpperX - MaxLowerX, MinStateY = MinUpperY - MaxLowerY, MinStateZ = MinUpperZ - MaxLowerZ;

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Error in open the file: " << fileName << endl;
      exit(0);
    }
  double*** realCoefficients = new double** [nbrStateX];
  double*** imaginaryCoefficients = new double** [nbrStateX];
  for (int i = 0; i < nbrStateX; ++i)
    {
      realCoefficients[i] = new double* [nbrStateY];
      imaginaryCoefficients[i] = new double* [nbrStateY];   
      for (int j = 0; j < nbrStateY; ++j)
	{
	  realCoefficients[i][j] = new double [nbrStateZ];
	  imaginaryCoefficients[i][j] = new double [nbrStateZ];
	  for (int k = 0; k < nbrStateZ; ++k)
	    File >> realCoefficients[i][j][k] >> imaginaryCoefficients[i][j][k];
	}
    }
  File.close();  


  realImpulsionX = 0.0;
  realImpulsionY = 0.0;
  realImpulsionZ = 0.0;
  imaginaryImpulsionX = 0.0;
  imaginaryImpulsionY = 0.0;
  imaginaryImpulsionZ = 0.0;
  double TmpRe = 0.0, TmpIm = 0.0;
  double* Re1; double* Im1; double* Re2; double* Im2;
  // x direction
  if ((this->LowerImpulsionY != lowerImpulsionY) || (this->LowerImpulsionX == lowerImpulsionX))
    {
      realImpulsionX = 0.0;
      imaginaryImpulsionX = 0.0;
    }
  else
    {
      for (int m = MaxLowerX; m < MinUpperX; ++m)
	{
	  TmpRe = 0.0; TmpIm = 0.0;
	  for (int n = MaxLowerY; n < MinUpperY; ++n)
	    {
	      Re1 = this->RealCoefficients[m - this->LowerImpulsionX][n - this->LowerImpulsionY];
	      Im1 = this->ImaginaryCoefficients[m - this->LowerImpulsionX][n - this->LowerImpulsionY];
	      Re2 = realCoefficients[m - lowerImpulsionX][n - lowerImpulsionY];
	      Im2 = imaginaryCoefficients[m - lowerImpulsionX][n - lowerImpulsionY];
	      for (int p = MaxLowerZ; p < MinUpperZ; ++p)
		{
		  TmpRe += (Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
		  TmpRe += (Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] - Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);
		}
	    }
	  realImpulsionX += (double(m) * TmpRe); imaginaryImpulsionX += (double(m) * TmpIm);
	}
      realImpulsionX = realImpulsionX / sizeX; imaginaryImpulsionX = imaginaryImpulsionX / sizeX;
    }

  // y direction
  if ((this->LowerImpulsionX != lowerImpulsionX) || (this->LowerImpulsionY == lowerImpulsionY))
    {
      realImpulsionY = 0.0;
      imaginaryImpulsionY = 0.0;
    }
  else
    {
      for (int m = MaxLowerX; m < MinUpperX; ++m)
	{
	  for (int n = MaxLowerY; n < MinUpperY; ++n)
	    {
	      TmpRe = 0.0; TmpIm = 0.0;
	      Re1 = this->RealCoefficients[m - this->LowerImpulsionX][n - this->LowerImpulsionY];
	      Im1 = this->ImaginaryCoefficients[m - this->LowerImpulsionX][n - this->LowerImpulsionY];
	      Re2 = realCoefficients[m - lowerImpulsionX][n - lowerImpulsionY];
	      Im2 = imaginaryCoefficients[m - lowerImpulsionX][n - lowerImpulsionY];
	      for (int p = MaxLowerZ; p < MinUpperZ; ++p)
		{
		  TmpRe += (Re1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ] + Im1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ]);
		  TmpRe += (Re1[p - this->LowerImpulsionZ] * Im2[p - lowerImpulsionZ] - Im1[p - this->LowerImpulsionZ] * Re2[p - lowerImpulsionZ]);
		}
	      realImpulsionY += (double(n) * TmpRe); imaginaryImpulsionY += (double(n) * TmpIm);
	    }
	}
      realImpulsionY = realImpulsionY / sizeY; imaginaryImpulsionY = realImpulsionY / sizeY;
    }

  double tr = 0.0, ti = 0.0;
  // z direction
  if ((this->LowerImpulsionX != lowerImpulsionX) || (this->LowerImpulsionY != lowerImpulsionY))
    {
      realImpulsionZ = 0.0;
      imaginaryImpulsionZ = 0.0;
    }
  else
    {
      /*
      for (int m = MaxLowerX; m < MinUpperX; ++m)
	{
	  for (int n = MaxLowerY; n < MinUpperY; ++n)
	    {
      */
      for (int m = 0; m < this->NbrStateX; ++m)
	{
	  for (int n = 0; n < this->NbrStateY; ++n)
	    {
	      Re1 = this->RealCoefficients[m][n];
	      Im1 = this->ImaginaryCoefficients[m][n];
	      Re2 = realCoefficients[m][n];
	      Im2 = imaginaryCoefficients[m][n];
	      for (int p = 0; p < this->NbrStateZ; ++p)
		{
		  tr += (Re1[p] * Re2[p] + Im1[p] * Im2[p]);
		  ti += (Re1[p] * Im2[p] - Im1[p] * Re2[p]);
		  realImpulsionZ += ((Re1[p] * Re2[p] + Im1[p] * Im2[p]) * double(p));
		  imaginaryImpulsionZ += ((Re1[p] * Im2[p] - Im1[p] * Re2[p]) * double(p));
		}
	    }
	}
      realImpulsionZ = realImpulsionZ / sizeZ; imaginaryImpulsionZ = realImpulsionZ / sizeZ;	    
      cout << tr << '\t' << ti << endl;
    }  
}
