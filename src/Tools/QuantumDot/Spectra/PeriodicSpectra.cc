////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                         class for periodic  spectra                        //
//                                                                            //
//                        last modification : 12/20/2003                      //
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


#include "Tools/QuantumDot/Spectra/PeriodicSpectra.h"

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
// space: Hilbert space describing the particle
// fileName: name of the state file

PeriodicSpectra::PeriodicSpectra(Periodic3DOneParticle* space, char* fileName)
{
  this->NbrStateX = space->GetNbrStateX();
  this->NbrStateY = space->GetNbrStateY();
  this->NbrStateZ = space->GetNbrStateZ();
  this->LowerImpulsionX = space->GetLowerImpulsionX();
  this->LowerImpulsionY = space->GetLowerImpulsionY();
  this->LowerImpulsionZ = space->GetLowerImpulsionZ();
  
  int LengthX = (this->NbrStateX - 1) * 2 + 1; int LengthY = (this->NbrStateY - 1) * 2 + 1; int LengthZ = (this->NbrStateZ - 1) * 2 + 1;
  int OriginX = this->NbrStateX - 1; int OriginY = this->NbrStateY - 1; int OriginZ = this->NbrStateZ - 1;
  this->RealWaveFunctionOverlapX = new double [LengthX];
  this->ImaginaryWaveFunctionOverlapX = new double [LengthX];
  this->RealSquareOverlapX = new double [LengthX];
  this->ImaginarySquareOverlapX = new double [LengthX];
  double Tmp = 0.0;
  for (int i = 0; i < LengthX; ++i)
    {
      if (i != OriginX)
	{
	  Tmp = 1.0 / (2.0 * M_PI * double(i - OriginX));
	  this->RealWaveFunctionOverlapX[i] = 0.0;
	  this->ImaginaryWaveFunctionOverlapX[i] = -Tmp;
	  this->RealSquareOverlapX[i] = Tmp * Tmp * 2.0;
	  this->ImaginarySquareOverlapX[i] = -Tmp;
	}
      else
	{
	  this->RealWaveFunctionOverlapX[i] = 0.5;
	  this->ImaginaryWaveFunctionOverlapX[i] = 0.0;	  
	  this->RealSquareOverlapX[i] = 1.0 / 3.0;
	  this->ImaginarySquareOverlapX[i] = 0.0;
	}
    }  
  this->RealWaveFunctionOverlapY = new double [LengthY];
  this->ImaginaryWaveFunctionOverlapY = new double [LengthY];
  this->RealSquareOverlapY = new double [LengthY];
  this->ImaginarySquareOverlapY = new double [LengthY];
  for (int i = 0; i < LengthY; ++i)
    {
      if (i != OriginY)
	{
	  Tmp =  1.0 / (2.0 * M_PI * double(i - OriginY));
	  this->RealWaveFunctionOverlapY[i] = 0.0;
	  this->ImaginaryWaveFunctionOverlapY[i] = -Tmp;
	  this->RealSquareOverlapY[i] = Tmp * Tmp * 2.0;
	  this->ImaginarySquareOverlapY[i] = -Tmp;
	}
      else
	{
	  this->RealWaveFunctionOverlapY[i] = 0.5;
	  this->ImaginaryWaveFunctionOverlapY[i] = 0.0;
	  this->RealSquareOverlapY[i] = 1.0 / 3.0;
	  this->ImaginarySquareOverlapY[i] = 0.0;
	}
    } 
  this->RealWaveFunctionOverlapZ = new double [LengthZ];
  this->ImaginaryWaveFunctionOverlapZ = new double [LengthZ];
  this->RealSquareOverlapZ = new double [LengthZ];
  this->ImaginarySquareOverlapZ = new double [LengthZ];
  for (int i = 0; i < LengthZ; ++i)
    {
      if (i != OriginZ)
	{
	  Tmp = 1.0 / (2.0 * M_PI * double(i - OriginZ));
	  this->RealWaveFunctionOverlapZ[i] = 0.0;
	  this->ImaginaryWaveFunctionOverlapZ[i] = -Tmp;
	  this->RealSquareOverlapZ[i] = Tmp * Tmp * 2.0;
	  this->ImaginarySquareOverlapZ[i] = -Tmp;
	}
      else
	{
	  this->RealWaveFunctionOverlapZ[i] = 0.5;
	  this->ImaginaryWaveFunctionOverlapZ[i] = 0.0;
 	  this->RealSquareOverlapZ[i] = 1.0 / 3.0;
	  this->ImaginarySquareOverlapZ[i] = 0.0;
	}
    } 

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

// get mean value in X direction
//
// return: position in 1.0 scale

double PeriodicSpectra::GetMeanValueX(double& squareX)
{
  int OriginX = this->NbrStateX - 1;
  double TmpRe = 0.0, TmpIm = 0.0;
  double Tmp = 0.0;
  squareX = 0.0;
  for (int m1 = 0; m1 < this->NbrStateX; ++m1)
    for (int m2 = 0; m2 < this->NbrStateX; ++m2)
      {
	TmpRe = 0.0; TmpIm = 0.0;
	for (int n = 0; n < this->NbrStateY; ++n)
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpRe += (RealCoefficients[m1][n][p] * RealCoefficients[m2][n][p] + ImaginaryCoefficients[m1][n][p] * ImaginaryCoefficients[m2][n][p]);
	      TmpIm += (RealCoefficients[m1][n][p] * ImaginaryCoefficients[m2][n][p] - ImaginaryCoefficients[m1][n][p] * RealCoefficients[m2][n][p]);
	    }
	Tmp += (TmpRe * this->RealWaveFunctionOverlapX[-m1 + m2 + OriginX] - TmpIm * this->ImaginaryWaveFunctionOverlapX[-m1 + m2 + OriginX]);   
  	squareX += (TmpRe * this->RealSquareOverlapX[m1 - m2 + OriginX] - TmpIm * this->ImaginarySquareOverlapX[-m1 + m2 + OriginX]);   
      }
  return Tmp;
}

// get mean value in Y direction
//
// return: position in 1.0 scale

double PeriodicSpectra::GetMeanValueY(double& squareY)
{
  int OriginY = this->NbrStateY - 1;
  double TmpRe = 0.0, TmpIm = 0.0;
  double Tmp = 0.0;
  squareY = 0.0;
  for (int n1 = 0; n1 < this->NbrStateY; ++n1)
    for (int n2 = 0; n2 < this->NbrStateY; ++n2)
      {
	TmpRe = 0.0; TmpIm = 0.0;
	for (int m = 0; m < this->NbrStateX; ++m)
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpRe += (RealCoefficients[m][n1][p] * RealCoefficients[m][n2][p] + ImaginaryCoefficients[m][n1][p] * ImaginaryCoefficients[m][n2][p]);
	      TmpIm += (RealCoefficients[m][n1][p] * ImaginaryCoefficients[m][n2][p] - ImaginaryCoefficients[m][n1][p] * RealCoefficients[m][n2][p]);
	    }
	Tmp += (TmpRe * this->RealWaveFunctionOverlapY[-n1 + n2 + OriginY] - TmpIm * this->ImaginaryWaveFunctionOverlapY[-n1 + n2 + OriginY]); 
 	squareY += (TmpRe * this->RealSquareOverlapY[-n1 + n2 + OriginY] - TmpIm * this->ImaginarySquareOverlapY[-n1 + n2 + OriginY]);        
      }
  return Tmp;
}

// get mean value in Z direction
//
// return: position in 1.0 scale

double PeriodicSpectra::GetMeanValueZ(double& squareZ)
{
  int OriginZ = this->NbrStateZ - 1;
  double TmpRe = 0.0, TmpIm = 0.0;
  double Tmp = 0.0;
  squareZ = 0.0;
  for (int p1 = 0; p1 < this->NbrStateZ; ++p1)
    for (int p2 = 0; p2 < this->NbrStateZ; ++p2)
      {
	TmpRe = 0.0; TmpIm = 0.0;
	for (int m = 0; m < this->NbrStateX; ++m)
	  for (int n = 0; n < this->NbrStateY; ++n)
	    {
	      TmpRe += (RealCoefficients[m][n][p1] * RealCoefficients[m][n][p2] + ImaginaryCoefficients[m][n][p1] * ImaginaryCoefficients[m][n][p2]);
	      TmpIm += (RealCoefficients[m][n][p1] * ImaginaryCoefficients[m][n][p2] - ImaginaryCoefficients[m][n][p1] * RealCoefficients[m][n][p2]);
	    }
	Tmp += (TmpRe * this->RealWaveFunctionOverlapZ[-p1 + p2 + OriginZ] - TmpIm * this->ImaginaryWaveFunctionOverlapZ[-p1 + p2 + OriginZ]);       
	squareZ += (TmpRe * this->RealSquareOverlapZ[-p1 + p2 + OriginZ] - TmpIm * this->ImaginarySquareOverlapZ[-p1 + p2 + OriginZ]);  
      }
  return Tmp;
}

// get the wave function value of a state at a given point
//
// x, y, z : the position of the point
// SizeX, SizeY, SizeZ : the 3D-sizes of the sample
// Real, Imaginary : references to the real and imaginary components of the wave function

void PeriodicSpectra::WaveFunctionValue(double x, double SizeX, double y, double SizeY, double z, double SizeZ, double& Real, double& Imaginary)
{
  double TmpX = 0.0; double TmpY = 0.0; double TmpZ = 0.0;
  double TmpRe = 0.0; double TmpIm = 0.0;
  double* TmpRealCoefficients; double* TmpImaginaryCoefficients;
  double TmpRe2 = 0.0; double TmpIm2 = 0.0;
  for (int m = 0; m < this->NbrStateX; ++m)
    {
      TmpX = double(m - this->LowerImpulsionX) * x / SizeX;
      for (int n = 0; n < this->NbrStateY; ++n)
	{
	  TmpY = double(n - this->LowerImpulsionY) * y /SizeY + TmpX;
	  TmpRealCoefficients = this->RealCoefficients[m][n];
	  TmpImaginaryCoefficients = this->ImaginaryCoefficients[m][n];
	  for (int p = 0; p < this->NbrStateZ; ++p)
	    {
	      TmpZ = 2 * M_PI * (TmpY + double(p - this->LowerImpulsionZ) * z / SizeZ);
	      TmpRe2 = cos(TmpZ); TmpIm2 = sin(TmpZ);
	      TmpRe += (TmpRealCoefficients[p] * TmpRe2 - TmpImaginaryCoefficients[p] * TmpIm2);
	      TmpRe += (TmpRealCoefficients[p] * TmpIm2 + TmpImaginaryCoefficients[p] * TmpRe2);
	    }
	}
    }
  Real = TmpRe; Imaginary = TmpIm;
}
