////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 2D tight binding model                  //
//                                                                            //
//                        last modification : 01/05/2012                      //
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
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "GeneralTools/Endian.h"

#include <fstream>
#include <iostream>
#include <sys/time.h>


using std::ofstream;
using std::endl;
using std::cout;
using std::max;
using std::min;


// default constructor
//

Abstract2DTightBindingModel::Abstract2DTightBindingModel()
{
}

// destructor
//

Abstract2DTightBindingModel::~Abstract2DTightBindingModel()
{
}

// write an header that describes the tight binding model
// 
// output = reference on the output stream
// return value  = reference on the output stream

ofstream& Abstract2DTightBindingModel::WriteHeader(ofstream& output)
{
  int Dimension = 2;
  int HeaderSize = ((2 * Dimension) * sizeof(double)) + ((Dimension + 1) * sizeof(int));
  WriteLittleEndian(output, HeaderSize);
  WriteLittleEndian(output, Dimension);
  WriteLittleEndian(output, this->NbrSiteX);
  WriteLittleEndian(output, this->KxFactor);
  WriteLittleEndian(output, this->GammaX);
  WriteLittleEndian(output, this->NbrSiteY);
  WriteLittleEndian(output, this->KyFactor);
  WriteLittleEndian(output, this->GammaY);
  return output; 
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool Abstract2DTightBindingModel::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(kx, ky);
	  File << kx << " " << ky; 
	  for (int i = 0; i < this->NbrBands; ++i)
	    File << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex];
	  File << endl;
	}
    }
  File.close();
  return true;
}

// write the full band structure information in an ASCII file
//
// fileName = name of the output file 
// return value = true if no error occured  

bool Abstract2DTightBindingModel::WriteBandStructureASCII(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  for (int i = 0; i < this->NbrBands; ++i)
    for (int j = 0; j < this->NbrBands; ++j)
      File <<  "    U_{" << i << ", " << j << "}";
  File << endl;
  Complex Tmp;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(kx, ky);
	  File << kx << " " << ky; 
	  for (int i = 0; i < this->NbrBands; ++i)
	    File << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex];
	  for (int i = 0; i < this->NbrBands; ++i)
	    for (int j = 0; j < this->NbrBands; ++j)
	      {
		this->GetOneBodyMatrix(LinearizedMomentumIndex).GetMatrixElement(i, j, Tmp);
		File <<  "    " << Tmp;
	      }
	  File << endl;
	}
    }

  File.close();
  return true;
}


// write the eigenvalues of the D matrix in an ASCII file
//
// fileName = name of the ASCII file 
//nbrOccupiedBands = nbr of occupied bands
// return value = true if no error occured

bool Abstract2DTightBindingModel::WriteAsciiDMatrixEigenValues(char* fileName, int nbrOccupiedBands)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# ky" ;
  for (int i = 0; i < nbrOccupiedBands; ++i)
    File <<  "   DEigenValue_" << i << "    Theta_" << i;
  File << endl;
  
  double distancePlus;
  double distanceMoins;
  double distanceMod2PiPlus;
  double distanceMod2PiMoins;
 
  double theta1;
  double theta2;
  
  Complex** Lambda = this->ComputeDMatrixEigenvalues(nbrOccupiedBands, 0, this->NbrSiteY - 1, this->NbrSiteY); 
  double** Theta = new double*[this->NbrSiteY];
  for (int ky = 0; ky < this->NbrSiteY; ++ky)
  {
   Theta[ky] = new double[2];
   for (int i = 0; i < 2; ++i)
   {
     theta1 = atan2(Lambda[ky][nbrOccupiedBands - 2].Im,Lambda[ky][nbrOccupiedBands - 2].Re);
     theta2 = atan2(Lambda[ky][nbrOccupiedBands - 1].Im,Lambda[ky][nbrOccupiedBands - 1].Re);
     Theta[ky][0] = max(theta1, theta2);
     Theta[ky][1] = min(theta1, theta2); 
   }
  }
  
  for (int ky = 0; ky < this->NbrSiteY - 1; ++ ky)
  {
    distancePlus = abs(Theta[ky][0] - Theta[ky + 1][0]);
    distanceMod2PiPlus = abs(Theta[ky][0] - Theta[ky + 1][1] - 2*M_PI);
    distanceMoins = abs(Theta[ky][1] - Theta[ky + 1][1]);
    distanceMod2PiMoins = abs(Theta[ky][1] - Theta[ky + 1][0] + 2*M_PI);
    
    if (distanceMod2PiPlus < distancePlus)
    {
     double Tmp = Theta[ky + 1][0];
     Theta[ky + 1][0] = Theta[ky + 1][1] + 2*M_PI;
     Theta[ky + 1][1] = Tmp;
    }
    
    if (distanceMod2PiMoins < distanceMoins)
    {
     double Tmp = Theta[ky + 1][1];
     Theta[ky + 1][1] = Theta[ky + 1][0] - 2*M_PI;
     Theta[ky + 1][0] = Tmp;
    }
  }
  for (int ky = 0; ky < this->NbrSiteY; ++ky)
    {
      File << ky; 
      File << " " << Lambda[ky][0] << " " << atan2(Lambda[ky][0].Im,Lambda[ky][0].Re) << " " << Lambda[ky][1] << " "<< atan2(Lambda[ky][1].Im,Lambda[ky][1].Re) << " " << Theta[ky][0] << " " << Theta[ky][1] ;
      File << endl;
    }
   
  File.close();
  return true;
}

// compute the Chern number of a given band
//
// band = band index
// return value = Chern number

double Abstract2DTightBindingModel::ComputeChernNumber(int band)
{
  if (this->HaveOneBodyBasis() == false)
    {
      cout << "error, the tight binding model does not provide the one body basis" << endl;
      return 0.0;
    }
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  Complex TmpChernNumber = 0.0;
  Complex Tmp1[4];
  Complex Tmp2[8];
  for (long LinearizedMomentumIndex = 0l; LinearizedMomentumIndex < this->NbrStatePerBand; ++LinearizedMomentumIndex)
    {
      int Kx;
      int Ky;
      this->GetLinearizedMomentumIndex(LinearizedMomentumIndex, Kx, Ky);
      int LinearizedMomentumIndexIncX = this->GetLinearizedMomentumIndex((Kx + 1) % this->NbrSiteX, Ky);
      int LinearizedMomentumIndexDecX;
      if (Kx > 0)
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex((Kx - 1) % this->NbrSiteX, Ky);
      else
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex(this->NbrSiteX - 1, Ky);
      int LinearizedMomentumIndexIncY = this->GetLinearizedMomentumIndex(Kx, (Ky + 1) % this->NbrSiteY);
      int LinearizedMomentumIndexDecY;
      if (Ky > 0)
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, (Ky - 1) % this->NbrSiteY);
      else
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, this->NbrSiteY - 1);

      ComplexMatrix& LocalBasis = this->OneBodyBasis[LinearizedMomentumIndex];
      ComplexMatrix& LocalBasisIncX = this->OneBodyBasis[LinearizedMomentumIndexIncX];
      ComplexMatrix& LocalBasisDecX = this->OneBodyBasis[LinearizedMomentumIndexDecX];
      ComplexMatrix& LocalBasisIncY = this->OneBodyBasis[LinearizedMomentumIndexIncY];
      ComplexMatrix& LocalBasisDecY = this->OneBodyBasis[LinearizedMomentumIndexDecY];  
      Tmp1[0] = 0.0;
      Tmp1[1] = 0.0;
      Tmp1[2] = 0.0;
      Tmp1[3] = 0.0;

      Tmp2[0] = 0.0;
      Tmp2[1] = 0.0;
      Tmp2[2] = 0.0;
      Tmp2[3] = 0.0;
      Tmp2[4] = 0.0;
      Tmp2[5] = 0.0;
      Tmp2[6] = 0.0;
      Tmp2[7] = 0.0;

      for (int i = 0; i < this->NbrBands; ++i)
	{
	  Tmp1[0] += LocalBasis[band][i] * Conj(LocalBasisIncX[band][i]);
	  Tmp1[1] += LocalBasis[band][i] * Conj(LocalBasisDecX[band][i]);
	  Tmp1[2] += LocalBasis[band][i] * Conj(LocalBasisIncY[band][i]);
	  Tmp1[3] += LocalBasis[band][i] * Conj(LocalBasisDecY[band][i]);

	  Tmp2[0] += Conj(LocalBasisIncX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[1] += Conj(LocalBasisDecX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[2] += Conj(LocalBasisIncX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[3] += Conj(LocalBasisDecX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[4] += Conj(LocalBasisIncY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[5] += Conj(LocalBasisDecY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[6] += Conj(LocalBasisIncY[band][i]) * LocalBasisDecX[band][i];
	  Tmp2[7] += Conj(LocalBasisDecY[band][i]) * LocalBasisDecX[band][i];
	}

      TmpChernNumber += (Tmp1[2] * Conj(Tmp1[0]) * Tmp2[0]);
      TmpChernNumber -= (Tmp1[2] * Conj(Tmp1[1]) * Tmp2[1]);
      TmpChernNumber -= (Tmp1[3] * Conj(Tmp1[0]) * Tmp2[2]);
      TmpChernNumber += (Tmp1[3] * Conj(Tmp1[1]) * Tmp2[3]);
	  
      TmpChernNumber -= (Tmp1[0] * Conj(Tmp1[2]) * Tmp2[4]);
      TmpChernNumber += (Tmp1[0] * Conj(Tmp1[3]) * Tmp2[5]);
      TmpChernNumber += (Tmp1[1] * Conj(Tmp1[2]) * Tmp2[6]);
      TmpChernNumber -= (Tmp1[1] * Conj(Tmp1[3]) * Tmp2[7]);

    }
  TmpChernNumber /= 8.0 * M_PI;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
  cout << "Chern number computed in  " << Dt << "s" << endl;
  return TmpChernNumber.Im;
}

// compute the Berry curvature  of a given band
//
// band = band index
// fileName = name of the output file 
// return value = Chern number

double Abstract2DTightBindingModel::ComputeBerryCurvature(int band, char* fileName)
{
  if (this->HaveOneBodyBasis() == false)
    {
      cout << "error, the tight binding model does not provide the one body basis" << endl;
      return 0.0;
    }
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  Complex TmpChernNumber = 0.0;
  Complex Tmp1[4];
  Complex Tmp2[8];
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky    Berry_curvature";
  double Fluctations = 0.0;
  for (long LinearizedMomentumIndex = 0l; LinearizedMomentumIndex < this->NbrStatePerBand; ++LinearizedMomentumIndex)
    {
      int Kx;
      int Ky;
      this->GetLinearizedMomentumIndex(LinearizedMomentumIndex, Kx, Ky);
      int LinearizedMomentumIndexIncX = this->GetLinearizedMomentumIndex((Kx + 1) % this->NbrSiteX, Ky);
      int LinearizedMomentumIndexDecX;
      if (Kx > 0)
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex((Kx - 1) % this->NbrSiteX, Ky);
      else
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex(this->NbrSiteX - 1, Ky);
      int LinearizedMomentumIndexIncY = this->GetLinearizedMomentumIndex(Kx, (Ky + 1) % this->NbrSiteY);
      int LinearizedMomentumIndexDecY;
      if (Ky > 0)
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, (Ky - 1) % this->NbrSiteY);
      else
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, this->NbrSiteY - 1);

      ComplexMatrix& LocalBasis = this->OneBodyBasis[LinearizedMomentumIndex];
      ComplexMatrix& LocalBasisIncX = this->OneBodyBasis[LinearizedMomentumIndexIncX];
      ComplexMatrix& LocalBasisDecX = this->OneBodyBasis[LinearizedMomentumIndexDecX];
      ComplexMatrix& LocalBasisIncY = this->OneBodyBasis[LinearizedMomentumIndexIncY];
      ComplexMatrix& LocalBasisDecY = this->OneBodyBasis[LinearizedMomentumIndexDecY];  
      Tmp1[0] = 0.0;
      Tmp1[1] = 0.0;
      Tmp1[2] = 0.0;
      Tmp1[3] = 0.0;

      Tmp2[0] = 0.0;
      Tmp2[1] = 0.0;
      Tmp2[2] = 0.0;
      Tmp2[3] = 0.0;
      Tmp2[4] = 0.0;
      Tmp2[5] = 0.0;
      Tmp2[6] = 0.0;
      Tmp2[7] = 0.0;

      for (int i = 0; i < this->NbrBands; ++i)
	{
	  Tmp1[0] += LocalBasis[band][i] * Conj(LocalBasisIncX[band][i]);
	  Tmp1[1] += LocalBasis[band][i] * Conj(LocalBasisDecX[band][i]);
	  Tmp1[2] += LocalBasis[band][i] * Conj(LocalBasisIncY[band][i]);
	  Tmp1[3] += LocalBasis[band][i] * Conj(LocalBasisDecY[band][i]);

	  Tmp2[0] += Conj(LocalBasisIncX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[1] += Conj(LocalBasisDecX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[2] += Conj(LocalBasisIncX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[3] += Conj(LocalBasisDecX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[4] += Conj(LocalBasisIncY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[5] += Conj(LocalBasisDecY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[6] += Conj(LocalBasisIncY[band][i]) * LocalBasisDecX[band][i];
	  Tmp2[7] += Conj(LocalBasisDecY[band][i]) * LocalBasisDecX[band][i];
	}

      Complex TmpCurvature = 0.0;
      TmpCurvature += (Tmp1[2] * Conj(Tmp1[0]) * Tmp2[0]);
      TmpCurvature -= (Tmp1[2] * Conj(Tmp1[1]) * Tmp2[1]);
      TmpCurvature -= (Tmp1[3] * Conj(Tmp1[0]) * Tmp2[2]);
      TmpCurvature += (Tmp1[3] * Conj(Tmp1[1]) * Tmp2[3]);
	  
      TmpCurvature -= (Tmp1[0] * Conj(Tmp1[2]) * Tmp2[4]);
      TmpCurvature += (Tmp1[0] * Conj(Tmp1[3]) * Tmp2[5]);
      TmpCurvature += (Tmp1[1] * Conj(Tmp1[2]) * Tmp2[6]);
      TmpCurvature -= (Tmp1[1] * Conj(Tmp1[3]) * Tmp2[7]);

      TmpCurvature *= 0.25;

      Fluctations += (TmpCurvature.Im - (2.0 * M_PI / ((double) this->NbrStatePerBand))) * (TmpCurvature.Im - (2.0 * M_PI / ((double) this->NbrStatePerBand)));

      File << Kx << " " << Ky << " " << TmpCurvature.Im << endl;

      TmpChernNumber += TmpCurvature;
    }
//  Fluctations *= ((double) this->NbrStatePerBand) ;
  cout << "Berry curvature fluctuations = " << Fluctations << " " <<  sqrt(Fluctations) << " " << (TmpChernNumber.Im * TmpChernNumber.Im) << " " << sqrt(Fluctations - (TmpChernNumber.Im * TmpChernNumber.Im)) 
       << "( " << (sqrt(Fluctations - (TmpChernNumber.Im * TmpChernNumber.Im)) / (2.0 * M_PI) )<< ") in 2 pi units" << endl;
  TmpChernNumber /= 2.0 * M_PI;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
  cout << "Chern number computed in  " << Dt << "s" << endl;
//  cout << "Berry curvature fluctuations = " << sqrt ()<< endl;

  File.close();

  return TmpChernNumber.Im;
}


//compute the complex eigenvalues of the D(ky) matrix (in order to compute the Z2 invariant)
//
//bandIndex = band index (corresponds to two bands that are related by time reversal symmetry)
//nbrOccupiedBands = dimension of the D matrix
//DMatrixEigenvalues = array of eigenvalues of the D Matrix, for all values of ky
//kyMin = minimal value of ky for which the D matrix has to be diagonalized
//kyMax = maximal value of ky for which the D matrix has to be diagonalized
//nbrKy = number of ky values for which the D matrix has to be diagonalized
//return value = array of eigenvalues of the D matrix
Complex** Abstract2DTightBindingModel::ComputeDMatrixEigenvalues(int nbrOccupiedBands, int kyMin, int kyMax, int nbrKy)
{
  Complex** DMatrixEigenvalues;
  DMatrixEigenvalues = new Complex*[this->NbrSiteY];
  for (int i = 0; i < this->NbrSiteY; ++i)
    DMatrixEigenvalues[i] = new Complex[nbrOccupiedBands];
  ComplexMatrix TmpDMatrix(nbrOccupiedBands, nbrOccupiedBands, true);
  ComplexMatrix FMatrix(nbrOccupiedBands, nbrOccupiedBands, true);
  
//   ComplexMatrix Rotation(this->NbrBands, this->NbrBands, true);
//   Rotation.SetMatrixElement(0, 0, M_SQRT1_2);
//   Rotation.SetMatrixElement(0, 1, Complex(0.0, M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 0, Complex(0.0, -1.0*M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 1, -1.0*M_SQRT1_2);
//   Rotation.SetMatrixElement(2, 2, 1.0);
//   Rotation.SetMatrixElement(3, 3, 1.0);
  
//   ComplexMatrix Rotation1(this->NbrBands, this->NbrBands, true);
//   Rotation.SetMatrixElement(0, 0, M_SQRT1_2);
//   Rotation.SetMatrixElement(0, 1, Complex(0.0, M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 0, Complex(0.0, -1.0*M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 1, -1.0*M_SQRT1_2);
//   Rotation1.SetMatrixElement(2, 2, 1.0);
//   Rotation1.SetMatrixElement(3, 3, 1.0);
  
  
  for (int ky = kyMin; ky <= kyMax; ++ky)
  {
    TmpDMatrix.SetToIdentity();
    for (int i = 0; i < this->NbrSiteX; ++i)
    {
      double KX = (double) i *2.0 * M_PI / ((double) this->NbrSiteX);
      double KX1 = (double) (i + 1) *2.0 * M_PI / ((double) this->NbrSiteX);
//       cout << cos(KX) << " " << sin(KX) << endl;
//       Rotation.SetMatrixElement(0, 0, Complex(M_SQRT1_2*cos(KX), M_SQRT1_2*sin(KX)));
//       Rotation.SetMatrixElement(0, 1, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 0, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 1, Complex(-1.0*M_SQRT1_2*cos(KX), M_SQRT1_2*sin(KX)));
// 
//       
//       Rotation1.SetMatrixElement(0, 0, Complex(M_SQRT1_2*cos(KX1), M_SQRT1_2*sin(KX1)));
//       Rotation1.SetMatrixElement(0, 1, M_SQRT1_2);
//       Rotation1.SetMatrixElement(1, 0, M_SQRT1_2);
//       Rotation1.SetMatrixElement(1, 1, Complex(-1.0*M_SQRT1_2*cos(KX1), M_SQRT1_2*sin(KX1)));
//       Rotation.SetMatrixElement(0, 0, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 1, -1.0*M_SQRT1_2);
//       Rotation.SetMatrixElement(0, 1, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 0, M_SQRT1_2);
 
	int LinearizedMomentumIndex1 = this->GetLinearizedMomentumIndex(i, ky);
	int LinearizedMomentumIndex2 = this->GetLinearizedMomentumIndex((i + 1) % this->NbrSiteX, ky);
// 	cout << i << " " << LinearizedMomentumIndex1 << " " << LinearizedMomentumIndex2 << " " << endl;
	ComplexMatrix& LocalBasis = this->OneBodyBasis[LinearizedMomentumIndex1];
	ComplexMatrix& LocalBasisIncX = this->OneBodyBasis[LinearizedMomentumIndex2];
// 	ComplexMatrix LocalBasis(this->NbrBands, this->NbrBands, true);
// 	ComplexMatrix LocalBasisIncX(this->NbrBands, this->NbrBands, true);
	
// 	LocalBasis = TmpLocalBasis*Rotation;
// 	LocalBasisIncX = TmpLocalBasisIncX*Rotation1;
	
// 	LocalBasis = TmpLocalBasis;
// 	LocalBasisIncX = TmpLocalBasisIncX;
	
	for (int n = 0; n < nbrOccupiedBands; ++n)
	  {
	    for (int m = 0; m < nbrOccupiedBands; ++m)
	      {
// 		Complex Tmp = 0.0;
// 		for (int alpha = 0; alpha < this->NbrBands; ++alpha)
// 		{
// 		  Tmp += Conj(LocalBasis[n][alpha]) * LocalBasisIncX[m][alpha];
// 		}
// 		FMatrix.SetMatrixElement(n, m, Tmp);
		FMatrix.SetMatrixElement(n, m, LocalBasis[n] * LocalBasisIncX[m]);
	      }
	  }
// 	  cout << i << endl;
// 	  cout << FMatrix << endl;
// 	  ComplexMatrix TmpMatrix = TmpDMatrix;
// 	  TmpDMatrix = TmpMatrix*FMatrix;
	  TmpDMatrix.Multiply(FMatrix);
      }
    
// 	  if (ky == 0)
// 	  {
// 	   cout << "ky = 0" << endl;
// 	   cout << TmpDMatrix << endl; 
// 	  }
    ComplexDiagonalMatrix TmpDiag(nbrOccupiedBands);
#ifdef __LAPACK__
    TmpDMatrix.LapackDiagonalize(TmpDiag);
#else
    TmpDMatrix.Diagonalize(TmpDiag);
#endif
//     cout << TmpDiag << endl;
    
    for (int j = 0; j < nbrOccupiedBands; ++j)
    {
      DMatrixEigenvalues[ky][j] = TmpDiag[j] ;
//       cout << DMatrixEigenvalues[ky][j] << endl;
    }
  }
  
  return DMatrixEigenvalues;
}

//compute the Z2 topological invariant for a system with time reversal symmetry
//
//nbrOccupiedBands = number of occupied bands
//return value = Z2 invariant
int Abstract2DTightBindingModel::ComputeZ2Invariant(int nbrOccupiedBands)
{
  int z2Invariant = 0;
  double referenceLine = 0.9267;
  
  double distancePlus;
  double distanceMoins;
  double distanceMod2PiPlus;
  double distanceMod2PiMoins;
  
  int ModPiPlus = 0;
  int ModPiMoins = 0;
  
  double theta1;
  double theta2;
  
  Complex** Lambda = this->ComputeDMatrixEigenvalues(nbrOccupiedBands, 0, this->NbrSiteY - 1, this->NbrSiteY); 
  double** Theta = new double*[this->NbrSiteY];
  for (int ky = 0; ky < this->NbrSiteY; ++ky)
  {
   Theta[ky] = new double[2];
   for (int i = 0; i < 2; ++i)
   {
     theta1 = atan2(Lambda[ky][nbrOccupiedBands - 2].Im,Lambda[ky][nbrOccupiedBands - 2].Re);
     theta2 = atan2(Lambda[ky][nbrOccupiedBands - 1].Im,Lambda[ky][nbrOccupiedBands - 1].Re);
     Theta[ky][0] = max(theta1, theta2);
     Theta[ky][1] = min(theta1, theta2); 
   }
  }
  
  for (int ky = 0; ky < this->NbrSiteY  - 1; ++ ky)
  {
    distancePlus = abs(Theta[ky][0] - Theta[ky + 1][0]);
    distanceMod2PiPlus = abs(Theta[ky][0] - Theta[ky + 1][1] - 2*M_PI);
    distanceMoins = abs(Theta[ky][1] - Theta[ky + 1][1]);
    distanceMod2PiMoins = abs(Theta[ky][1] - Theta[ky + 1][0] + 2*M_PI);
    
    if (distanceMod2PiPlus < distancePlus)
    {
     ModPiPlus += 1;
     double Tmp = Theta[ky + 1][0];
     Theta[ky + 1][0] = Theta[ky + 1][1] + 2*M_PI;
     Theta[ky + 1][1] = Tmp;
    }
    
    if (distanceMod2PiMoins < distanceMoins)
    {
     ModPiMoins += 1;
     double Tmp = Theta[ky + 1][1];
     Theta[ky + 1][1] = Theta[ky + 1][0] - 2*M_PI;
     Theta[ky + 1][0] = Tmp;
    }
  }
//   cout << ModPiPlus << " " << ModPiMoins << endl;
  for (int ky = 1 ; ky < this->NbrSiteY/2 - 1; ++ky)
  {
//     cout << ky << " " << Theta[ky][0] << " " << Theta[ky][1] << endl;
    for (int i = 0; i <= ModPiPlus; ++i)
    {
//       cout << Theta[ky + 1][0] - (referenceLine + i*2*M_PI) << " " << Theta[ky][0] - (referenceLine + i*2*M_PI) << " " << (Theta[ky + 1][1] - (referenceLine - i*2*M_PI)) << " " << (Theta[ky][1] - (referenceLine - i*2*M_PI)) << endl;
      if ((Theta[ky + 1][0] - (referenceLine + i*2*M_PI)) * (Theta[ky][0] - (referenceLine + i*2*M_PI)) < 0)
	z2Invariant += 1;
    }
    for (int i = 0; i <= ModPiMoins; ++i)
    {
      if ((Theta[ky + 1][1] - (referenceLine - i*2*M_PI)) * (Theta[ky][1] - (referenceLine - i*2*M_PI)) < 0)
	z2Invariant += 1;
    }
  }
  return (z2Invariant % 2); 
}