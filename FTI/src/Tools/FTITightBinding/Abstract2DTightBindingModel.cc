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

#include <fstream>
#include <iostream>
#include <sys/time.h>


using std::ofstream;
using std::endl;
using std::cout;


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

