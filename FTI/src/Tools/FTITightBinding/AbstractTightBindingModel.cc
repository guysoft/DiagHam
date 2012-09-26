////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of abstract tight binding model                    //
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
#include "Tools/FTITightBinding/AbstractTightBindingModel.h"
#include "GeneralTools/Endian.h"

#include <fstream>


using std::ofstream;
using std::ios;
using std::endl;


// destructor
//

AbstractTightBindingModel::~AbstractTightBindingModel()
{
}

// write an ASCII header that describes the tight binding model
// 
// output = reference on the output stream
// commentChar = optional ASCII character to add in front of each header line
// return value  = reference on the output stream

ostream& AbstractTightBindingModel::WriteASCIIHeader(ostream& output, char commentChar)
{
  return output;
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool AbstractTightBindingModel::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# index";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrStatePerBand; ++kx)
    {
      File << kx; 
      for (int i = 0; i < this->NbrBands; ++i)
	File << " " << this->GetEnergy(i, kx);
      File << endl;
    }
  File.close();
  return true;
}

// compute the spread of a given band
// 
// band = band index
// return value = band spread 

double AbstractTightBindingModel::ComputeBandSpread(int band)
{
  if ((band < 0) || (band > this->NbrBands))
    {
      return -1.0;
    }
  double MinBand = this->GetEnergy(band, 0);
  double MaxBand = this->GetEnergy(band, 0);
  for (int i = 1; i < this->NbrStatePerBand; ++i)
    {
      if (this->GetEnergy(band, i) > MaxBand)
	{
	  MaxBand = this->GetEnergy(band, i);
	}
      else
	{
	  if (this->GetEnergy(band, i) < MinBand)
	    {
	      MinBand = this->GetEnergy(band, i);
	    }
	}
    }
  return (MaxBand - MinBand);
}

// compute the direct gap between two bands
// 
// band1 = index of the lower band 
// band2 = index of the upper band (-1 if it has to be set to band1 + 1)
// return value =  direct band gap

double AbstractTightBindingModel::ComputeDirectBandGap(int band1, int band2)
{
  if ((band1 < 0) || (band1 > this->NbrBands))
    {
      return -1.0;
    }
  if (band2 < 0)
    band2 = band1 + 1;
  if ((band2 < 0) || (band2 > this->NbrBands) || (band2 <= band1))
    {
      return -1.0;
    }
  double Gap = this->GetEnergy(band2, 0) - this->GetEnergy(band1, 0);
  for (int i = 1; i < this->NbrStatePerBand; ++i)
    {
      double TmpGap = this->GetEnergy(band2, i) - this->GetEnergy(band1, i);
      if (TmpGap < Gap)
	Gap = TmpGap;
    }
  return Gap;
}

// write the full band structure information in a binary file
//
// fileName = name of the output file 
// return value = true if no error occured  

bool AbstractTightBindingModel::WriteBandStructure(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  WriteLittleEndian(File, this->NbrBands);
  WriteLittleEndian(File, this->NbrStatePerBand);
  for (int kx = 0; kx < this->NbrStatePerBand; ++kx)
    for (int i = 0; i < this->NbrBands; ++i)
      {
	double Tmp = this->GetEnergy(i, kx);
	WriteLittleEndian(File, Tmp);
      }
  if (this->HaveOneBodyBasis() == true)
    {
      for (int kx = 0; kx < this->NbrStatePerBand; ++kx)
	this->GetOneBodyMatrix(kx).WriteMatrix(File);
    }
  File.close();
  return true;
}

// write the full band structure information in an ASCII file
//
// fileName = name of the output file 
// return value = true if no error occured  

bool AbstractTightBindingModel::WriteBandStructureASCII(char* fileName)
{
  ofstream File;
  this->WriteASCIIHeader(File, '#');
  File << "# index";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  for (int i = 0; i < this->NbrBands; ++i)
    for (int j = 0; j < this->NbrBands; ++j)
      File <<  "    U_{" << i << ", " << j << "}";
  File << endl;
  Complex Tmp;
  for (int LinearizedMomentumIndex = 0; LinearizedMomentumIndex < this->NbrStatePerBand; ++LinearizedMomentumIndex)
    {
      File << LinearizedMomentumIndex; 
      for (int i = 0; i < this->NbrBands; ++i)
	File << " " <<  this->GetEnergy(i, LinearizedMomentumIndex);
      for (int i = 0; i < this->NbrBands; ++i)
	for (int j = 0; j < this->NbrBands; ++j)
	  {
	    this->GetOneBodyMatrix(LinearizedMomentumIndex).GetMatrixElement(i, j, Tmp);
	    File <<  "    " << Tmp;
	  }
      File << endl;
    }

  File.close();
  return true;
}

// compute the band structure
//

void AbstractTightBindingModel::ComputeBandStructure()
{
  this->CoreComputeBandStructure(0l, 0l);
}

