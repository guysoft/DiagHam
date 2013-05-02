////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 1D tight binding model                  //
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
#include "Tools/FTITightBinding/Abstract1DTightBindingModel.h"
#include "GeneralTools/Endian.h"

#include <fstream>


using std::ofstream;
using std::endl;


// default constructor
//

Abstract1DTightBindingModel::Abstract1DTightBindingModel()
{
  this->EnergyBandStructure = 0;
  this->OneBodyBasis = 0;
  this->EmbeddingX = RealVector();
}

// destructor
//

Abstract1DTightBindingModel::~Abstract1DTightBindingModel()
{
  if (this->OneBodyBasis != 0)
    {
      delete[] this->OneBodyBasis;
    }
  if (this->EnergyBandStructure != 0)
    {
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  delete[] this->EnergyBandStructure[i];
	}
      delete[] this->EnergyBandStructure;
    }
}

// write an header that describes the tight binding model
// 
// output = reference on the output stream
// return value  = reference on the output stream

ofstream& Abstract1DTightBindingModel::WriteHeader(ofstream& output)
{
  int Dimension = 1;
  int HeaderSize = (((this->NbrBands + 2) * Dimension) * sizeof(double)) + ((Dimension + 1) * sizeof(int));
  WriteLittleEndian(output, HeaderSize);
  WriteLittleEndian(output, Dimension);
  WriteLittleEndian(output, this->NbrSiteX);
  WriteLittleEndian(output, this->KxFactor);
  WriteLittleEndian(output, this->GammaX);
  if (this->EmbeddingX.GetVectorDimension() != this->NbrBands)
  {
      double Tmp = 0.0;
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, Tmp);
  }
  else
  {
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, this->EmbeddingX[i]);
  }
  return output; 
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool Abstract1DTightBindingModel::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      File << kx; 
      for (int i = 0; i < this->NbrBands; ++i)
	File << " " << this->EnergyBandStructure[i][kx];
      File << endl;
    }
  File.close();
  return true;
}
