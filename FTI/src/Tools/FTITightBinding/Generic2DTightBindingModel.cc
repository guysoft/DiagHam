////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of generic 2D tight binding model                 //
//                                                                            //
//                        last modification : 03/10/2012                      //
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
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;


// default constructor
//
// fileName = name of the binary file that contains the band structure information

Generic2DTightBindingModel::Generic2DTightBindingModel(char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return;
    }
  File.seekg (0l, ios::end);
  unsigned long FileSize = File.tellg ();
  File.close();

  File.open(fileName, ios::binary | ios::in);
  ReadLittleEndian(File, this->NbrBands);
  ReadLittleEndian(File, this->NbrStatePerBand);
  int HeaderSize = -1;
  ReadLittleEndian(File, HeaderSize);
  if (HeaderSize >= (3 * sizeof(int) + 4 * sizeof(double)))
    {
      int TmpDimension = -1;
      ReadLittleEndian(File, TmpDimension);
      if (TmpDimension >= 2)
	{
	  ReadLittleEndian(File, this->NbrSiteX);
	  ReadLittleEndian(File, this->KxFactor);
	  ReadLittleEndian(File, this->GammaX);	 
	  ReadLittleEndian(File, this->NbrSiteY);
	  ReadLittleEndian(File, this->KyFactor);
	  ReadLittleEndian(File, this->GammaY);	  
	}
      HeaderSize -= (3 * sizeof(int) + 4 * sizeof(double));
      if (HeaderSize > 0) 
	File.seekg (HeaderSize, ios::cur);
    }
  else
    {
      this->NbrSiteX = this->NbrStatePerBand;
      this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
      this->GammaX = 0.0;
      this->NbrSiteY = 1;
      this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
      this->GammaY = 0.0;
      File.seekg (HeaderSize, ios::cur);
    }
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
      for (int j = 0; j < this->NbrStatePerBand; ++j)
	{
	  ReadLittleEndian(File, this->EnergyBandStructure[i][j]);
	}
    }
  if (FileSize == ((sizeof(double) * this->NbrStatePerBand * this->NbrBands) + sizeof(long) + sizeof(int) + sizeof(int) + HeaderSize))
    {
      this->OneBodyBasis = 0;
    }
  else
    {
      this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
      for (int j = 0; j < this->NbrStatePerBand; ++j)	
	{
	  this->OneBodyBasis[j].ReadMatrix(File);
	}     
    }
  File.close();
}

// destructor
//

Generic2DTightBindingModel::~Generic2DTightBindingModel()
{
}
