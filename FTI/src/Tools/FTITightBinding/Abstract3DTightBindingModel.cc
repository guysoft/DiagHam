////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 3D tight binding model                  //
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
#include "Tools/FTITightBinding/Abstract3DTightBindingModel.h"

#include <fstream>


using std::ofstream;
using std::endl;


// default constructor
//

Abstract3DTightBindingModel::Abstract3DTightBindingModel()
{
}

// destructor
//

Abstract3DTightBindingModel::~Abstract3DTightBindingModel()
{
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool Abstract3DTightBindingModel::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky    kz";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  for (int kz = 0; kz < this->NbrSiteZ; ++kz)
	    {
	      int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(kx, ky, kz);
	      File << kx << " " << ky << " " << kz; 
	      for (int i = 0; i < this->NbrBands; ++i)
		File << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex];
	      File << endl;
	    }
	}
    }
  File.close();
  return true;
}
