////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                            class for 2D potential                          //
//                                                                            //
//                        last modification : 09/15/2003                      //
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


#include "Tools/QuantumDot/Potential/TwoDPotential.h"

#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>

using std::ostream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;

// default constructor
//

TwoDPotential::TwoDPotential()
{
  this->Potential = 0;
  this->NumberX = 0;
  this->NumberY = 0;
  this->NumberZ = 0;
  this->Alloy = 0;
}

// construct a frame potential
//
// x, y, z: three dimensions in respective direction

TwoDPotential::TwoDPotential(int x, int y, int z)
{
  this->Potential = new double*[y];
  for (int j = 0; j < y; ++j)
    this->Potential[j] = new double[x];
  this->NumberX = x;
  this->NumberY = y;
  this->NumberZ = z;
  this->Alloy = new bool** [z];
  for (int k = 0; k < z; ++k)
    {
      this->Alloy[k] = new bool* [y];
	for (int j = 0; j < y; ++j)
	  this->Alloy[k][j] = new bool [x];
    }
}  

// destructor
//

TwoDPotential::~TwoDPotential()
{
  delete[] this->Potential;
  delete[] this->Alloy;
}

// read 2D potential from a file
//
// FileName = name of the data file
// convention: x columns and y lines 

void TwoDPotential::ReadTwoDPotential(char* FileName)
{
  ifstream InputFile;
  InputFile.open(FileName, ios::in);
  if (!InputFile.is_open())
    {
      cout << "Error in reading the file: " << FileName << ". Exit now" << endl;
      exit(0);
    }
  for (int j = 0; j < NumberY; ++j)
    for (int i = 0; i < NumberX; ++i)
      InputFile >> this->Potential[j][i];
  InputFile.close();

}

// read 2D potential from a file when the number of lines and columns are known
//
// FileName = nam of the data file
// convention: x columns and y lines

void TwoDPotential::ReadTwoDPotential(char* FileName, int NX, int NY, int NZ)
{
  this->NumberX = NX;
  this->NumberY = NY;
  this->NumberZ = NZ;
  this->Potential = new double* [NY];
  for (int i = 0; i < NX; ++i)
    this->Potential[i] = new double [NX];

  ifstream InputFile;
  InputFile.open(FileName, ios::in);
  if (!InputFile.is_open())
    {
      cout << "Error in reading the file: " << FileName << ". Exit now" << endl;
      exit(0);
    }

  for (int j = 0; j < NY; ++j)
    for (int i = 0; i < NX; ++i)
      InputFile >> this->Potential[j][i];
  InputFile.close();
}

// print the potential in a output
//
// output = output (file, screen, ...)
// return: the ostream
// convention: there are x columns and y lines

ostream& TwoDPotential::PrintPotential(ostream& output)
{
  for (int j = 0; j < this->NumberY; ++j)
    {
      for (int i = 0; i < this->NumberX; ++i)
	output << this->Potential[j][i] << " ";
      output << '\n';
    }
  return output;
}

// construct a well potential with uniform probability
//
// proportion = proportion of the alloy
// offset = offset of the alloy and the host meterials, withour virtual crystal approximation
// weight = weight of potential in each mono-layer
// scratch = true if constructed from scratch, false if from existing diagram

void TwoDPotential::UniformWell(double proportion, double offset, double* weight, bool scratch)
{
  if (scratch)
    {
      srand(time(NULL));
      for (int k = 0; k < this->NumberZ; ++k)
	{
	  for (int j = 0; j < this->NumberY; ++j)
	    {
	      for (int i = 0; i < this->NumberX; ++i)
		{
		  if ((double(rand())/RAND_MAX) < proportion)
		    {
		      this->Alloy[k][j][i] = true;
		    }
		  else 
		    {
		      this->Alloy[k][j][i] = false;
		    }
		}  
	    }
	}
    }
  double tmp = 0.0; double pro1 = 1.0 - proportion;
  for (int i = 0; i < this->NumberX; ++i)
    for (int j = 0; j < this->NumberY;++j)
      {
	tmp = 0.0; 
	for (int k = 0; k < this->NumberZ; ++k)
	  if (this->Alloy[k][j][i])
	    tmp -= pro1 * weight[k];
	  else
	    tmp += proportion * weight[k];
	this->Potential[j][i] = tmp * offset;    
      }
}

// construct a dot potential with uniform probability
void TwoDPotential::UniformPyramidDot(double proportion)
{
  
}


