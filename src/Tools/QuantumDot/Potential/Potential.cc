////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                           base class for potential                         //
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


#include "Tools/QuantumDot/Potential/Potential.h"
#include "BitmapPicture/BmpFormat.h"
#include "Color/PicRGB.h"

#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>

using std::ifstream;
using std::cout;
using std::endl;

// default constructor
//

Potential::Potential():NumberX(0),NumberY(0),NumberZ(0)
{
  this->Alloy = 0;
}

// constructor from a 3D dimension
// 
// x, y, z: three dimensions in respective direction

Potential::Potential(int x, int y, int z):NumberX(x),NumberY(y),NumberZ(z)
{
  this->Alloy = new bool** [z];
  for (int k = 0; k < z; ++k)
    {
      this->Alloy[k] = new bool* [y];
	for (int j = 0; j < y; ++j)
	  this->Alloy[k][j] = new bool [x];
    }
}

// virtual destructor
//

Potential::~Potential()
{
  delete[] this->Alloy;
}

// print the potential dimension
//
// output = an ostream to display
// display Nx, Ny and Nz in an output

ostream& Potential::PrintDimension(ostream& output)
{
  output << "Potential dimension is: " << this->NumberX << ", " << this->NumberY << ", " << this->NumberZ << " in x, y, z respectively.\n";
  return output;
}

// print the diagram of atomic distribution in an output
//
// output = an ostream to display
// display Alloy arry as z blocks, each block contains x columns and y lines

ostream& Potential::PrintDiagram(ostream& output)
{
  for (int k = 0; k < this->NumberZ; ++k)
    {
      for (int j = 0; j < this->NumberY; ++j)
	{
	  for (int i = 0; i < this->NumberX; ++i)
	    output << this->Alloy[k][j][i] << " ";
	  output << '\n';
	}
      output << '\n';
    }
  return output;
}

// read the atom diagram from an input, which has z blocks each block contains x columns and y lines
//
// fileName = name of the storage file

void Potential::ReadDiagram(char* fileName)
{
  ifstream input(fileName);
  if (! input.is_open())
    {
      cout << "Error when open the diagram file: " << fileName << " . Exit now" << endl;
      exit(1);
    }
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	input >> this->Alloy[k][j][i];
  input.close();
}

// generate an arbitrary distribution of atoms with a given proportion
//
// proportion = nominal proportion

void Potential::ArbitraryDistribution(double proportion)
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

// calculate the mean first neighbors in a potential diagram
//
// Total: reference to the total number of InN 
// return = the mean number of first neighbors (6)

double Potential::MeanFirstNeighbors(double& Total)
{
  Total = 0.0; int Pairs = 0;

  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	if (this->Alloy[k][j][i])
	  {
	    ++Total;
	    if (i > 0)
	      Pairs += this->Alloy[k][j][i - 1];
	    if (i < (this->NumberX - 1))
	      Pairs += this->Alloy[k][j][i + 1];
	    if (j > 0)
	      Pairs += this->Alloy[k][j - 1][i];
	    if (j < (this->NumberY - 1))
	      Pairs += this->Alloy[k][j + 1][i];
	    if (k > 0)
	      Pairs += this->Alloy[k - 1][j][i];
	    if (k < (this->NumberZ - 1))
	      Pairs += this->Alloy[k + 1][j][i];
	  }
  return double(Pairs) / double(Total);      
}

// save the diagram to a bitmap file
//
// u: the number of monolayers omitted under the considered structure
// a: the number of monolayers omitted above the considered structure
// startX: point to start in diagram in X direction
// endX: point to end in diagram in X direction
// startY: point to start in diagram in Y direction
// endY: point to end in diagram in Y direction
// choice: choice of orientation: 1 for XY, 2 for XZ, 3 for YZ
// sizeX: the size in X direction of each cell in pixel
// sizeY: the size in Y direction of each cell in pixel  
// InN:   color (in RGB definition) of InN cell
// GaN:   color (in RGB definition) of GaN cell
// background : color (in RGB definition) of background
// NbrX: number of cell displayed in X direction
// fileName: name of the file to store the picture

bool Potential::SaveBmpPicture(int u, int a, int startX, int endX, int startY, int endY, int choice, int sizeX, int sizeY, PicRGB& InN, PicRGB& GaN, PicRGB& background, int NbrX, char* fileName)
{
  int Layers = 0;
  switch(choice)
    {
    case 1: 
      Layers = this->NumberZ - u - a;
      break;
    case 2: 		  
      Layers = this->NumberY - u - a;
      break;
    default:
      Layers = this->NumberX - u - a;
    }
  int NbrY = (Layers - 1) / NbrX + 1;
  int Border = 5;
  int tmpLx = (endX - startX) * sizeX + Border;
  int tmpLy = (endY - startY) * sizeY + Border;
  int Lx = tmpLx * NbrX - Border;
  int Ly = tmpLy * NbrY - Border;

  BmpFormat* picture = new BmpFormat(Lx, Ly, background);

  int tmpX, tmpY, tmpZ, originX, originY, tmpOriginX;
  bool dopage;
  for (int k = 0; k < Layers; ++k)
    {
      tmpY = k / NbrX;
      tmpX = k - tmpY * NbrX;
      tmpZ = k + u;
      originX = tmpX * tmpLx;
      originY = tmpY * tmpLy;

      for (int j = startY; j < endY; ++j)
	{	
	  tmpOriginX = originX;
	  for (int i = startX; i < endX; ++i)
	    {
	      switch(choice)
		{
		case 1: 
		  dopage = this->Alloy[tmpZ][j][i];
		  break;
		case 2: 		  
		  dopage = this->Alloy[j][tmpZ][i];
		default:
		  dopage = this->Alloy[j][i][tmpZ];
		}
	      if (dopage)
		CellFill(tmpOriginX, sizeX, originY, sizeY, InN, picture);
	      else
		CellFill(tmpOriginX, sizeX, originY, sizeY, GaN, picture);
      
	      tmpOriginX += sizeX;
	    }
	  originY += sizeY;
	}
    }
  
  picture->SavePicture(fileName);

  return true;
}

// fill a cell with a given color
//
// startX: X coordination to start
// sizeX:  size in X direction
// startY: Y coordination to start
// sizeY:  size in Y direction
// Col: color in RGB definition
// picture: pointer to picture which is filled

void CellFill(int startX, int sizeX, int startY, int sizeY, PicRGB& Col, AbstractBitmapPicture* picture)
{
  for (int i = startX; i < startX + sizeX; ++i)
    for (int j = startY; j < startY + sizeY; ++j)
      picture->SetPixel(i, j, Col);
}
