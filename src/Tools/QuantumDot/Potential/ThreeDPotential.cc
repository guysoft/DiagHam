////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                            class for 3D potential                          //
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


#include "Tools/QuantumDot/Potential/ThreeDPotential.h"

#include <time.h>
#include <stdlib.h>
#include <math.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ofstream;

// default constructor
//
ThreeDPotential::ThreeDPotential()
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

ThreeDPotential::ThreeDPotential(int x, int y, int z)
{
  this->NumberX = x;
  this->NumberY = y;
  this->NumberZ = z;
  this->Alloy = new bool** [z];
  this->Potential = new double** [z];
  for (int k = 0; k < z; ++k)
    {
      this->Alloy[k] = new bool* [y];
      this->Potential[k] = new double* [y];
      for (int j = 0; j < y; ++j)	
	{  
	  this->Alloy[k][j] = new bool [x];
	  this->Potential[k][j] = new double[x];
	}
    }
}

// construct a frame potential with more geometric parameter
//
// u and a: number of mono-layers under and above the structure respectively

ThreeDPotential::ThreeDPotential(int x, int y, int z, int u, int a)
{
  this->NumberX = x;
  this->NumberY = y;
  this->NumberZ = z;
  this->under = u;
  this->above = a;
  this->Alloy = new bool** [z];
  this->Potential = new double** [z - u - a];
  this->UnderSize = new double[u];
  this->UnderPotential = new double[u];
  this->AboveSize = new double[a];
  this->AbovePotential = new double[a];
  for (int k = 0; k < z - u - a; ++k)
    {
      this->Potential[k] = new double* [y];
	for (int j = 0; j < y; ++j)
	  this->Potential[k][j] = new double [x];	  
    }
  for (int k = 0; k < z; ++k)
    {
      this->Alloy[k] = new bool* [y];
      for (int j = 0; j < y; ++j)	  
	this->Alloy[k][j] = new bool [x];
    }
  this->Concentration = 0;
}

// destructor
//
ThreeDPotential::~ThreeDPotential()
{
  delete[] this->Potential;
  delete[] this->Alloy;
}

// construct a well with arbitrary distribution of alloy in the well
//
// down_field = electric field under the considered structure
// field = electric field in the considered structure
// up_field = electric field above the considered structure
// cell = the height of a monolayer
// proportion = proportion of alloy in the considered structure
// offset = offset between host material and alloy
// scratch = true if construct from scratch
// fileName = name of the file to store parameters
// 0 reference: the monolayer just before the downside of the well

void ThreeDPotential::UniformWell(double down_field, double field, double up_field, double cell, double proportion, double offset, bool scratch, char* fileName = 0)
{
  if (fileName != 0)
    {
      ofstream File;
      File.open(fileName, ios::out | ios::app);
      File << "Uniform quantum well parameters\n";
      File << "Dimension (cells): X = " << this->NumberX << ", Y = " << this->NumberY << ", Z = " << this->NumberZ << '\n';
      File << "Number of cells under and above the well: Under = " << this->under << ", Above = " << this->above << '\n';
      File << "Proportion of InN: " << proportion << '\n';
      File << "Electric fields under, in and above the well: " << down_field << ", " << field << ", " << up_field << '\n';
      File << "Band offset: " << offset << '\n';
      File << "Growth direction lattice constant: " << cell << endl;
      File.close();
    }

  srand(time(NULL));
  if (scratch)
    {
      for (int k = 0; k < this->under; ++k)
	for (int j = 0; j < this->NumberY; ++j)
	  for (int i = 0; i < this->NumberX; ++i)
	    this->Alloy[k][j][i] = false;

      for (int k = this->under; k < ((this->NumberZ) - this->above); ++k)
	for (int j = 0; j < this->NumberY; ++j)
	  for (int i = 0; i < this->NumberX; ++i)
	    {
	      if ((double(rand())/RAND_MAX) < proportion)
		this->Alloy[k][j][i] = true;
	      else
		this->Alloy[k][j][i] = false;
	    }

      for (int k = (this->NumberZ) - (this->above); k < this->NumberZ; ++k)
	for (int j = 0; j < this->NumberY; ++j)
	  for (int i = 0; i < this->NumberX; ++i)
	    this->Alloy[k][j][i] = false;
    }
  double tmp = 0.0;
  for (int k = 0; k < this->under; ++k)
    {
      this->UnderPotential[k] = down_field * (k + 1 - this->under) * cell;
      this->UnderSize[k] = cell;
    }
  double tmp2 = 0.0; int tmpInd = 0;
  for (int k = this->under; k < ((this->NumberZ) - this->above); ++k)
    {
      tmp = field * (k + 1 - this->under) * cell;
      tmp2 = tmp - offset;
      tmpInd = k - this->under;
      for (int j = 0; j < this->NumberY; ++j)
	for (int i = 0; i < this->NumberX; ++i)
	  {
	    if (this->Alloy[k][j][i])
	      this->Potential[tmpInd][j][i] = tmp2;
	    else
	      this->Potential[tmpInd][j][i] = tmp;
	  }
    }
  tmpInd = NumberZ - above - under;
  for (int k = 0; k < this->above; ++k)
    {
      this->AbovePotential[k] = field * tmpInd + up_field * (k + 1) * cell;
      this->AboveSize[k] = cell;
    }
}

// construct a well with a pocket alloy upon the well
//
// height = height of the pocket
// width = width of the pocket
// down_field = electric field under the considered structure
// field = electric field in the considered structure
// up_field = electric field above the considered structure
// cell = the height of a monolayer
// offset = offset between host material and alloy
// scratch = true if construct from scratch
// fileName = name of the file to store parameters
// 0 reference: the monolayer just before the downside of the well
// ATTENTION TO Resize of Potential

void ThreeDPotential::FluctuatedWell(int height, int width, double down_field, double field, double up_field, double cell, double offset, char* fileName = 0)
{
  if (fileName != 0)
    {
      ofstream File;
      File.open(fileName, ios::out | ios::app);
      File << "Width fluctuated quantum well parameters\n";
      File << "Dimension (cells): X = " << this->NumberX << ", Y = " << this->NumberY << ", Z = " << this->NumberZ << '\n';
      File << "Number of cells under and above the well: Under = " << this->under << ", Above = " << this->above << '\n';
      File << "Height and width of fluctuation " << height << ", " << width << '\n';
      File << "Electric fields under, in and above the well: " << down_field << ", " << field << ", " << up_field << '\n';
      File << "Band offset: " << offset << '\n';
      File << "Growth direction lattice constant: " << cell << endl;
      File.close();
    }
  delete[] this->Potential;
  this->Potential = 0;
  this->Potential = new double ** [NumberZ - under - above + height];
  for (int k = 0; k < NumberZ - under - above + height; ++k)
    {
      this->Potential[k] = new double* [this->NumberY];
      for (int j = 0; j < this->NumberY; ++j)
	this->Potential[k][j] = new double[this->NumberX];
    }
  double tmp = 0.0; double tmp1 = 0.0;
  for (int k = 0; k < this->under; ++k)
    {
      tmp = ((k - this->under) * down_field - (this->NumberZ - this->under - this->above) * field) * cell;
      for (int j = 0; j < this->NumberY; ++j)
	for (int i = 0; i < this->NumberX; ++i)
	  {
	    this->Alloy[k][j][i] = false;
	    this->UnderSize[k] = cell;
	    this->UnderPotential[k] = tmp;
	  }
    }
  int tmpInd = 0;
  for (int k = this->under; k < this->NumberZ - this->above; ++k)
    {
      tmpInd = k - this->under;
      tmp = cell * (k + this->above - this->NumberZ) * field - offset;
      for (int j = 0; j < this->NumberY; ++j)
	for (int i = 0; i < this->NumberX; ++i)
	  {
	    this->Alloy[k][j][i] = true;
	    this->Potential[tmpInd][j][i] = tmp;
	  }
    }
  int X1 = (this->NumberX - width)/2; int X2 = (this->NumberX + width)/2;
  int Y1 = (this->NumberY - width)/2; int Y2 = (this->NumberY + width)/2;
  for (int k = this->NumberZ - this->above; k < this->NumberZ - this->above + height; ++k)
    {
      tmp = (k - this->NumberZ + this->above) * up_field * cell;
      tmp1 = tmp - offset;
      tmpInd = k - this->under;
      for (int j = 0; j < Y1; ++j)
	for (int i = 0; i < this->NumberX; ++i)
	  {
	    this->Alloy[k][j][i] = false;
	    this->Potential[tmpInd][j][i] = tmp;
	  }
      for (int j = Y1; j < Y2; ++j)
	{
	  for (int i = 0; i < X1; ++i)
	    {
	      this->Alloy[k][j][i] = false;
	      this->Potential[tmpInd][j][i] = tmp;
	    }
	  for (int i = X1; i < X2; ++i)
	    {
	      this->Alloy[k][j][i] = true;	  
	      this->Potential[tmpInd][j][i] = tmp1;
	    }
	  for (int i = X2; i < this->NumberX; ++i)
	    {
	      this->Alloy[k][j][i] = false;
	      this->Potential[tmpInd][j][i] = tmp;
	    }
	}
      for (int j = Y2; j < this->NumberY; ++j)
	for (int i = 0; i < this->NumberX; ++i)
	  {
	    this->Alloy[k][j][i] = false;      
	    this->Potential[tmpInd][j][i] = tmp;
	  }
    }
  for (int k = this->NumberZ - this->above + height; k < this->NumberZ; ++k)
    {
      tmp = (k - this->NumberZ + this->above) * up_field * cell;
      tmpInd = k - this->NumberZ + this->above - height;
      for (int j = 0; j < this->NumberY; ++j)
	for (int i = 0; i < this->NumberX; ++i)
	  {
	    this->Alloy[k][j][i] = false;
	    this->AboveSize[tmpInd] = cell;
	    this->AbovePotential[tmpInd] = tmp;
	  }
    }
  this->above -= height;
}

// potential 0 reference: upside monolayer
// ATTENTION TO Resize of Potential, UnderSize, UnderPotential, AboveSize & AbovePotential
void ThreeDPotential::IEFFluctuatedWell(double p, double down_field, double field, double up_field, double cell, double offset, char* fileName = 0)
{
  srand(time(NULL));

  if (this->Potential != 0)
    delete[] this->Potential;
  if (this->UnderSize != 0)
    delete[] this->UnderSize;
  if (this->UnderPotential != 0)
    delete[] this->UnderPotential;
  if (this->AboveSize != 0)
    delete[] this->AboveSize;
  if (this->AbovePotential != 0)
    delete[] this->AbovePotential;

  int oldUnder = this->under;
  int oldAbove = this->above;
  int wellWidth = this->NumberZ - this->under - this->above;

  this->above -= 1;
  this->under = this->NumberZ - this->above - 2;
  this->UnderSize = new double[this->under];
  this->UnderPotential = new double[this->under];
  this->AboveSize = new double[this->above];
  this->AbovePotential = new double[this->above];
  this->Potential = new double**[2];
  for (int k = 0; k < 2 ; ++k)
    {
      this->Potential[k] = new double*[this->NumberY];
      for (int j = 0; j < this->NumberY; ++j)	
	this->Potential[k][j] = new double[this->NumberX];	
    }

  for (int k = 0; k < 2; ++k)    
    for (int j = 0; j < this->NumberY; ++j)
      {
	if ((double(rand())/RAND_MAX) < p)
	  this->Alloy[k + this->under][j][0] = true;
	else
	  this->Alloy[k + this->under][j][0] = false;
	for (int i = 1; i < this->NumberX; ++i)
	  {
	    if ((double(rand())/RAND_MAX) < p)
	      this->Alloy[k + this->under][j][i] = this->Alloy[oldAbove - 1][j][i - 1];
	    else
	      this->Alloy[k + this->under][j][i] = !(this->Alloy[oldAbove - 1][j][i - 1]);
	  }
      }

  for (int k = 0; k < oldUnder; ++k)
    {
      this->UnderSize[k] = cell;
      this->UnderPotential[k] = (-field * wellWidth + down_field * (k - oldUnder)) * cell;    
      for (int j = 0; j < this->NumberY; ++j)
	for (int i = 0; i < this->NumberX; ++i)
	  this->Alloy[k][j][i] = false;  
    }
  for (int k = oldUnder; k < this->under ;++k)
    {
      this->UnderSize[k] = cell;
      this->UnderPotential[k] = field * (k - this->under - 1) * cell - offset;
      for (int j = 0; j < this->NumberY; ++j)
	for (int i = 0; i < this->NumberX; ++i)
	  this->Alloy[k][j][i] = true;  
    }
 
  for (int j = 0; j < this->NumberY; ++j)
    for (int i = 0; i < this->NumberX; ++i)
      {
	if (Alloy[this->under][j][i])
	  this->Potential[0][j][i] = -field * cell - offset;
	else
	  this->Potential[0][j][i] = -field * cell;;
	}
  for (int j = 0; j < this->NumberY; ++j)
    for (int i = 0; i < this->NumberX; ++i)
      {
	if (Alloy[this->under + 1][j][i])
	  this->Potential[1][j][i] = -offset;
	else
	  this->Potential[1][j][i] = 0.0;
      }
  for (int k = 0; k < this->above; ++k)
    {
      this->AboveSize[k] = cell;
      this->AbovePotential[k] = up_field * (k + 1) * cell;
      for (int j = 0; j < this->NumberY; ++j)
	for (int i = 0; i < this->NumberX; ++i)
	  this->Alloy[k + this->NumberZ - this->above][j][i] = false;  	       
    }
}

// construct a truncated hexagonal pyramid dot used for nitride materials simulation with uniform distribution
//
// proportion = proportion of alloy materials
// Rb = base radius, Rt = top radius
// w = thickness of wetting layer
// down_field, wetting_field, dot_field, up_field = the electric field in respective part
// offset = offset bulk-alloy material
// c: height of a monolayer
// scratch = true if constructed from nothing, else from existing Diagram
// fileName = file to store parameters

void ThreeDPotential::ArbitraryPyramidDot(double proportion, int Rb, int Rt, int w, double down_field, double wetting_field, double dot_field, double up_field, double offset, double c, bool scratch, char* fileName = 0)
{
  srand(time(NULL));
  if (fileName != 0)
    {
      ofstream File;
      File.open(fileName, ios::out | ios::app);
      File << "Pyramid quantum dot parameters\n";
      File << "Dimension (cells): X = " << this->NumberX << ", Y = " << this->NumberY << ", Z = " << this->NumberZ << '\n';
      File << "Number of cells under and above the dot: Under = " << this->under << ", Above = " << this->above << '\n';
      File << "Proportion of InN: " << proportion << '\n';
      File << "Base and top radius: Base = " << Rb << ", Top = " << Rt << '\n';
      File << "Wetting layer thickness: " << w << '\n';
      File << "Electric fields under the dot, in the wetting layers, in and above the dot: " << down_field << ", " << wetting_field << ", " << dot_field << ", " << up_field << '\n';
      File << "Band offset: " << offset << '\n';
      File << "Growth direction lattice constant: " << c << endl;
      File.close();
    }
  //position of base center
  int Cx = this->NumberX/2;
  int Cy = this->NumberY/2;

  double tmp1 = 0.0, tmp2 = 0.0;

  // under the wetting layer
  for (int k = 0; k < under; ++k)
    {
      tmp1 = down_field * (k + 1 - under) * c;
      if (scratch)
        for (int j = 0; j < NumberY; ++j)
	  for (int i = 0; i < NumberX; ++i)
	    this->Alloy[k][j][i] = false;

      this->UnderPotential[k] = tmp1;
      this->UnderSize[k] = c;
    }

  // wetting layer
  for (int k = under; k < under + w; ++k)
    {
      tmp1 = (k + 1 - under) * wetting_field * c;
      tmp2 = tmp1 - offset;
      for (int j = 0; j < NumberY; ++j)
	for (int i = 0; i < NumberX; ++i)
	  {
	    if (scratch)
	      if ((double(rand())/RAND_MAX) < proportion)	      
		Alloy[k][j][i] = true;	      
	      else	      
		Alloy[k][j][i] = false;	      
	    if (Alloy[k][j][i])
	      Potential[k - under][j][i] = tmp2;
	    else
	      Potential[k - under][j][i] = tmp1;	      
	  }
    }

  // quantum dot
  double t1 = 0.0; double in1 = 0;
  double t2 = 0.0; double in2 = 0;
  double Rk = 0.0;
  bool base[NumberY][NumberX]; bool interface[NumberY][NumberX];

  for (int j = 0; j < NumberY; ++j)
    for (int i = 0; i < NumberX; ++i)
      {
	interface[j][i] = false;
	base[j][i] = false;
      }

  // the first mono-layer of the quantum dot
  tmp1 = (wetting_field * w + dot_field) * c; tmp2 = tmp1 - offset; double tmp3 = wetting_field * w * c;
  for (int j = 0; j < NumberY; ++j)
    for (int i = 0; i < NumberX; ++i)
      {
	t1 = fabs((double)(i - Cx));
	t2 = fabs((double)(j - Cy));
	in1 = t1 - Rb + t2/sqrt(3.0);
	in2 = t2 - Rb*sqrt(3.0)/2.0;
	// in the dot
	if ((in1 < 0) && (in2 < 0))
	  {
	    if (scratch)
	      if ((double(rand())/RAND_MAX) < proportion)	      
		Alloy[under + w][j][i] = true;	      
	      else	      
		Alloy[under + w][j][i] = false;

	    if (Alloy[under + w][j][i])
	      Potential[w][j][i] = tmp2;
	    else
	      Potential[w][j][i] = tmp1;
	    base[j][i] = true;
	    interface[j][i] = true;
	  }
	// outside the dot
	else
	  {
	    Potential[w][j][i] = tmp3;
	    if (scratch)
	      Alloy[under + w][j][i] = false;
	  }
      }

  // the following mono-layers in the dot
  for (int k = under + w + 1; k < NumberZ - above; ++k)
    {
      tmp1 = (wetting_field * w + dot_field * (k + 1 - under - w)) * c; tmp2 = tmp1 - offset;
      Rk = Rb - (k - under - w + 1)* 1.0 * (Rb - Rt)/(NumberZ - w -under - above);
      for (int j = 0; j < NumberY; ++j)
	for (int i = 0; i < NumberX; ++i)
	  {
	    // outside the base
	    if (!base[j][i])
	      {
		if (scratch)
		  Alloy[k][j][i] = false;
		Potential[k - under][j][i] = tmp3;
	      }
	    // inside the base
	    else
	      {
		t1 = fabs((double)(i - Cx));
		t2 = fabs((double)(j - Cy));
		in1 = t1 - Rk + t2/sqrt(3.0);
		in2 = t2 - Rk * sqrt(3.0)/2.0;
		// in the dot
		if ((in1 < 0) && (in2 < 0))
		  {
		    if (scratch)		      
		      if ((double(rand())/RAND_MAX) < proportion)			
			Alloy[k][j][i] = true;		      
		      else		      
			Alloy[k][j][i] = false;		      

		    if (Alloy[k][j][i])
		      Potential[k - under][j][i] = tmp2;
		    else
		      Potential[k - under][j][i] = tmp1;
		  }
		// between the dot and the base
		else
		  {
		    if (interface[j][i])
		      {
			interface[j][i] = false;
			if (scratch)
			  if ((double(rand())/RAND_MAX) < proportion)			    
			    Alloy[k][j][i] = true;			  
			  else			  
			    Alloy[k][j][i] = false;
			
			if (Alloy[k][j][i])
			  Potential[k - under][j][i] = tmp2;
			else
			  Potential[k - under][j][i] = tmp1;
		      }
		    else
		      {
			Potential[k - under][j][i] = Potential[k - under - 1][j][i];
			if (scratch)
			  Alloy[k][j][i] = false;
		      }
		  }
	      }
	  }
    }

  // above the dot
  int tmpIndice = this->NumberZ - this->above;
  for (int k = 0; k < above; ++k)
    {
      tmp1 = (wetting_field * w + dot_field * (NumberZ - w - under - above) + up_field * (k + 1)) * c ;
      if (scratch)
	for (int j = 0; j < NumberY; ++j)
	  for (int i = 0; i < NumberX; ++i)
	    Alloy[tmpIndice + k][j][i] = false;
      this->AbovePotential[k] = tmp1;
      this->AboveSize[k] = c;
    }
}

// print the 3-D potential in an output
//
// output = output (file, screen, ...)
// output has z block, each block has y lines and x columns

ostream& ThreeDPotential::PrintPotential(ostream& output)
{
  for (int k = 0; k < this->NumberZ; ++k)
    {
      //if ((k >= this->under) && (k < (this->NumberZ - this->above)))
	//cout << "Dot" << '\n';
      for (int j = 0; j < this->NumberY; ++j)
	{
	  for (int i = 0; i < this->NumberX; ++i)
	    output << this->Potential[k][j][i] << '\t';
	  output << '\n';
	}
      output << '\n';
    }
  return output;
}

// print the 3-D potential with potential array along growth direction in an output
//  
// output : file, screen, ...
// output has z block, each block has y lines and x columns

ostream& ThreeDPotential::PrintPotentialWithField(ostream& output)
{
  output << this->under << " ";
  for (int k = 0; k < this->under; ++k)
    output << this->UnderPotential[k] << " " << this->UnderSize[k] << " ";
  output << '\n';
  output << this->above << " ";
  for (int k = 0; k < this->above; ++k)
    output << this->AbovePotential[k] << " " << this->AboveSize[k] << " ";
  output << '\n';
  for (int k = 0; k < (this->NumberZ - this->above - this->under); ++k)
    {
      for (int j = 0; j < this->NumberY; ++j)
	{
	  for (int i = 0; i < this->NumberX; ++i)
	    output << this->Potential[k][j][i] << '\t';
	  output << '\n';
	}
      output << '\n';
    }
  return output; 
}

// construct a potential from a file
//
// input = data from a file
// convention: z blocks, each block has y lines and x columns
// order of potential: zyx

void ThreeDPotential::ReadPotential(char* fileName)
{
  ifstream input(fileName);
  if (! input.is_open())
    {
      cout << "Error when open the potential file. Exit now" << endl;
      exit(1);
    }
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	input >> this->Potential[k][j][i];
}

// construct a potential from a file with electric field
//
// fileName = name of the storage file
// convention: z blocks, each block has y lines and x columns
// order of potential: zyx

void ThreeDPotential::ReadPotentialWithField(char* fileName)
{
  ifstream input(fileName);
  if (! input.is_open())
    {
      cout << "Error when open the potential file. Exit now" << endl;
      exit(1);
    }
  input >> this->under;
  for (int k = 0; k < this->under; ++k)
    input >> this->UnderPotential[k] >> this->UnderSize[k];
  input >> this->above;
  for (int k = 0; k < this->above; ++k)
    input >> this->AbovePotential[k] >> this->AboveSize[k];

  for (int k = 0; k < (this->NumberZ - this->above - this->under); ++k)    
    for (int j = 0; j < this->NumberY; ++j)      
      for (int i = 0; i < this->NumberX; ++i)
	input >> this->Potential[k][j][i];
  input.close();
}
