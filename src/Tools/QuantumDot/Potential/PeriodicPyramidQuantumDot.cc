////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                    class for periodic pyramid quantum dot                  //
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

#include "Tools/QuantumDot/Potential/PeriodicPyramidQuantumDot.h"

#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <math.h>

using std::cout;
using std::endl;


// constructor
//
// scratch = true if construct from scratch

PeriodicPyramidQuantumDotPotential::PeriodicPyramidQuantumDotPotential(int NbrCellX, int NbrCellY, int NbrCellZ, double Lz, int u, int a, int rb, int rt, int w, double offset, double concentration, double piezofield, bool scratch, char* logfile)
{
  this->NumberX = NbrCellX;
  this->NumberY = NbrCellY;
  this->NumberZ = NbrCellZ;
  this->Alloy = new bool** [NbrCellZ];
  this->Potential = new double** [NbrCellZ];
  for (int k = 0; k < NbrCellZ; ++k)
    {
      this->Alloy[k] = new bool* [NbrCellY];
      this->Potential[k] = new double* [NbrCellY];
      for (int j = 0; j < NbrCellY; ++j)	
	{  
	  this->Alloy[k][j] = new bool [NbrCellX];
	  this->Potential[k][j] = new double[NbrCellX];
	}
    }
  this->under = u;
  this->above = a;
  this->BaseRadius = rb;
  this->TopRadius = rt;
  this->WettingWidth = w;
  if (!this->GeneratePotential(concentration, piezofield, Lz, offset, scratch, logfile))
    cout << "Potential cannot be generated" << endl;
}

// destructor
// 

PeriodicPyramidQuantumDotPotential::~PeriodicPyramidQuantumDotPotential()
{

}

bool PeriodicPyramidQuantumDotPotential::GeneratePotential(double C, double F, double Lz, double offset, bool scratch, char* filename)
{
  srand(time(NULL));
  int** Height = new int* [this->NumberX];
  int temp;
  for (int j = 0; j < this->NumberY; ++j)
    {
      Height[j] = new int [this->NumberX];
      for (int i = 0; i < this->NumberX; ++i)
	{
	  temp = 0;
	  for (int k = this->under + this->WettingWidth; this->InTheDot(i, j, k); ++k)
	    ++temp;
	  Height[j][i] = temp + this->WettingWidth;	  
	}
    }
  
  double** Barrier = new double* [this->NumberY];
  double** Dot = new double* [this->NumberY];

  double*** ReferencePotential = new double** [this->NumberZ];
  for (int k = 0; k < this->NumberZ; ++k)
    {
      ReferencePotential[k] = new double* [this->NumberY];
      for (int j = 0; j < this->NumberY; ++j)	
	ReferencePotential[k][j] = new double[this->NumberX];
    }

  for (int j = 0; j < this->NumberY; ++j)
    {
      Barrier[j] = new double [this->NumberX];
      Dot[j] = new double [this->NumberX];
      for (int i = 0; i < this->NumberX; ++i)
	{
	  Barrier[j][i] = F * Height[j][i] / double(this->NumberZ);
	  Dot[j][i] = Barrier[j][i] - F;
	  Barrier[j][i] *= Lz;
	  Dot[j][i] *= Lz;	  
	}
    }

  for (int j = 0; j < this->NumberY; ++j)
    for (int i = 0; i < this->NumberX; ++i)
      {
	if (scratch)
	  this->Alloy[0][j][i] = false;
	ReferencePotential[0][j][i] = 0.0;
	this->Potential[0][j][i] = 0.0;
      }

  for (int k = 1; k < this->under; ++k)    
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)  
	{ 
	  if (scratch)
	    this->Alloy[k][j][i] = false; 	  
	  ReferencePotential[k][j][i] = ReferencePotential[k - 1][j][i] + Barrier[j][i];  
	  this->Potential[k][j][i] = ReferencePotential[k - 1][j][i] + Barrier[j][i];         
	}
  
  double TmpPotential;
  for (int k = this->under; k < (this->under + this->WettingWidth); ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)  
	  { 
	    ReferencePotential[k][j][i] = ReferencePotential[k - 1][j][i] + Dot[j][i];  
	    TmpPotential = ReferencePotential[k - 1][j][i] + Dot[j][i];
	    if (scratch)	    
	      if ((double(rand())/RAND_MAX) < C)	      	
		this->Alloy[k][j][i] = true;  	      
	      else	      	
		this->Alloy[k][j][i] = false;  	      
	    if (this->Alloy[k][j][i])
	      this->Potential[k][j][i] = TmpPotential - offset;
	    else
	      this->Potential[k][j][i] = TmpPotential;
	  }

  for (int k = this->under + this->WettingWidth; k < (this->NumberZ - this->above); ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i) 
	{
	  if (InTheDot(i, j, k))
	    {
	      ReferencePotential[k][j][i] = ReferencePotential[k - 1][j][i] + Dot[j][i];  
	      TmpPotential = ReferencePotential[k - 1][j][i] + Dot[j][i];	      
	      if (scratch)	    
		if ((double(rand())/RAND_MAX) < C)	      	
		  this->Alloy[k][j][i] = true;  	      
		else	      	
		  this->Alloy[k][j][i] = false;  	      
	      if (this->Alloy[k][j][i])
		this->Potential[k][j][i] = TmpPotential - offset;
	      else
		this->Potential[k][j][i] = TmpPotential; 
	    }	      	
	  else
	    {
	      ReferencePotential[k][j][i] = ReferencePotential[k - 1][j][i] + Barrier[j][i];  
	      TmpPotential = ReferencePotential[k - 1][j][i] + Barrier[j][i];	   
	      if (scratch)
		this->Alloy[k][j][i] = false; 
	      this->Potential[k][j][i] = TmpPotential;	
 	      
	    }
	}

  for (int k = this->NumberZ - this->above; k < this->NumberZ; ++k)    
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)  
	{ 
	  if (scratch)	
	    this->Alloy[k][j][i] = false;  
	  ReferencePotential[k][j][i] = ReferencePotential[k - 1][j][i] + Barrier[j][i];  
	  TmpPotential = ReferencePotential[k - 1][j][i] + Barrier[j][i];	   
	  this->Potential[k][j][i] = TmpPotential;     
	}
  
  delete[] ReferencePotential; delete[] Barrier; delete[] Dot; delete[] Height;
  return true;
}


// determine if a cell is in the dot (wetting layer is included)
// x: x coordinate of the cell
// y: y coordinate of the cell
// z: z coordinate of the cell
//
// return : true if the cell is in the dot, false otherwise
bool PeriodicPyramidQuantumDotPotential::InTheDot(int x, int y, int z)
{
  //position of base center
  int Cx = this->NumberX/2;
  int Cy = this->NumberY/2;

  if ((z >= this->under) && (z < (this->under + this->WettingWidth)))
    return true; // in the wetting layer
  
  if (z >= (this->under + this->WettingWidth) && (z < (this->NumberZ - this->above)))
    {
      double t1 = fabs((double)(x - Cx));
      double t2 = fabs((double)(y - Cy));
      double Rz = double(this->BaseRadius) - double(z - this->under - this->WettingWidth + 1) * double(this->BaseRadius - this->TopRadius) / double(this->NumberZ - this->WettingWidth - this->under - this->above);
      double in1 = t1 - Rz + t2 / sqrt(3.0);
      double in2 = t2 - Rz * sqrt(3.0) / 2.0;
      // in the dot
      if ((in1 < 0) && (in2 < 0))
	return true;
      else
	return false;
    }
  return false;
}
