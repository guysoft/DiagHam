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


// construct a well with segragation distribution of alloy
// down_field = electric field under the considered structure
// field = electric field in the considered structure
// up_field = electric field above the considered structure
// cell = the height of a monolayer
// proportion1 = probability without InN neighborhood
// proportion2 = probability with InN neighborhood
// offset = offset between host material and alloy
// scratch = true if construct from scratch

void ThreeDPotential::SegregationWell (double down_field, double field, double up_field, double cell, double proportion1, double proportion2, double offset, bool scratch, char* fileName)
{
  int counter = 0;
  int counterbis = 0;
  int X1 = 100, X2 = 400, Y1 = 100, Y2 = 400, Z1 = 30, Z2 = 90;
  for (int k = 0; k < this->NumberZ; ++k)
    for (int j = 0; j < this->NumberY; ++j)
      for (int i = 0; i < this->NumberX; ++i)
	this->Alloy[k][j][i] = false;
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
	      if (this->NeighborhoodBis(i, j, k))
		{
		  if ((double(rand())/RAND_MAX) < proportion2)
		    {
		      this->Alloy[k][j][i] = true;
		      ++counter;
		      if ((i <= X2) && (i > X1) && (j <= Y2) && (j > Y1) && (k <= Z2) && (k > Z1))
			++counterbis;
		    }
		  else
		    this->Alloy[k][j][i] = false;
		}
	      else
		{
		  if ((double(rand())/RAND_MAX) < proportion1)
		    {
		      this->Alloy[k][j][i] = true;
		      ++counter;
		      if ((i <= X2) && (i > X1) && (j <= Y2) && (j > Y1) && (k <= Z2) && (k > Z1))
			++counterbis;
		    }
		  else
		    this->Alloy[k][j][i] = false;
		}
	    }

      for (int k = (this->NumberZ) - (this->above); k < this->NumberZ; ++k)
	for (int j = 0; j < this->NumberY; ++j)
	  for (int i = 0; i < this->NumberX; ++i)
	    this->Alloy[k][j][i] = false;
    }
  double tmpC = counterbis / (double(X2 - X1) * double(Y2 - Y1) * double(Z2 - Z1));
  cout << "Concentration := " << tmpC << endl;
  cout << "Total number := " << counterbis << endl;
  ofstream test;
  test.open("test.txt", ios::out | ios::app);
  test << proportion1 << '\t' << proportion2 << '\t' << tmpC << endl;

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
      this->AbovePotential[k] = (field * tmpInd + up_field * (k + 1)) * cell;
      this->AboveSize[k] = cell;
    }

  double percent = double(counter)/double(this->NumberX * this->NumberY * (this->NumberZ - this->under - this->above));

  if (fileName != 0)
    {
      ofstream File;
      File.open(fileName, ios::out | ios::app);
      File << "Uniform quantum well parameters\n";
      File << "Dimension (cells): X = " << this->NumberX << ", Y = " << this->NumberY << ", Z = " << this->NumberZ << '\n';
      File << "Number of cells under and above the well: Under = " << this->under << ", Above = " << this->above << '\n';
      File << "Probability given: " << proportion1 << ", " << proportion2 << '\n';
      File << "Real proportion: " << percent << '\n';
      File << "Electric fields under, in and above the well: " << down_field << ", " << field << ", " << up_field << '\n';
      File << "Band offset: " << offset << '\n';
      File << "Growth direction lattice constant: " << cell << endl;
      File.close();
    }
  cout << "Real proportion: " << percent << endl;
}

// construct a truncated hexagonal pyramid dot used for nitride materials simulation with uniform distribution
//
// proportion1 = probability without InN neighborhood
// proportion2 = probability with InN neighborhood 
// Rb = base radius, Rt = top radius
// w = thickness of wetting layer
// down_field, wetting_field, dot_field, up_field = the electric field in respective part
// offset = offset bulk-alloy material
// c: height of a monolayer
// scratch = true if constructed from nothing, else from existing Diagram
// fileName = file to store parameters

void ThreeDPotential::SegregationPyramidDot(double proportion1, double proportion2, int Rb, int Rt, int w, double down_field, double wetting_field, double dot_field, double up_field, double offset, double c, bool scratch, char* fileName = 0)
{
  srand(time(NULL));
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
  int counter = 0;
  // wetting layer
  for (int k = under; k < under + w; ++k)
    {
      tmp1 = (k + 1 - under) * wetting_field * c;
      tmp2 = tmp1 - offset;
      for (int j = 0; j < NumberY; ++j)
	for (int i = 0; i < NumberX; ++i)
	  {
	    if (scratch)
	      {
		if (this->NeighborhoodBis(i, j, k))
		  if ((double(rand())/RAND_MAX) < proportion2)	
		    {      
		      Alloy[k][j][i] = true;	      
		      ++counter;		      
		    }
		  else	      
		    Alloy[k][j][i] = false;	      
		else
		  if ((double(rand())/RAND_MAX) < proportion1)	      
		    {
		      Alloy[k][j][i] = true;	      
		      ++counter;
		    }
		  else	      
		    Alloy[k][j][i] = false;
	      }
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
	      {
		if (this->NeighborhoodBis(i, j, under + w))
		  if ((double(rand())/RAND_MAX) < proportion2)	      
		    {
		      Alloy[under + w][j][i] = true;	      
		      ++counter;
		    }
		  else	      
		    Alloy[under + w][j][i] = false;	      
		else
		  if ((double(rand())/RAND_MAX) < proportion1)	      
		    {
		      Alloy[under + w][j][i] = true;	      
		      ++counter;
		    }
		  else	      
		    Alloy[under + w][j][i] = false;
	      }

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
		      {
			if (this->NeighborhoodBis(i, j, k))
			  if ((double(rand())/RAND_MAX) < proportion2)			
			    {
			      Alloy[k][j][i] = true;		      
			      ++counter;
			    }
			  else		      
			    Alloy[k][j][i] = false;		      
			else			  
			  if ((double(rand())/RAND_MAX) < proportion1)			
			    {
			      Alloy[k][j][i] = true;		      
			      ++counter;
			    }
			  else		      
			    Alloy[k][j][i] = false;
		      }
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
			  {
			    if (this->NeighborhoodBis(i, j, k))
			      if ((double(rand())/RAND_MAX) < proportion2)			
				{
				  Alloy[k][j][i] = true;		      
				  ++counter;
				}
			      else		      
				Alloy[k][j][i] = false;		      
			    else			  
			      if ((double(rand())/RAND_MAX) < proportion1)			
				{
				  Alloy[k][j][i] = true;		      
				  ++counter;
				}
			      else		      
				Alloy[k][j][i] = false;
			  }
			
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

  int tmpIndice = this->NumberZ - this->above;
  // above the dot
  for (int k = 0; k < this->above; ++k)
    {
      tmp1 = (wetting_field * w + dot_field * (NumberZ - w - under - above) + up_field * (k + 1)) * c;
      if (scratch)
	for (int j = 0; j < NumberY; ++j)
	  for (int i = 0; i < NumberX; ++i)
	    Alloy[tmpIndice][j][i] = false;
      ++tmpIndice;
      this->AbovePotential[k] = tmp1;
      this->AboveSize[k] = c;

    }
  
  double dot_volume = this->NumberX * this->NumberY * w + (this->NumberZ - this->under - this->above - w) * (Rb * Rb + Rt * Rt) * 3.0 * sqrt(3.0)/4; 
  double percent =  double(counter)/dot_volume;
  this->Concentration = percent;
  if (fileName != 0)
    {
      ofstream File;
      File.open(fileName, ios::out | ios::app);
      File << "Pyramid quantum dot parameters\n";
      File << "Dimension (cells): X = " << this->NumberX << ", Y = " << this->NumberY << ", Z = " << this->NumberZ << '\n';
      File << "Number of cells under and above the dot: Under = " << this->under << ", Above = " << this->above << '\n';
      File << "Probability given: " << proportion1 << ", " << proportion2 << '\n';
      File << "Real proportion: " << percent << '\n';      
      File << "Base and top radius: Base = " << Rb << ", Top = " << Rt << '\n';
      File << "Wetting layer thickness: " << w << '\n';
      File << "Electric fields under the dot, in the wetting layers, in and above the dot: " << down_field << ", " << wetting_field << ", " << dot_field << ", " << up_field << '\n';
      File << "Band offset: " << offset << '\n';
      File << "Growth direction lattice constant: " << c << endl;
      File.close();
    }

  cout << "Real proportion: " << percent << endl;
}


bool ThreeDPotential::Neighborhood(int m, int n, int h)
{
  int m1 = m - 1;
  int m2 = m + 1;
  int n1 = n - 1;
  int n2 = n + 1;
  int h1 = h - 1;
  int h2 = h + 1;
  if (m1 < 0) m1 = 0;
  if (m2 >= this->NumberX) m2 = this->NumberX - 1;
  if (n1 < 0) n1 = 0;
  if (n2 >= this->NumberY) n2 = this->NumberY - 1;
  if (h1 < 0) h1 = 0;
  if (h2 >= this->NumberZ) h2 = this->NumberZ - 1;

  for (int i = m1; i <= m2; ++i)
    for (int j = n1; j <= n2; ++j)
      for (int k = h1; k <= h2; ++k)
	if (this->Alloy[k][j][i])
	  return true;

  return false;
}

bool ThreeDPotential::NeighborhoodBis(int m, int n, int h)
{
  int m1 = m - 1;
  int m2 = m + 1;
  int n1 = n - 1;
  int n2 = n + 1;
  int h1 = h - 1;
  int h2 = h + 1;
  if (m1 < 0) m1 = 0;
  if (m2 >= this->NumberX) m2 = this->NumberX - 1;
  if (n1 < 0) n1 = 0;
  if (n2 >= this->NumberY) n2 = this->NumberY - 1;
  if (h1 < 0) h1 = 0;
  if (h2 >= this->NumberZ) h2 = this->NumberZ - 1;

  for (int i = m1; i <= m2; ++i)
    if (this->Alloy[h][n][i])
      return true;
  for (int j = n1; j <= n2; ++j)
    if (this->Alloy[h][j][m])
      return true;
  for (int k = h1; k <= h2; ++k)
    if (this->Alloy[k][n][m])
      return true;

  return false;
}
