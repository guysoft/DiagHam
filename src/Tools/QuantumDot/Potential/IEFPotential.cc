#include "Potential/ThreeDPotential.h"

#include <time.h>
#include <stdlib.h>
#include <math.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ofstream;

void ThreeDPotential::IEFFluctuateWelll(int height, int width, double down_field, double field, double up_field, double cell, double offset, char* fileName = 0)
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
      tmp = (k - this->under) * down_field - (this->NumberZ - this->under - this->above) * field;
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
      tmp = (k + this->above - this->NumberZ) * field - offset;
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
      tmp = (k - this->NumberZ + this->above) * up_field;
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
      tmp = (k - this->NumberZ + this->above) * up_field;
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
