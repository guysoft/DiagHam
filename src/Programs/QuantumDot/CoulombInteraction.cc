#include "Vector/RealVector.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;

double FormFunction(RealVector& Electron, RealVector& Hole, double WaveVector, double Step);
void PlanarIntegral();

int main(int argc, char** argv)
{
  ifstream Electron("/home/phuong/Helios/DF/12E.txt");
  if (!Electron.is_open())
    {
      cout << "Error when open the electron state file. Exit now" << endl;
      exit(1);
    }
  ifstream Hole("/home/phuong/Helios/DF/12H.txt");
  if (!Hole.is_open())
    {
      cout << "Error when open the hole state file. Exit now" << endl;
      exit(1);
    }
  int N = 1268; double tmp = 0.0; 
  RealVector electron(1268); RealVector hole(1268); 
  for (int i = 0; i < N; ++i)
    {
      Electron >> tmp; 
      Electron >> tmp; 
      electron[i] = tmp * tmp;
      Hole >> tmp;
      Hole >> tmp;
      hole[i] = tmp * tmp;
    }

  PlanarIntegral();

  /*
  ofstream FormFile("FormZ.txt");
  if (!FormFile.is_open())
    {
      cout << "Error when open the form file to write. Exit now" << endl;
      exit(1);
    }    
  double Q = 0.0;
  for (; Q < 3.0; Q += 0.1)
    {
      FormFile << Q << '\t';
      FormFile << FormFunction(electron, hole, Q, 0.2) << '\n';
    }
  */

  Electron.close(); Hole.close();
  return 0;
}

// Electron: probability density
// Hole: probability density
// Step: step in z direction
double FormFunction(RealVector& Electron, RealVector& Hole, double WaveVector, double Step)
{
  double Tmp = 0.0; double TmpTotal = 0.0;
  int N = Electron.GetVectorDimension();
  WaveVector = -WaveVector * Step;

  for (int z1 = 0; z1 < N; ++z1)
    {
      Tmp = 0.0;
      for (int z2 = 0; z2 < z1; ++z2)	
	Tmp += (Hole[z2] * exp(WaveVector * (z1 - z2)));
      for (int z2 = z1; z2 < N ; ++z2)
	Tmp += (Hole[z2] * exp(WaveVector * (z2 - z1)));	
      TmpTotal += (Electron[z1] * Tmp );
    }
  return Step * Step * TmpTotal;
}

void PlanarIntegral()
{
  int NbrStateE = 40; int NbrStateH = 40;
  int NbrStateX = 120; int NbrStateY = 120;
  int DemiNbrQ = 60; int NbrQ = 1 + DemiNbrQ * 2;

  double*** InterM = new double** [NbrQ];
  double*** InterN = new double** [NbrQ];
  double*** DeltaM = new double** [NbrQ];
  double*** DeltaN = new double** [NbrQ];
  for (int q = 0; q < NbrQ; ++q)
    {
      InterM[q] = new double* [NbrStateX];
      DeltaM[q] = new double* [NbrStateX];
      InterN[q] = new double* [NbrStateY];
      DeltaN[q] = new double* [NbrStateY];
      for (int m1 = 0; m1 < NbrStateX; ++m1)
	{
	  InterM[q][m1] = new double [NbrStateX];
	  DeltaM[q][m1] = new double [NbrStateX];
	  for (int m2 = 0; m2 < NbrStateX; ++m2)
	    {
	      InterM[q][m1][m2] = 0.0; DeltaM[q][m1][m2] = 0.0;
	    }	  
	}
      for (int n1 = 0; n1 < NbrStateY; ++n1)
	{
	  InterN[q][n1] = new double [NbrStateY];
	  DeltaN[q][n1] = new double [NbrStateY];
	  for (int n2 = 0; n2 < NbrStateY; ++n2)
	    {
	      InterN[q][n1][n2] = 0.0; DeltaN[q][n1][n2] = 0.0;
	    }	  
	}
    }
  // q = 0
  for (int m = 0; m < NbrStateX; ++m)    
    DeltaM[0][m][m] = 0.25;    
  for (int n = 0; n < NbrStateY; ++n)
    DeltaM[0][n][n] = 0.25;

  // positive q
  int tmpQ = 0; int tmpIndex1 = 0, tmpIndex2 = 0; double tmpInter = 0.0;
  for (int q = 1; q <= DemiNbrQ; ++q)
    {
      tmpQ = 2 * q;
      for (int m1 = 0; m1 < NbrStateX; ++m1)
	{ 
	  tmpIndex1 = m1 + tmpQ; tmpIndex2 = m1 - tmpQ; 
	  for (int m2 = 0; m2 < NbrStateX; ++m2)
	    {
	      if ((m2 == tmpIndex1) || (m2 == -tmpIndex1) || (m2 == tmpIndex2) || (m2 == -tmpIndex2))
		DeltaM[q][m1][m2] = 0.25;
	      
	    }
	}
      for (int n1 = 0; n1 < NbrStateY; ++n1)
	{ 
	  tmpIndex1 = n1 + tmpQ; tmpIndex2 = n1 - tmpQ; 
	  for (int n2 = 0; n2 < NbrStateY; ++n2)
	    if ((n2 == tmpIndex1) || (n2 == -tmpIndex1) || (n2 == tmpIndex2) || (n2 == -tmpIndex2))
	      DeltaN[q][n1][n2] = 0.25;
	}
      for (int m1 = 1; m1 < NbrStateX; m1 += 2)
	{
	  tmpIndex1 = tmpQ + m1; tmpIndex2 = tmpQ - m1; tmpInter = 8 * q * m1 / M_PI;
	  for (int m2 = 2; m2 < NbrStateX; m2 += 2)
	    {
	      InterM[q][m1][m2] = tmpInter * m2 / ((tmpIndex1 + m2) * (tmpIndex1 - m2) * (tmpIndex2 + m2) * (tmpIndex2 - m2));
	      InterM[q][m2][m1] = InterM[q][m1][m2];
	    }
	}
      for (int n1 = 1; n1 < NbrStateY; n1 += 2)
	{
	  tmpIndex1 = tmpQ + n1; tmpIndex2 = tmpQ - n1; tmpInter = 8 * q * n1 / M_PI;
	  for (int n2 = 2; n2 < NbrStateY; n2 += 2)
	    {
	      InterM[q][n1][n2] = tmpInter * n2 / ((tmpIndex1 + n2) * (tmpIndex1 - n2) * (tmpIndex2 + n2) * (tmpIndex2 - n2));
	      InterM[q][n2][n1] = InterM[q][n1][n2];
	    }
	}
    }

  // negative q
  for (int q = 1 + DemiNbrQ; q < NbrQ; ++q)
    {
      tmpQ = 2 * (DemiNbrQ - q);
      for (int m1 = 0; m1 < NbrStateX; ++m1)
	{ 
	  tmpIndex1 = m1 + tmpQ; tmpIndex2 = m1 - tmpQ; 
	  for (int m2 = 0; m2 < NbrStateX; ++m2)
	    if ((m2 == tmpIndex1) || (m2 == -tmpIndex1) || (m2 == tmpIndex2) || (m2 == -tmpIndex2))
	      DeltaM[q][m1][m2] = 0.25;
	}
      for (int n1 = 0; n1 < NbrStateY; ++n1)
	{ 
	  tmpIndex1 = n1 + tmpQ; tmpIndex2 = n1 - tmpQ; 
	  for (int n2 = 0; n2 < NbrStateY; ++n2)
	    if ((n2 == tmpIndex1) || (n2 == -tmpIndex1) || (n2 == tmpIndex2) || (n2 == -tmpIndex2))
	      DeltaN[q][n1][n2] = 0.25;
	}
      for (int m1 = 1; m1 < NbrStateX; m1 += 2)
	{
	  tmpIndex1 = tmpQ + m1; tmpIndex2 = tmpQ - m1; tmpInter = 8 * q * m1 / M_PI;
	  for (int m2 = 2; m2 < NbrStateX; m2 += 2)
	    {
	      InterM[q][m1][m2] = tmpInter * m2 / ((tmpIndex1 + m2) * (tmpIndex1 - m2) * (tmpIndex2 + m2) * (tmpIndex2 - m2));
	      InterM[q][m2][m1] = InterM[q][m1][m2];
	    }
	}
      for (int n1 = 1; n1 < NbrStateY; n1 += 2)
	{
	  tmpIndex1 = tmpQ + n1; tmpIndex2 = tmpQ - n1; tmpInter = 8 * q * n1 / M_PI;
	  for (int n2 = 2; n2 < NbrStateY; n2 += 2)
	    {
	      InterM[q][n1][n2] = tmpInter * n2 / ((tmpIndex1 + n2) * (tmpIndex1 - n2) * (tmpIndex2 + n2) * (tmpIndex2 - n2));
	      InterM[q][n2][n1] = InterM[q][n1][n2];
	    }
	}
    }

  double*** FunctionE = new double** [NbrStateE];
  double*** FunctionH = new double** [NbrStateH];
  char* FileNameE = new char [100];
  ifstream FileE(FileNameE);
  if (!FileE.is_open())
    {
      cout << "Error when open the planar electron state file. Exit now" << endl;
      exit(1);
    }
  for (int e = 0; e < NbrStateE; ++e)
    {
      FunctionE[e] = new double* [NbrStateY];    
      for (int n = 0; n < NbrStateY; ++n)
	{	  
	  FunctionE[e][n] = new double [NbrStateX];
	  for (int m = 0; m < NbrStateX; ++m)	 
	    FileE >> FunctionE[e][n][m];	    
	}
    }
  char* FileNameH = new char [100];
  ifstream FileH(FileNameH);
  if (!FileH.is_open())
    {
      cout << "Error when open the planar hole state file. Exit now" << endl;
      exit(1);
    }
  for (int h = 0; h < NbrStateH; ++h)
    {
      FunctionH[h] = new double* [NbrStateY];    
      for (int n = 0; n < NbrStateY; ++n)
	{	  
	  FunctionE[h][n] = new double [NbrStateX];
	  for (int m = 0; m < NbrStateX; ++m)
	    FileH >> FunctionH[h][n][m];	    
	}
    }
  double**** CosE = new double***[NbrStateE];
  double**** SinE = new double***[NbrStateE];
  double**** CosH = new double***[NbrStateH];
  double**** SinH = new double***[NbrStateH];

  double TmpSin = 0.0; double TmpCos = 0.0; double TmpTotalSin = 0.0; double TmpTotalCos = 0.0;
  for (int e1 = 0; e1 < NbrStateE; ++e1)
    {
      CosE[e1] = new double** [e1 + 1];
      SinE[e1] = new double** [e1 + 1];
      for (int e2 = 0; e2 <= e1; ++e2)
	{
	  CosE[e1][e2] = new double* [NbrQ];
	  SinE[e1][e2] = new double* [NbrQ];
	  for (int q1 = 0; q1 < NbrQ; ++q1)
	    {
	      CosE[e1][e2][q1] = new double[DemiNbrQ + 1];
	      SinE[e1][e2][q1] = new double[DemiNbrQ + 1];
	      for (int q2 = 0; q2 < (DemiNbrQ + 1); ++q2)
		{
		  TmpTotalSin = 0.0; TmpTotalCos = 0.0;
		  for (int n1 = 0; n1 < NbrStateY; ++n1)
		    for (int m1 = 0; m1 < NbrStateX; ++m1)
		      {
			TmpSin = 0.0; TmpCos = 0.0;
			for (int n2 = 0; n2 < NbrStateY; ++n2)
			  for (int m2 = 0; m2 < NbrStateX; ++m2)			    
			    {
			      TmpCos += FunctionE[e2][n2][m2] * (DeltaM[q1][m1][m2] * DeltaN[q2][n1][n2] - InterM[q1][m1][m2] * InterN[q2][n1][n2]);
			      TmpSin += FunctionE[e2][n2][m2] * (DeltaM[q1][m1][m2] * InterN[q2][n1][n2] + DeltaN[q2][n1][n2] * InterM[q1][m1][m2]);
			    }
			TmpTotalCos += FunctionE[e1][n1][m1] * TmpCos;
			TmpTotalSin += FunctionE[e1][n1][m1] * TmpSin;			
		      }
		  CosE[e1][e2][q1][q2] = TmpTotalCos * 4.0;
		  SinE[e1][e2][q1][q2] = TmpTotalSin * 4.0;
		}
	    }
	}
    }
  for (int h1 = 0; h1 < NbrStateH; ++h1)
    {
      CosH[h1] = new double** [h1 + 1];
      SinH[h1] = new double** [h1 + 1];
      for (int h2 = 0; h2 <= h1; ++h2)
	{
	  CosH[h1][h2] = new double* [NbrQ];
	  SinH[h1][h2] = new double* [NbrQ];
	  for (int q1 = 0; q1 < NbrQ; ++q1)
	    {
	      CosH[h1][h2][q1] = new double[DemiNbrQ + 1];
	      SinH[h1][h2][q1] = new double[DemiNbrQ + 1];
	      for (int q2 = 0; q2 < (DemiNbrQ + 1); ++q2)
		{
		  TmpTotalSin = 0.0; TmpTotalCos = 0.0;
		  for (int n1 = 0; n1 < NbrStateY; ++n1)
		    for (int m1 = 0; m1 < NbrStateX; ++m1)
		      {
			TmpSin = 0.0; TmpCos = 0.0;
			for (int n2 = 0; n2 < NbrStateY; ++n2)
			  for (int m2 = 0; m2 < NbrStateX; ++m2)			    
			    {
			      TmpCos += FunctionH[h2][n2][m2] * (DeltaM[q1][m1][m2] * DeltaN[q2][n1][n2] - InterM[q1][m1][m2] * InterN[q2][n1][n2]);
			      TmpSin += FunctionH[h2][n2][m2] * (DeltaM[q1][m1][m2] * InterN[q2][n1][n2] + DeltaN[q2][n1][n2] * InterM[q1][m1][m2]);
			    }
			TmpTotalCos += FunctionH[h1][n1][m1] * TmpCos;
			TmpTotalSin += FunctionH[h1][n1][m1] * TmpSin;			
		      }
		  CosH[h1][h2][q1][q2] = TmpTotalCos * 4.0;
		  SinH[h1][h2][q1][q2] = TmpTotalSin * 4.0;
		}
	    }
	}
    }
}
