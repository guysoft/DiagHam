#include "Tools/QuantumDot/Spectra/DOSSpectra.h"
#include "Vector/RealVector.h"

#include <fstream>
#include <math.h>
#include <stdlib.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;

// constructor from a set of energy files. Each peak is assimilated to a Lorentzian function.
//
// FileNumber: number of files, Files: name of files
// StateNumber: integer array containing number of states in each file
// Gamma: FWHM
// Emin, Emax, dE: two energy bounds and the step in energy
DOSSpectra::DOSSpectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE)
{
  int N = (Emax - Emin) / dE;
  double * Energy = new double [N];
  double * DOS = new double [N]; double tmp1; double tmp2 = 0.0; double g = Gamma * Gamma / 4;

  for (int i = 0; i < N; ++i)
    {
      Energy[i] = Emin + dE * i;
      DOS[i] = 0.0;
    }

  for (int i = 0; i < FileNumber; ++i)
    {
      int n = StateNumber[i];
      double* tmp = new double [n];
      ifstream file;
      file.open(Files[i],ios::out);
      if (!file.is_open())
        {
	  cout << "Error in open the file: " << Files[i] << "Exit now" << endl;
	  exit(0);
	}
      for (int j = 0; j < n; ++j)
	file >> tmp[j];
      file.close();

      for (int j = 0; j < N; ++j)
	{
	  tmp1 = Energy[j]; tmp2 = 0.0;
	  for (int k = 0; k < n; ++k)
	    tmp2 += 1 / ((tmp1 - tmp[k]) * (tmp1 - tmp[k]) + g);
	  DOS[j] += tmp2;
	}
      delete[] tmp;
      tmp = 0;
    }
  double tmp3 = Gamma/(2 * M_PI * FileNumber);
  for (int i = 0; i < N; ++i)
    DOS[i] *= tmp3;

  this->AxeX = new RealVector(Energy, N);
  this->AxeY = new RealVector(DOS, N);
  this->PointNumber = N;
}
