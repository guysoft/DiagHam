#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFileOption ('\0', "lz-file", "name of the file containing the lz sorted spectrum", 0);
  SingleDoubleOption PrecisionOption ('\n', "precision", "precision used to compare two energy values", 1e-7);
  BooleanOption FullOption ('f', "full", "indicates that all Lz value are contained in the file (if not L output will have one less value that the Lz counterpart)");
  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &InputFileOption;
  OptionList += &PrecisionOption;
  OptionList += &FullOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type LzToL -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }

  char* InputName = InputFileOption.GetString();
  double Precision = PrecisionOption.GetDouble();
  ifstream InputFile;
  InputFile.open(InputName, ios::binary | ios::in);
  InputFile.seekg(0, ios::end);
  int FileSize = InputFile.tellg();
  --FileSize;
  InputFile.seekg(0, ios::beg);
  int NbrValue = 1;
  int TotalSize = 0;
  int CurrentLzValue;
  int TmpLzValue;
  double Dummy;
  InputFile >> CurrentLzValue >> Dummy;
  while ((InputFile.tellg() < FileSize) && (InputFile.tellg() >= 0))
    {
      InputFile >> TmpLzValue >> Dummy;
      if (TmpLzValue != CurrentLzValue)
        {
          CurrentLzValue = TmpLzValue;
          ++NbrValue;
        }
      ++TotalSize;
    }
  InputFile.close();

  int* Dimensions = new int [NbrValue];
  double** Eigenvalues = new double* [NbrValue];
  int Pos = 0;
  ifstream InputFile2;
  InputFile2.open(InputName, ios::binary | ios::in);
  InputFile2 >> CurrentLzValue >> Dummy;
  Dimensions[Pos] = 1;  
  while ((InputFile2.tellg() < FileSize) && (InputFile2.tellg() >= 0))
    {
      InputFile2 >> TmpLzValue >> Dummy;
      if (TmpLzValue != CurrentLzValue)
        {
          CurrentLzValue = TmpLzValue;
          Eigenvalues[Pos] = new double [Dimensions[Pos]];
          ++Pos;
          Dimensions[Pos] = 0;
        }
      ++Dimensions[Pos];
    }
  InputFile2.close();
  Eigenvalues[Pos] = new double [Dimensions[Pos]];

  Pos = 0;
  int Pos2 = 0;
  ifstream InputFile3;
  InputFile3.open(InputName, ios::binary | ios::in);
  InputFile3 >> CurrentLzValue >> Dummy;
  Eigenvalues[Pos][Pos2] = Dummy;
  ++Pos2;
  while ((InputFile3.tellg() < FileSize) && (InputFile3.tellg() >= 0))
    {
      InputFile3 >> TmpLzValue >> Dummy;
      if (TmpLzValue != CurrentLzValue)
        {
          ++Pos;
          CurrentLzValue = TmpLzValue;
          Pos2 = 0;
        }
      Eigenvalues[Pos][Pos2] = Dummy;
      ++Pos2;
    }
  InputFile3.close();


  int SpectrumSize = 0;
  double* Spectrum = new double [TotalSize];
  bool* Degeneracy = new bool [TotalSize];
  bool Flag;
  double TmpEigenvalue;
  if (FullOption.GetBoolean() == false)
    for (int i = 0; i < Dimensions[NbrValue - 1]; ++i)
      {
	Spectrum[i] = Eigenvalues[NbrValue - 1][i];
	++SpectrumSize;
      }
  else
    {
      for (int i = 0; i < Dimensions[NbrValue - 1]; ++i)
	{
	  Spectrum[i] = Eigenvalues[NbrValue - 1][i];
	  cout << (NbrValue - 1) << " " << Eigenvalues[NbrValue - 1][i] << endl;
	  ++SpectrumSize;
	}
    }
  for (int  L = NbrValue - 2; L >= 0; --L)
    {
      for (int j = 0; (j < SpectrumSize); ++j)
        Degeneracy[j] = false;
      for (int i = 0; i < Dimensions[L]; ++i)
        {
          Flag = false;
          TmpEigenvalue = Eigenvalues[L][i];
          for (int j = 0; ((j < SpectrumSize) && (Flag == false)); ++j)
            if ((Degeneracy[j] == false) && (fabs((TmpEigenvalue - Spectrum[j]) / TmpEigenvalue) < Precision))
              {
                Flag = true;
                Degeneracy[j] = true;
              }
          if (Flag == false)
            {
              Spectrum[SpectrumSize] = TmpEigenvalue;
	      cout << L << " " << TmpEigenvalue << endl;
              Degeneracy[SpectrumSize] = false;
              ++SpectrumSize;
            }
        }
       delete[] Eigenvalues[L];
    }
  delete[] Degeneracy;
  delete[] Spectrum;
  delete[] Eigenvalues;
  delete[] Dimensions;

  return 0;
}

