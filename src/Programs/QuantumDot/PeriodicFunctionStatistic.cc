#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/QuantumDot/Spectra/Spectra.h"
#include "Tools/QuantumDot/Spectra/PeriodicSpectra.h"

#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;

int compare_doubles (const void * a, const void * b);

int main(int argc, char** argv)
{
 
  // some running options and help 
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFile('\n', "input", "name of the input file", 0);
  SingleIntegerOption XCell('X', "Xcells", "number of the cells in X direction", 61);
  SingleIntegerOption YCell('Y', "Ycells", "number of the cells in Y direction", 61);
  SingleIntegerOption ZCell('Z', "Zcells", "number of the cells in Z direction", 60);  
  SingleIntegerOption ScaleOption('S', "scale", "number of logarithmic scake", 14);
  SingleStringOption Output('r', "output", "name of the output file", "");

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &InputFile;
  OptionList += &XCell;
  OptionList += &YCell;
  OptionList += &ZCell;
  OptionList += &ScaleOption;
  OptionList += &Output;

  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "bad options" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }

  int M, N, H;
  char* FileName;
  int Scale;

  FileName = InputFile.GetString();
  M = XCell.GetInteger();
  N = YCell.GetInteger();
  H = ZCell.GetInteger();
  Scale = ScaleOption.GetInteger();
  char * out = Output.GetString();


  Periodic3DOneParticle* Space = new Periodic3DOneParticle(M, M / 2, N, N / 2, H / 2 + 1, H / 4);

  int NbrStateX = Space->GetNbrStateX();
  int NbrStateY = Space->GetNbrStateY();
  int NbrStateZ = Space->GetNbrStateZ();
  int Dimension = NbrStateX * NbrStateY * NbrStateZ; 

  ifstream File(FileName);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << FileName << endl;
      exit(1);
    }

  File.precision(14);

  double* Probability = new double [Dimension];
  double Re, Im;
  int Index = 0;
  double Norm = 0.0;
  for (int i = 0; i < NbrStateX; ++i)
    for (int j = 0; j < NbrStateY; ++j)
      for (int k = 0; k < NbrStateZ; ++k)
	{
	  File >> Re >> Im;
	  Probability[Index] = Re * Re + Im * Im;
	  Norm += Probability[Index];
	  ++Index;
	}  
  
  qsort(Probability, Dimension, sizeof(double), compare_doubles);

  double max = Probability[Dimension - 1];
  cout << "Maximal probability: " << max << endl;
  cout << "Norm: " << Norm << endl;

  double* LogScale = new double [Scale];

  LogScale[Scale - 1] = max;
  for (int i = Scale - 2; i > 0; --i)
    LogScale[i] = LogScale[i + 1] / 10.0;
  
  int tNbr = 0; Index = 0;
  for (int i = 0; i < Scale; ++i)
    {
      tNbr = 0;      
      while (Probability[Index] < LogScale[i])
	{
	  ++tNbr;
	  ++Index;
	}
      cout << i - Scale + 1 << '\t' << tNbr << endl;
    }
  File.close();

  return 0;
}

int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}
