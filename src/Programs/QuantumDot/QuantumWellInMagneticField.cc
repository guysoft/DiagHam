#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/HermitianMatrix.h"

#include "Hamiltonian/QuantumDotHamiltonian/QuantumWellHamiltonianInMagneticField.h"


#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ostream;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  BooleanOption LanczosOption ('l', "lanczos", "enable lanczos diagonalization algorithm", true);
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption VerboseOption ('v', "verbose", "verbose mode", true);
  BooleanOption EigenstateOption ('e', "eigenstate", "evaluate eigenstates", true);
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 6);
  SingleIntegerOption MemoryOption ('\n', "memory", "amount of memory that can be used for precaching (in Mb)", 1000);
  SingleIntegerOption MValueOption ('M', "M-cell", "number of cells in the x direction", 50);
  SingleIntegerOption NValueOption ('N', "N-cell", "number of cells in the y direction", 50);
  SingleIntegerOption HValueOption ('H', "H-cell", "number of cells in the z direction", 30);
  SingleIntegerOption LeftSizeOption ('\n', "left-size", "size of the leftmost part in the z direction with constant null potential (in cell unit)", 6);
  SingleIntegerOption RightSizeOption ('\n', "right-size", "size of the rightmost part in the z direction with constant potential (in cell unit)", 10);
  SingleDoubleOption CellXSizeOption ('X', "cell-xsize", "cell size in the x direction in Angstrom", 2.97);
  SingleDoubleOption CellYSizeOption ('Y', "cell-ysize", "cell size in the y direction in Angstrom", 2.97);
  SingleDoubleOption CellZSizeOption ('Z', "cell-zsize", "cell size in the z direction in Angstrom", 2.64);  
  SingleDoubleOption MassOption ('\n', "mass", "electron effective mass (in bare electron mass unit)", 0.050);
  SingleDoubleOption BFieldOption ('\n', "bfield", "B field (in Tesla unit)", 31.09);
  SingleStringOption CoefficientFileNameOption('\n', "coefficients", "name of the file where interaction coeffcients are stored", 
					       "/home/regnault/development/DMRG/DiagHam/potentiel_10_10_10_2");
  BooleanOption CarrierTypeOption('c', "carrier", "carrier type, true for hole, false for electron", true);

  List<AbstractOption*> OptionList;
  OptionList += &LanczosOption;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &VerboseOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &EigenstateOption;
  OptionList += &MValueOption;
  OptionList += &NValueOption;
  OptionList += &HValueOption;
  OptionList += &CellXSizeOption;
  OptionList += &CellYSizeOption;
  OptionList += &CellZSizeOption;
  OptionList += &MassOption;
  OptionList += &BFieldOption;
  OptionList += &CoefficientFileNameOption;
  OptionList += &LeftSizeOption;
  OptionList += &RightSizeOption;
  OptionList += &MemoryOption;
  OptionList += &CarrierTypeOption; 

  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type ExplicitMatrixExample -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }

  int Memory = MemoryOption.GetInteger();
  bool LanczosFlag = LanczosOption.GetBoolean();
  bool SMPFlag = SMPOption.GetBoolean();
  bool VerboseFlag = VerboseOption.GetBoolean();
  bool EigenstateFlag = EigenstateOption.GetBoolean();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  char* CoefficientFileName = CoefficientFileNameOption.GetString();
  double Lx = CellXSizeOption.GetDouble();
  double Ly = CellYSizeOption.GetDouble();
  double Lz = CellZSizeOption.GetDouble();
  double Mass = MassOption.GetDouble();
  double BField = BFieldOption.GetDouble();
  int LeftSize = LeftSizeOption.GetInteger();
  int RightSize = RightSizeOption.GetInteger();
  bool Carrier = CarrierTypeOption.GetBoolean();

  QuantumWellHamiltonianInMagneticField Hamiltonian (1800.0, 1800.0, 58.7, Mass, BField, 0.0, 143.0, 2, 0, 5.87, 600, 0.53);
  cout << Hamiltonian.GetHilbertSpaceDimension() << endl;
  HermitianMatrix HamiltonianRepresentation;
  Hamiltonian.GetHamiltonian(HamiltonianRepresentation);
  double TmpElement;
  for (int i = 0; i < HamiltonianRepresentation.GetNbrRow(); ++i)
    {
      HamiltonianRepresentation.GetMatrixElement(i, i, TmpElement);
      cout << TmpElement << " ";
    }
  cout << endl;
  return 0;
}

