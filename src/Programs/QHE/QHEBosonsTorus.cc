#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Hamiltonian/ParticleOnTorusCoulombHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusDeltaHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "QuantumNumber/AbstractQuantumNumber.h"

#include "GeneralTools/ListIterator.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 40);
  BooleanOption GroundOption ('g', "ground", "restrict to the largest subspace");
  SingleIntegerOption NbrBosonOption ('p', "nbr-particles", "number of particles", 5);
  SingleIntegerOption MaxMomentumOption ('l', "max-momentum", "maximum momentum for a single particle", 15);
  SingleIntegerOption MomentumOption ('m', "momentum", "constraint on the total momentum modulo the maximum momentum (negative if none)", -1);
  SingleDoubleOption RatioOption ('r', "ratio", "ratio between the two torus lengths", 0.57735026919);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &GroundOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &NbrBosonOption;
  OptionList += &MaxMomentumOption;
  OptionList += &MomentumOption;
  OptionList += &RatioOption;

  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsTorus -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }


  bool GroundFlag = GroundOption.GetBoolean();
  bool SMPFlag = SMPOption.GetBoolean();
  int NbrProcessor = SMPNbrProcessorOption.GetInteger();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  int NbrBosons = NbrBosonOption.GetInteger();
  int MaxMomentum = MaxMomentumOption.GetInteger();
  int Momentum = MomentumOption.GetInteger();
  double XRatio = RatioOption.GetDouble();

  int L = 0;
  double GroundStateEnergy = 0.0;

  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "bosons_torus_delta_n_%d_2s_%d_ratio_%f.dat", NbrBosons, MaxMomentum, XRatio);
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);

  int Max = (MaxMomentum - 1);
  if (Momentum < 0)
    Momentum = 0;
  else
    Max = Momentum;

  for (; Momentum <= Max; ++Momentum)
    {     
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
      BosonOnTorus TotalSpace (NbrBosons, MaxMomentum, Momentum);
      cout << " Hilbert space dimension = " << TotalSpace.GetHilbertSpaceDimension() << endl;
      for (int i = 0; i < TotalSpace.GetHilbertSpaceDimension(); ++i)
	{
	  cout << i << " = ";
	  TotalSpace.PrintState(cout, i) << endl;
	}
/*	List<AbstractQuantumNumber*> QuantumNumbers ( TotalSpace.GetQuantumNumbers());
      ListIterator<AbstractQuantumNumber*> QuantumNumberIter (QuantumNumbers);
      AbstractQuantumNumber** TmpQuantumNumber;
      while ((TmpQuantumNumber = QuantumNumberIter()))
	{
	  cout << "momentum = " << (**TmpQuantumNumber) << endl;
	  SubspaceSpaceConverter Converter;
	  BosonOnTorus* Space = (BosonOnTorus*) TotalSpace.ExtractSubspace (**TmpQuantumNumber, Converter);
	  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	    {
	      cout << i << " = ";
	      Space->PrintState(cout, i) << endl;
	    }
	  cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;
	}*/

      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
//      AbstractHamiltonian* Hamiltonian = new ParticleOnTorusCoulombHamiltonian (&TotalSpace, NbrBosons, MaxMomentum, XRatio);
      AbstractHamiltonian* Hamiltonian = new ParticleOnTorusDeltaHamiltonian (&TotalSpace, NbrBosons, MaxMomentum, XRatio, Architecture);

      if (Hamiltonian->GetHilbertSpaceDimension() < 300)
	{
	  RealSymmetricMatrix HRep (Hamiltonian->GetHilbertSpaceDimension());
	  Hamiltonian->GetHamiltonian(HRep);
	  if (Hamiltonian->GetHilbertSpaceDimension() > 1)
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian->GetHilbertSpaceDimension());
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      if (L == 0)
		GroundStateEnergy = TmpTriDiag.DiagonalElement(0);
	      for (int j = 0; j < Hamiltonian->GetHilbertSpaceDimension() ; j++)
		{
		  File << Momentum << " " << TmpTriDiag.DiagonalElement(j) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      File << Momentum << " " << HRep(0, 0) << endl;
	    }
	}
      else
	{
	  AbstractLanczosAlgorithm* Lanczos = 0;
	  if (NbrEigenvalue > 1)
	    Lanczos = new FullReorthogonalizedLanczosAlgorithm (Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	  else
	    Lanczos = new BasicLanczosAlgorithm (Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	    
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = NbrEigenvalue + 3;
	  Lanczos->SetHamiltonian(Hamiltonian);
	  Lanczos->InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Precision > 1e-14) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
	    {
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue);
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
	    }
	  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
	    {
	      cout << "too much Lanczos iterations" << endl;
	      File << "too much Lanczos iterations" << endl;
	      File.close();
	      exit(0);
	    }
	  GroundStateEnergy = Lowest;
	  cout << endl;
	  cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	       << CurrentNbrIterLanczos << endl;
	  for (int i = 0; i <= NbrEigenvalue; ++i)
	    {
	      cout << TmpMatrix.DiagonalElement(i) << " ";
	      File << Momentum << " " << TmpMatrix.DiagonalElement(i) << endl;
	    }
	  cout << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	  delete Lanczos;
	}
      cout << "----------------------------------------------------------------" << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrBosons) << endl;
      delete Hamiltonian;
    }
  File.close();

  return 0;
}
