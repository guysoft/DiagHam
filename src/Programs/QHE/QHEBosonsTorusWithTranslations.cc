#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/QHEHilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/QHEHilbertSpace/BosonOnTorusState.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnTorusCoulombHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 40);
  BooleanOption GroundOption ('g', "ground", "restrict to the largest subspace");
  SingleIntegerOption NbrFermionOption ('p', "nbr-particles", "number of particles", 6);
  SingleIntegerOption MaxMomentumOption ('l', "max-momentum", "maximum momentum for a single particle", 9);
  SingleIntegerOption MomentumOption ('m', "momentum", "constraint on the total momentum modulo the maximum momentum (negative if none)", -1);
  SingleIntegerOption MaxFullDiagonalizationOption ('f', "max-full", "maximum hilbert space size allowed to use full diagonalization", 300);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &GroundOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &NbrFermionOption;
  OptionList += &MaxMomentumOption;
  OptionList += &MomentumOption;
  OptionList += &MaxFullDiagonalizationOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsTorus -h" << endl;
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
  int NbrFermions = NbrFermionOption.GetInteger();
  int MaxMomentum = MaxMomentumOption.GetInteger();
  int Momentum = MomentumOption.GetInteger();
  int MaxFullDiagonalization = MaxFullDiagonalizationOption.GetInteger();
  double XRatio = NbrFermions / 4.0;

  int InvNu = 3;
//  int MaxMomentum = InvNu * NbrFermions;
  int L = 0;
  double GroundStateEnergy = 0.0;

  int NbrState = 9;
  int ReducedNbrState = NbrState >> 2;
  int NbrStateRemainder = NbrState - (ReducedNbrState << 2);
  if (NbrStateRemainder == 0)
    {
      NbrStateRemainder = 4;
      --ReducedNbrState;
    }
  for (int k = 0; k < NbrState; ++k)
    {
      BosonOnTorusState State (ReducedNbrState + 1);
      for (int i = 0; i < NbrState; ++i)
	State.SetOccupation(i, 65 | (i << 1));
      BosonOnTorusState TmpState(State, ReducedNbrState + 1);
      State.LeftShiftState(ReducedNbrState, NbrStateRemainder, k);
      for (int i = 0; i < NbrState; ++i)
	if (((i >=  k) && (State.GetOccupation(i - k) != TmpState.GetOccupation(i)))
	    || ((i <  k) && (State.GetOccupation(NbrState + i -  k) != TmpState.GetOccupation(i))))
	  cout << "error " << i << endl;
      State.PrintState(cout, ReducedNbrState, NbrStateRemainder) << endl;
    }
  return 0;

  char* OutputNameLz = new char [512];
  sprintf (OutputNameLz, "fermions_torus_coulomb_n_%d_2s_%d_ratio_%f.dat", NbrFermions, MaxMomentum, XRatio);
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);


  
  AbstractArchitecture* Architecture = 0;
  if (SMPFlag == false)
    Architecture = new MonoProcessorArchitecture;
  else
    Architecture = new SMPArchitecture(NbrProcessor);

  int Max = (MaxMomentum - 1);
  if (Momentum < 0)
    Momentum = 0;
  else
    Max = Momentum;

  int MomentumModulo = FindGCD(NbrFermions, MaxMomentum);
  MomentumModulo = 1;
//  for (; Momentum <= Max; ++Momentum)
  for (int x = 0; x < MomentumModulo; ++x)
  for (int y = 0; y < MomentumModulo; ++y)
    {     
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
//      FermionOnTorus TotalSpace (NbrFermions, MaxMomentum, Momentum);
      FermionOnTorusWithMagneticTranslations TotalSpace (NbrFermions, MaxMomentum, x, y);
      cout << " Total Hilbert space dimension = " << TotalSpace.GetHilbertSpaceDimension() << endl;
//      cout << "momentum = " << Momentum << endl;
      cout << "momentum = (" << x << "," << y << ")" << endl;
      for (int i = 0; i < TotalSpace.GetHilbertSpaceDimension(); ++i)
	{
	  cout << i << " = ";
	  TotalSpace.PrintState(cout, i) << endl;
	}
      cout << endl << endl;
      for (int i = 0; i < TotalSpace.GetHilbertSpaceDimension(); ++i)
	{
	  cout << "---------------------------------------------" << endl;
	  cout << i << " = " << endl;;
	  for (int m1 = 0; m1 < MaxMomentum; ++m1)
	    for (int m2 = 0; m2 < m1; ++m2)
	      for (int m3 = 0; m3 < MaxMomentum; ++m3)
		{
		  int m4 = m1 + m2 - m3;
		  if (m4 < 0)
		    m4 += MaxMomentum;
		  else
		    if (m4 >= MaxMomentum)
		      m4 -= MaxMomentum;
		  if (m3 > m4)
		    {
		      double Coefficient = 0.0;
		      TotalSpace.AdAdAA(i, m1, m2, m3, m4, Coefficient);
		    }
		}
	}
/*	  
      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
      AbstractHamiltonian* Hamiltonian = new ParticleOnTorusCoulombHamiltonian (&TotalSpace, NbrFermions, MaxMomentum, XRatio);
      if (Hamiltonian->GetHilbertSpaceDimension() < MaxFullDiagonalization)
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
		  File << x << " " << y << " " << TmpTriDiag.DiagonalElement(j) << endl;
		  cout << x << " " << y << " " << TmpTriDiag.DiagonalElement(j) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      File << x << " " << y << " " << HRep(0, 0) << endl;
	    }
	}
      else
	{
	  FullReorthogonalizedLanczosAlgorithm Lanczos(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = NbrEigenvalue + 3;
	  Lanczos.SetHamiltonian(Hamiltonian);
	  Lanczos.InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  Lanczos.RunLanczosAlgorithm(NbrEigenvalue + 2);
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Precision > 1e-14) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
	    {
	      Lanczos.RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos.GetDiagonalizedMatrix());
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
	      File << x << " " << y << " " << TmpMatrix.DiagonalElement(i) << endl;
	    }
	  cout << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	}
      cout << "----------------------------------------------------------------" << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrFermions) << endl;
      delete Hamiltonian;
*/
    }
  File.close();

  return 0;
}
