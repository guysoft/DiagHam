#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/QuantumDot/Spectra/Spectra.h"
#include "Tools/QuantumDot/Spectra/OverlapSpectra.h"
#include "Tools/QuantumDot/Spectra/AverageSpectra.h"
#include "Tools/QuantumDot/Spectra/DOSSpectra.h"
#include "Tools/QuantumDot/Spectra/PeriodicSpectra.h"
#include "Tools/QuantumDot/Spectra/XYReflexionSymmetricPeriodicSpectra.h"
#include "Tools/QuantumDot/Spectra/HardBoxSpectra.h"
#include "Tools/QuantumDot/Spectra/CylinderInMagneticFieldSpectra.h"
#include "Tools/QuantumDot/Spectra/CylinderQuantumDotSpectra.h"

#include "HilbertSpace/QuantumDotHilbertSpace/Confined3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/XYReflexionSymmetricPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/VerticalPeriodicParticleInMagneticField.h"

#include "Tools/QuantumDot/Potential/QuantumDotThreeDConstantCylinderPotential.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;

int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("DefaultCylinderSymmetryOscllatorForce" , "0.01");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* FileGroup =  new OptionGroup ("File options");
  OptionGroup* SampleGroup =  new OptionGroup ("Sample options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += HilbertSpaceGroup;
  Manager += FileGroup;
  Manager += SampleGroup;
  Manager += MiscGroup;

  (*HilbertSpaceGroup) += new SingleIntegerOption ('R', "R-states", "number of states in plane", 50);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('Z', "Z-states", "number of cells in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('m', "momentum", "quantum number of kinetic momentum in the plane of the second state", 0);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('b', "magnetic", "magnetic field in Z direction (in Tesla unit)", 30);

  (*FileGroup) += new SingleStringOption ('\n', "state", "name of the file containing the first state", "eigenvector.0");
  (*FileGroup) += new SingleStringOption ('\n', "energy", "name of the file energy of the first state", "eigenvalues");
  (*FileGroup) += new SingleIntegerOption ('\n', "last", "number of the last destination state", 1);

  (*SampleGroup) += new SingleDoubleOption ('\n', "sizeZ", "size of the sample (in Angstrom unit)", 500);

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");
  

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrStateR = ((SingleIntegerOption*) Manager["R-states"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["Z-states"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  int NumberM = ((SingleIntegerOption*) Manager["momentum"])->GetInteger();
  double MagneticField = ((SingleDoubleOption*) Manager["magnetic"])->GetDouble();

  char* State = ((SingleStringOption*) Manager["state"])->GetString();
  char* Energy = ((SingleStringOption*) Manager["energy"])->GetString();
  int Last = ((SingleIntegerOption*) Manager["last"])->GetInteger();
 
  double SizeZ = ((SingleDoubleOption*) Manager["sizeZ"])->GetDouble();

   // define Hilbert space   
  VerticalPeriodicParticleInMagneticField* Space = new VerticalPeriodicParticleInMagneticField(0, NbrStateR, NbrStateZ, LowImpulsionZ); 

  VerticalPeriodicParticleInMagneticField* Space2 = new VerticalPeriodicParticleInMagneticField(NumberM, NbrStateR, NbrStateZ, LowImpulsionZ); 

  CylinderInMagneticFieldSpectra spectra (Space, State, MagneticField);

  double ReX, ImX, ReY, ImY, ReZ, ImZ;
  // polarization in the Z direction
  if (NumberM == 0)
    {
      ofstream polarization ("Polarization.txt");  
      ofstream PX("PolarizationX.txt");
      ofstream PY("PolarizationY.txt");  
      ofstream PZ("PolarizationZ.txt");
      ifstream energy ("eigenvalues");

      char** Files = new char* [Last + 1];
      double tmpE;  
      double fundamental;
      energy >> fundamental;
      for (int i = 1; i <= Last; ++i)
	{
	  Files[i] = new char[80];
	  AddString(Files[i], "eigenvector.", i, "");
	  spectra.GetMeanPosition(Space2, Files[i], SizeZ, ReX, ImX, ReY, ImY, ReZ, ImZ);
	  energy >> tmpE;
	  polarization << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) * (tmpE - fundamental) << '\t' << ((ReY * ReY) + (ImY * ImY)) * (tmpE - fundamental) << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) * (tmpE - fundamental) << '\n';
	  PX << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) * (tmpE - fundamental) << '\n';
	  PY << tmpE - fundamental << '\t' << ((ReY * ReY) + (ImY * ImY)) * (tmpE - fundamental) << '\n';
	  PZ << tmpE - fundamental << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) * (tmpE - fundamental) << '\n';	  
	}
      PX.close(); PY.close(); PZ.close();
      energy.close(); polarization.close();
      delete[] Files;
    }
  // polarization in the plane
  else
    {
      ofstream polarization ("Polarization.txt");  
      ofstream PX("PolarizationX.txt");
      ofstream PY("PolarizationY.txt");  
      ofstream PZ("PolarizationZ.txt");
      ifstream energybis (Energy);
      double fundamental;
      energybis >> fundamental;
      energybis.close();
      ifstream energy ("eigenvalues");

      char** Files = new char* [Last + 1];
      double tmpE;  
      for (int i = 0; i <= Last; ++i)
	{
	  Files[i] = new char[80];
	  AddString(Files[i], "eigenvector.", i, "");
	  spectra.GetMeanPosition(Space2, Files[i], SizeZ, ReX, ImX, ReY, ImY, ReZ, ImZ);
	  energy >> tmpE;
	  polarization << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) * (tmpE - fundamental) << '\t' << ((ReY * ReY) + (ImY * ImY)) * (tmpE - fundamental) << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) * (tmpE - fundamental) << '\n';
	  PX << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) * (tmpE - fundamental) << '\n';
	  PY << tmpE - fundamental << '\t' << ((ReY * ReY) + (ImY * ImY)) * (tmpE - fundamental) << '\n';
	  PZ << tmpE - fundamental << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) * (tmpE - fundamental) << '\n';	  
	}
      PX.close(); PY.close(); PZ.close();
      energy.close(); polarization.close();
      delete[] Files;
    }

  return 1;
}
