#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/Options.h"

#include "MainTask/FQHEOnTorusMainTask.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <cstdio>
#include <fstream>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEFermionsTorusWithTranslation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);
  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "max-momentum", "maximum momentum for a single particle", 18);
  (*SystemGroup) += new SingleIntegerOption  ('x', "x-momentum", "constraint on the total momentum in the x direction (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "constraint on the total momentum in the y direction (negative if none)", -1);
  (*SystemGroup) += new SingleDoubleOption   ('R', "ratio", 
					      "ratio between lengths along the x and y directions (-1 if has to be taken equal to nbr-particles/4)", 
					      -1);
  (*SystemGroup) += new SingleIntegerOption  ('L', "landau-level", "Landau-level to be simulated", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "interaction-file", "file describing the interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "all-points", "calculate all points", false);
  (*SystemGroup) += new BooleanOption  ('\n', "no-wigner", "do not consider the energy contribution from the Wigner crystal", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsTorusWithTranslation -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int MaxMomentum = ((SingleIntegerOption*) Manager["max-momentum"])->GetInteger();
  int XMomentum = ((SingleIntegerOption*) Manager["x-momentum"])->GetInteger();
  int YMomentum = ((SingleIntegerOption*) Manager["y-momentum"])->GetInteger();
  char *LoadPrecalculationFile=Manager.GetString("load-precalculation");
  int LandauLevel=0;
  int NbrPseudopotentials=0;
  double *Pseudopotentials=NULL;
  double HaveCoulomb=false;
  char *InteractionName=NULL;
  if (Manager.GetString("interaction-file")!=NULL)
    {
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  exit(-1);
	}
      if (InteractionDefinition["CoulombLandauLevel"] != NULL)
	{
	  LandauLevel = atoi(InteractionDefinition["CoulombLandauLevel"]);
	  HaveCoulomb=true;
	}
      if (InteractionDefinition["Name"] == NULL)
	{
	  if ((InteractionDefinition["CoulombLandauLevel"] != NULL) && (InteractionDefinition["Pseudopotentials"] == NULL))
	    {
	      InteractionName = new char[18];
	      if (LandauLevel>=0)
		sprintf(InteractionName,"coulomb_l_%d",LandauLevel);
	      else
		sprintf(InteractionName,"graphene_l_%d",-LandauLevel);
	    }
	  else
	    {
	      cout << "Attention, using unnamed interaction! Please include a line 'Name = ...'" << endl;
	      InteractionName = new char[10];
	      sprintf(InteractionName,"unnamed");
	    }
	}
      else
	{
	  InteractionName = new char[strlen(InteractionDefinition["Name"])+1];
	  strcpy(InteractionName, InteractionDefinition["Name"]);
	}
      InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', Pseudopotentials, NbrPseudopotentials);
    }
  else
    {
      LandauLevel = Manager.GetInteger("landau-level");
      InteractionName = new char[1];
      InteractionName[0]='\0';
      HaveCoulomb=true;
    }
  double XRatio = NbrFermions / 4.0;
  if (((SingleDoubleOption*) Manager["ratio"])->GetDouble() > 0)
    {
      XRatio = ((SingleDoubleOption*) Manager["ratio"])->GetDouble();
    }
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* OutputName = new char [512];
  if (NbrPseudopotentials>0)
    {
      sprintf (OutputName, "fermions_torus_%s_n_%d_2s_%d_ratio_%f.dat", InteractionName, NbrFermions, MaxMomentum, XRatio);
    }
  else
    {
      if (LandauLevel>0)
	sprintf (OutputName, "fermions_torus_coulomb_l_%d_n_%d_2s_%d_ratio_%f.dat", LandauLevel, NbrFermions, MaxMomentum, XRatio);
      else
	if (LandauLevel<0)
	  sprintf (OutputName, "fermions_torus_graphene_l_%d_n_%d_2s_%d_ratio_%f.dat", -LandauLevel, NbrFermions, MaxMomentum, XRatio);
	else
	  sprintf (OutputName, "fermions_torus_coulomb_n_%d_2s_%d_ratio_%f.dat", NbrFermions, MaxMomentum, XRatio);
    }

  int MomentumModulo = FindGCD(NbrFermions, MaxMomentum);
  int XMaxMomentum = (MomentumModulo - 1);
  bool GenerateMomenta = false;
  if ((XMomentum < 0)||(YMomentum < 0))
    GenerateMomenta = true;
  if (XMomentum < 0)
    XMomentum = 0;
  else
    XMaxMomentum = XMomentum;
  int YMaxMomentum = (MaxMomentum - 1);
  if (YMomentum < 0)
    YMomentum = 0;
  else
    YMaxMomentum = YMomentum;

  int NbrMomenta;
  int *XMomenta;
  int *YMomenta;
  int *Multiplicities = NULL;
  int CenterX=0, CenterY=0;

  if (GenerateMomenta==false)
    {
      NbrMomenta=1;
      XMomenta = new int[1];
      YMomenta = new int[1];
      XMomenta[0]=XMomentum;
      YMomenta[0]=YMomentum;
    }
  else
    {
      if (Manager.GetBoolean("all-points"))
	{
	  int Pos=0;
	  NbrMomenta = (XMaxMomentum-XMomentum+1)*(YMaxMomentum-YMomentum+1);
	  XMomenta = new int[NbrMomenta];
	  YMomenta = new int[NbrMomenta];
	  for (; XMomentum <= XMaxMomentum; ++XMomentum)
	    for (int YMomentum2 = YMomentum; YMomentum2<= YMaxMomentum; ++YMomentum2)
	      {
		XMomenta[Pos]=XMomentum;
		YMomenta[Pos]=YMomentum2;
		++Pos;
		cout << "Pos="<<Pos<<endl;
	      }
	}
      else // determine inequivalent states in BZ
	{
	  if (NbrFermions&1)
	    {
	      CenterX=0;
	      CenterY=0;
	    }
	  else
	    {
	      if ((NbrFermions/MomentumModulo*MaxMomentum/MomentumModulo)&1) // p*q odd?
		{
		  CenterX=MomentumModulo/2;
		  CenterY=MomentumModulo/2;
		}
	      else
		{
		  CenterX=0;
		  CenterY=0;
		}
	    }
	  if (XRatio == 1.0)
	    {
	      NbrMomenta=0;
	      for (int Kx = CenterX; Kx<=CenterX+MomentumModulo/2; ++Kx)
		for (int Ky= (Kx-CenterX)+CenterY; Ky<=CenterY+MomentumModulo/2; ++Ky)
		  {
		    ++NbrMomenta;
		  }
	      int Pos=0;
	      XMomenta = new int[NbrMomenta];
	      YMomenta = new int[NbrMomenta];
	      Multiplicities = new int[NbrMomenta];
	      for (int Kx = 0; Kx<=MomentumModulo/2; ++Kx)
		for (int Ky= Kx; Ky<=MomentumModulo/2; ++Ky, ++Pos)
		  {
		    XMomenta[Pos]=CenterX+Kx;
		    YMomenta[Pos]=CenterY+Ky;
		    if (Kx==0)
		      {
			if (Ky==0)
			  Multiplicities[Pos]=1; // BZ center
			else if (Ky==MomentumModulo/2)
			  Multiplicities[Pos]=2;
			else Multiplicities[Pos]=4;
		      }
		    else if (Kx==MomentumModulo/2)
		      {
			Multiplicities[Pos]=1; // BZ corner
		      }
		    else
		      {
			if (Ky==Kx) // diagonal ?
			  {
			    Multiplicities[Pos]=4; 
			  }
			else
			  {
			    if (Ky==MomentumModulo/2)
			      Multiplicities[Pos]=4;
			    else
			      Multiplicities[Pos]=8;
			  }
		      }
		  }
	    }
	  else // rectangular torus
	    {
	      NbrMomenta=(MomentumModulo/2+1)*(MomentumModulo/2+1);
	      int Pos=0;
	      XMomenta = new int[NbrMomenta];
	      YMomenta = new int[NbrMomenta];
	      Multiplicities = new int[NbrMomenta];
	      for (int Kx = 0; Kx<=MomentumModulo/2; ++Kx)
		for (int Ky= 0; Ky<=MomentumModulo/2; ++Ky, ++Pos)
		  {
		    XMomenta[Pos]=CenterX+Kx;
		    YMomenta[Pos]=CenterY+Ky;
		    if (Kx==0)
		      {
			if (Ky==0)
			  Multiplicities[Pos]=1; // BZ center
			else // on Gamma->X]
			  Multiplicities[Pos]=2;
		      }
		    else
		      {
			if (Ky==0)
			  Multiplicities[Pos]=2;
			else
			  {
			    if (Kx==MomentumModulo/2)
			      {
				if (Ky==MomentumModulo/2) // BZ corner?
				  Multiplicities[Pos]=1;
				else
				  Multiplicities[Pos]=2;
			      }
			    else
			      {
				if (Ky==MomentumModulo/2) // on edge?
				  Multiplicities[Pos]=2;
				else
				  Multiplicities[Pos]=4;
			      }
			  }
		      }
		  }
	    }
	}
    }
  
  bool FirstRun=true;
  for (int Pos=0;Pos<NbrMomenta; ++Pos)
    {
      XMomentum=XMomenta[Pos];
      YMomentum=YMomenta[Pos];
      
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
      //	FermionOnTorus TotalSpace (NbrFermions, MaxMomentum, y);
      FermionOnTorusWithMagneticTranslations *TotalSpace = new FermionOnTorusWithMagneticTranslations(NbrFermions, MaxMomentum, XMomentum, YMomentum);
      //cout << " Total Hilbert space dimension = " << TotalSpace->GetHilbertSpaceDimension() << endl;
      //cout << "momentum = (" << XMomentum << "," << YMomentum << ")" << endl;
      Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());

      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian (TotalSpace, 
													   NbrFermions, MaxMomentum, XMomentum, XRatio, HaveCoulomb, LandauLevel, NbrPseudopotentials, Pseudopotentials, Manager.GetBoolean("no-wigner"),
													   Architecture.GetArchitecture(), 
													   Memory, LoadPrecalculationFile);
      
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate"))	
	{
	  EigenvectorName = new char [512];
	  char *TmpName=RemoveExtensionFromFileName(OutputName, ".dat");
	  sprintf (EigenvectorName, "%s_kx_%d_ky_%d", TmpName, XMomentum, YMomentum);
	  delete [] TmpName;
	}
      double Shift=0.0;
      FQHEOnTorusMainTask Task (&Manager, TotalSpace, Hamiltonian, YMomentum, Shift, OutputName, FirstRun, EigenvectorName);
      Task.SetKxValue(XMomentum);
      if (Multiplicities!=0)
	Task.SetMultiplicity(Multiplicities[Pos]);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      if (FirstRun == true)
	FirstRun = false;
      delete Hamiltonian;
      delete TotalSpace;
    }

  delete [] XMomenta;
  delete [] YMomenta;
  if (Multiplicities!=0)
    delete [] Multiplicities;
  return 0;
}
