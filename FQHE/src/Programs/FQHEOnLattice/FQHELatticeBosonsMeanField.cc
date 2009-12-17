#include "Tools/FQHESpectrum/LatticePhases.h"
#include "Tools/FQHEWaveFunction/GrossPitaevskiiOnLatticeState.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeBosonsMeanField" , "0.01");  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OptimizationGroup = new OptionGroup ("optimization options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  LatticePhases::AddOptionGroup(&Manager);
  Manager += OptimizationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleDoubleOption  ('m', "mu", "chemical potential", 1.0);
  (*SystemGroup) += new SingleStringOption  ('e', "interaction-file", "use definition of two-body interactions from a file");
  (*SystemGroup) += new SingleStringOption  ('E', "interaction-name", "descriptor of external interaction (if in use)","ext");
  (*SystemGroup) += new SingleStringOption  ('f', "potential-file", "use definition of one-body interactions from a file");
  (*SystemGroup) += new SingleStringOption  ('F', "potential-name", "descriptor of external single particle potential (if in use)","");
  
  (*OptimizationGroup) += new SingleDoubleOption('\n', "tolerance", "tolerance for variational parameters in condensate",1e-6);
  (*OptimizationGroup) += new SingleIntegerOption('i', "nbr-iter", "number of iterations for optimization",10000);
  (*OptimizationGroup) += new SingleStringOption('\n', "parameters", "file with initial parameters");
  
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  double ChemicalPotential = Manager.GetDouble("mu");

  if ((Manager.GetString("interaction-file")==NULL)&&(Manager.GetString("potential-file")==NULL))
    {
      cout << "An external definition of the hopping or interaction matrix elements is required. Use option -e and/or -f"<<endl;
      exit(1);
    }

  // get the lattice geometry
  LatticePhases *Lattice = new LatticePhases();
  char* LatticeName = Lattice->GeometryString();
  int NbrSites = Lattice->GetNbrSites();
  
  char* OutputName;
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      char* TmpChar = new char[1024];
      OutputName = new char [1024];
      sprintf(OutputName,"bosons_lattice_mf_%s",LatticeName);
      if (Manager.GetString("interaction-file")!=NULL)
	{
	  if (Manager.GetString("interaction-name")!=NULL)
	    sprintf(TmpChar,"%s_%s",OutputName,Manager.GetString("interaction-name"));
	  else
	    sprintf(TmpChar,"%s_%s",OutputName,Manager.GetString("interaction-file"));
	  strcpy(OutputName,TmpChar);
	}
      if (Manager.GetString("potential-file")!=NULL)
	{
	  if (Manager.GetString("potential-name")!=NULL)
	    sprintf(TmpChar,"%s_%s",OutputName,Manager.GetString("potential-name"));
	  else
	    sprintf(TmpChar,"%s_%s",OutputName,Manager.GetString("potential-file"));
	  strcpy(OutputName,TmpChar);
	}
      sprintf(TmpChar,"%s_mu_%g.dat",OutputName,ChemicalPotential);
      strcpy(OutputName,TmpChar);
      delete [] TmpChar;
    }

  char* BaseName = RemoveExtensionFromFileName(OutputName,".dat");
  int Counter=0;
  char* ParameterName = GetUniqueFileName(BaseName,Counter,".par");

  RealVector *InitialParameters = NULL;

  if (Manager.GetString("parameters")!=NULL)
    {
      InitialParameters = new RealVector;
      if (InitialParameters->ReadVector(Manager.GetString("parameters"))==false)
	{
	  cout << "Could not read vector of initial parameters" <<Manager.GetString("parameters")<<endl;
	  exit(1);
	}
    }
  
  GrossPitaevskiiOnLatticeState MeanFieldState(NbrSites, Manager.GetString("potential-file"), Manager.GetString("interaction-file"), Lattice, InitialParameters);
  MeanFieldState.SetChemicalPotential(ChemicalPotential);
  if (InitialParameters==NULL)
    MeanFieldState.SetToRandomPhase();
  int MaxEval = 2*NbrSites*Manager.GetInteger("nbr-iter");
  double Energy=MeanFieldState.Optimize(Manager.GetDouble("tolerance"), MaxEval);
  RealVector Parameters=MeanFieldState.GetVariationalParameters();
  Parameters.WriteVector(ParameterName);
  cout << "Found mean field state with energy: "<<Energy<<" and density "<< MeanFieldState.GetNbrParticles()/NbrSites<<endl<<ParameterName<<endl;
  ofstream File ( OutputName, ios::app );
  File.precision(10);
  if (Counter==0)
    File << "#Count\tE_tot\tDensity\tE_int\tParameters"<<endl;
  File << Counter << "\t" << Energy << "\t" << MeanFieldState.GetNbrParticles()/NbrSites << "\t" << Energy/MeanFieldState.GetNbrParticles() - ChemicalPotential;
  File.precision(nearbyint(-log(Manager.GetDouble("tolerance"))/log(10)));
  for (int i=0; i<Parameters.GetVectorDimension(); ++i)
    File << "\t" << Parameters[i];
  File << "\t" << ParameterName << endl;
  delete [] BaseName;
  delete [] OutputName;
  delete [] ParameterName;
  delete Lattice;
  delete [] LatticeName;
  if (InitialParameters != NULL)
    delete InitialParameters;
  return 0;
}
