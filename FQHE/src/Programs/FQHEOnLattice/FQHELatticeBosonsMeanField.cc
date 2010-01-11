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

  (*OptimizationGroup) += new BooleanOption('\n', "gradient", "Use optimization based on gradients");
  (*OptimizationGroup) += new SingleDoubleOption('\n', "tolerance", "tolerance for variational parameters in condensate",1e-6);
  (*OptimizationGroup) += new SingleIntegerOption('i', "nbr-iter", "number of iterations for optimization",10000);
  (*OptimizationGroup) += new SingleStringOption('\n', "parameters", "file with initial parameters");

  (*OptimizationGroup) += new SingleIntegerOption('a', "nbr-attempts", "number of separate attempts to optimize a configuration",1);
  (*OptimizationGroup) += new SingleIntegerOption('s', "nbr-save", "maximum number of (distinct) configurations to be saved",10);
  (*OptimizationGroup) += new SingleDoubleOption('\n', "overlap-same", "overlap above which two conf's are considered identical", 0.99);
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file", NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  double ChemicalPotential = Manager.GetDouble("mu");

  int NbrAttempts = Manager.GetInteger("nbr-attempts");
  int NbrToSave = Manager.GetInteger("nbr-save");
  double IdentityThreshold = Manager.GetDouble("overlap-same");
  
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

  ComplexVector *InitialParameters = NULL;

  if ((Manager.GetString("parameters")!=NULL)&&(NbrAttempts==1))
    {
      InitialParameters = new ComplexVector;
      if (InitialParameters->ReadVector(Manager.GetString("parameters"))==false)
	{
	  cout << "Could not read vector of initial parameters" <<Manager.GetString("parameters")<<endl;
	  exit(1);
	}
    }

  ofstream File;
  ifstream TestFile;
  TestFile.open(OutputName, ios::in);
  if (TestFile.is_open())
    {
      TestFile.close();
      File.open(OutputName, ios::app );
    }
  else
    {
      File.open(OutputName, ios::out );
      File << "#E_tot\tDensity\tE_int\tParameters"<<endl;
    }
    
  ComplexVector *OptimalWaveFunctions = new ComplexVector[NbrToSave];
  int NbrFound=0;
  double *LowestEnergies = new double[NbrToSave];
  
  for (int i=0; i<NbrAttempts; ++i)
    {
      GrossPitaevskiiOnLatticeState MeanFieldState(NbrSites, Manager.GetString("potential-file"), Manager.GetString("interaction-file"), Lattice, InitialParameters);
      MeanFieldState.SetChemicalPotential(ChemicalPotential);
      if (InitialParameters==NULL)
	MeanFieldState.SetToRandomPhase();
      int MaxEval = 2*NbrSites*Manager.GetInteger("nbr-iter");
      double Energy;
      if (Manager.GetBoolean("gradient"))
	Energy=MeanFieldState.GradientOptimize(Manager.GetDouble("tolerance"), MaxEval, /*initialStep*/ 0.01, /*lineMinPar*/ Manager.GetDouble("tolerance")/10.0);
      else
	Energy=MeanFieldState.Optimize(Manager.GetDouble("tolerance"), MaxEval);
      RealVector Parameters=MeanFieldState.GetVariationalParameters();
      bool Recorded=false;
      ComplexVector TmpWaveFunction=MeanFieldState.GetWaveFunction();
      for (int k=0; (k<NbrFound) && (Recorded==false); ++k)
	if (Energy <= LowestEnergies[k])
	  {
	    if (Norm(TmpWaveFunction*OptimalWaveFunctions[k])>IdentityThreshold*TmpWaveFunction.Norm()*OptimalWaveFunctions[k].Norm())
	      { // same configuration: simply replace with the one of lower energy
		OptimalWaveFunctions[k]=TmpWaveFunction;
		LowestEnergies[k] = Energy;
		Recorded=true;
	      }
	    else
	      {
		int UpperLimit;
		if (NbrFound < NbrToSave)
		  {
		    LowestEnergies[NbrFound]=LowestEnergies[NbrFound-1];
		    OptimalWaveFunctions[NbrFound]=OptimalWaveFunctions[NbrFound-1];
		    UpperLimit=NbrFound-1;
		    NbrFound++;
		  }
		else UpperLimit = NbrToSave-1;
		for (int s=UpperLimit; s>k; --s)
		  {
		    LowestEnergies[s]=LowestEnergies[s-1];
		    OptimalWaveFunctions[s]=OptimalWaveFunctions[s-1];
		  }
		OptimalWaveFunctions[k]=TmpWaveFunction;
		LowestEnergies[k] = Energy;
		Recorded=true;
	      }
	    }
	else
	  {
	    if (Norm(TmpWaveFunction*OptimalWaveFunctions[k])>IdentityThreshold*TmpWaveFunction.Norm()*OptimalWaveFunctions[k].Norm())
	      Recorded=true;
	  }
      if ((Recorded==false)&&(NbrFound<NbrToSave))
	{
	  OptimalWaveFunctions[NbrFound]=MeanFieldState.GetWaveFunction();
	  LowestEnergies[NbrFound] = Energy;
	  ++NbrFound;
	}
      
      cout << "Found mean field state with energy: "<<Energy<<" and density "<< MeanFieldState.GetNbrParticles()/NbrSites<<endl;

      File.precision(10);
      File << Energy << "\t" << MeanFieldState.GetNbrParticles()/NbrSites << "\t" << Energy/MeanFieldState.GetNbrParticles() + ChemicalPotential;
      File.precision(nearbyint(-log(Manager.GetDouble("tolerance"))/log(10)));
      for (int i=0; i<Parameters.GetVectorDimension(); ++i)
	File << "\t" << Parameters[i];
      File << endl;
    }
  File.close();

  if (NbrAttempts>1)
    {
      char* SelectedName = AddExtensionToFileName(BaseName,"select.dat");
      ofstream SelectFile ( SelectedName, ios::app );
      int CellPosition[2];
      int Sub;
      RealVector SitePosition;

      for (int i=0; i<NbrFound; ++i)
	{
	  int Counter=0;
	  char* ParameterName = GetUniqueFileName(BaseName,Counter,".par");
	  char* FieldName = ReplaceExtensionToFileName(ParameterName,"par","wf");
	  ofstream VectorField;
	  VectorField.open(FieldName,ios::out);
	  if (Counter==0)
	    SelectFile << "#Count\tE_tot\tFilename"<<endl;
	  SelectFile.precision(12);
	  SelectFile << Counter << "\t" << LowestEnergies[i] << "\t" << ParameterName << endl;
	  OptimalWaveFunctions[i].WriteVector(ParameterName);
	  for (int s=0; s<NbrSites; ++s)
	    {
	      Lattice->GetSiteCoordinates(s, CellPosition, Sub);
	      SitePosition = Lattice->GetSitePosition(CellPosition,Sub);
	      VectorField << SitePosition[0] << "\t" << SitePosition[1]
			  << "\t" << (OptimalWaveFunctions[i])[s].Re
			  << "\t" << (OptimalWaveFunctions[i])[s].Im << endl;
	    }
	  VectorField.close();
	  delete [] ParameterName;
	}
      SelectFile.close();
      delete [] SelectedName;
    }

  delete [] BaseName;
  delete [] OutputName;
  delete Lattice;
  delete [] LatticeName;
  if (InitialParameters != NULL)
    delete InitialParameters;
  return 0;
}
