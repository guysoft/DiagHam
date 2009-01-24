#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "Tools/FQHEMonteCarlo/AbstractMCSamplingFunction.h"
#include "Tools/FQHEMonteCarlo/QHESamplingFunctionManager.h"
#include "Tools/FQHEMonteCarlo/SimpleMonteCarloOnSphereAlgorithm.h"
#include "Tools/FQHEMonteCarlo/SphereCoulombEnergy.h"
#include "Tools/FQHEMonteCarlo/SphereGeneralEnergy.h"

#include "MCObservables/RealObservable.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/Endian.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereMonteCarlo" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");      
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
    
  QHEWaveFunctionManager WaveFunctionManager(QHEWaveFunctionManager::SphereGeometry);
  QHESamplingFunctionManager SamplingFunctionManager(QHESamplingFunctionManager::SphereGeometry);
  
  Manager += SystemGroup;
  SimpleMonteCarloOnSphereAlgorithm::AddOptionGroup(&Manager);
  WaveFunctionManager.AddOptionGroup(&Manager);
  SamplingFunctionManager.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum angular momentum per particle", 15);
  (*SystemGroup) += new SingleStringOption('\n',"interaction-params","File containing parameters of interaction");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new BooleanOption ('\n', "list-samplingfunctions", "list all available sampling-fuctions");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout); // new feature to process options and display help in 1 line!

  if (Manager.GetBoolean("list-wavefunctions") == true)
    {
      WaveFunctionManager.ShowAvalaibleWaveFunctions(cout);
      return 0;
    }
  
  SamplingFunctionManager.TestForShow();
  

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int NbrFlux = Manager.GetInteger("lzmax");
  
  Abstract1DComplexFunction* TestWaveFunction = WaveFunctionManager.GetWaveFunction();
  
  AbstractMCSamplingFunction* SamplingFunction = SamplingFunctionManager.GetSamplingFunction();

  cout << "Function: " << WaveFunctionManager.GetDescription()<<endl;
  cout << "Sampler:  " << SamplingFunctionManager.GetDescription()<<endl;

  SimpleMonteCarloOnSphereAlgorithm MonteCarloRoutine(NbrParticles, TestWaveFunction, SamplingFunction,
						      &Manager);
  
  if (Manager.GetString("interaction-params")==0)
    {
      // add Coulomb energy in LLL
      SphereCoulombEnergy *Energy=new SphereCoulombEnergy(/*Flux */ NbrFlux);
      MonteCarloRoutine.AddObservable(Energy);
    }
  else
    {
      // add a general energy function 
      SphereGeneralEnergy *Energy=new SphereGeneralEnergy(NbrFlux, Manager.GetString("interaction-params"));
      MonteCarloRoutine.AddObservable(Energy);
    }

  // run simulation
  MonteCarloRoutine.Simulate();

  // print final results:
  cout << "Final results:" << endl;
  MonteCarloRoutine.WriteObservations(cout);  
  
}

