#include "Options/Options.h"

#include "Vector/RealVector.h"
 
//#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
//#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "GeneralTools/FilenameTools.h"

#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Hamiltonian/ParticleOnSphereTwoLandauLevelDeltaHamiltonian.h"

#include "Tools/FQHESpectrum/PseudoPotentials.h"
#include "Tools/FQHESpectrum/AbstractZDensityProfile.h"

#define		my_min(X, Y)  ((X) < (Y) ? (X) : (Y))

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <stdlib.h>
#include <math.h>

using std::ofstream;
using std::ios;
using std::cout;
using std::endl;
using std::string;

int main(int argc, char** argv)
{
  OptionManager Manager ("Delta2LLPseudopotentials" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (i.e. twice the maximum momentum for a single particle on LLL)", 8);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-type", "type of interaction to generate pseudopotentials for (delta or coulomb support at pressent). Default: delta", "delta");
  
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Delta2LLPseudopotentials -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");
  char *InteractionType = Manager.GetString("interaction-type");
  
  int LzMaxUp = NbrFluxQuanta + 2;
  int LzMaxDown = NbrFluxQuanta + 1;
  char* OutputFile = Manager.GetFormattedString("pseudopotential_%interaction-type%_s_%nbr-flux%.dat");
  
  if ( strcmp(InteractionType, "delta") != 0 && strcmp(InteractionType, "coulomb") != 0 )
    {
      cout << InteractionType << " unsupported at pressent. Please choose either delta or coulomb." << endl;
      return -1;
    }
  
  ClebschGordanCoefficients ClebschDownDown (LzMaxDown - 1, LzMaxDown - 1);
  ClebschGordanCoefficients ClebschUpUp (LzMaxUp, LzMaxUp);
  ClebschGordanCoefficients ClebschUpDown (LzMaxUp, LzMaxDown - 1);
  ClebschGordanCoefficients ClebshDownUp (LzMaxDown - 1, LzMaxUp);
  
  // these are the labels of the arrays as they will be in the file.
  string PseudoLabels[9] = {"PseudopotentialsUpUpUpUp","PseudopotentialsUpUpDownDown","PseudopotentialsUpUpUpDown",
			    "PseudopotentialsDownDownUpUp","PseudopotentialsDownDownDownDown","PseudopotentialsDownDownUpDown",
			    "PseudopotentialsUpDownUpUp","PseudopotentialsUpDownDownDown","PseudopotentialsUpDownUpDown"};
  
  // these are the lenghts of the arrays corresponding to the labels above. 			    
  int PseudoLengths[9] = { LzMaxUp + 1, LzMaxUp - 1 , LzMaxUp - 1, LzMaxUp - 1, LzMaxUp - 1, LzMaxUp - 2, LzMaxUp - 1, LzMaxUp - 2, LzMaxUp - 1}; 
      
  if ( strcmp(InteractionType, "delta") == 0 ) 
    {
      ofstream File;
      File.open(OutputFile, ios::binary | ios::out);
      File.precision(14);
      
      double **PseudoPotentialArrays;
      PseudoPotentialArrays = Evaluate2LLSphereDeltaPseudopotentials(NbrFluxQuanta, true);
    
      for ( int j = 0; j < 9 ; j++ ) 
	{
	  File << PseudoLabels[j]<< "=";	
	  for (int i = 0; i < PseudoLengths[j] ; ++i)
	      File << " " << PseudoPotentialArrays[j][i];	
	  File << endl ;
      }
	    
      File << endl;      
      File.close();
      
      for ( int i = 0 ; i < 9 ; i++ ) 
	{
	  delete [] PseudoPotentialArrays[i];
	}
      delete [] PseudoPotentialArrays;      
    }
  else if ( strcmp(InteractionType, "coulomb") == 0 ) 
    {
      ofstream File;
      File.open(OutputFile, ios::binary | ios::out);
      File.precision(14);
      
      double **PseudoPotentialArrays;
      PseudoPotentialArrays = Evaluate2LLSphereCoulombPseudopotentials(NbrFluxQuanta, true);
           	  
      for ( int j = 0; j < 9 ; j++ )
	{
	File << PseudoLabels[j]<< "=";        
	  for (int i = 0; i < PseudoLengths[j] ; ++i)
	      File << " " << PseudoPotentialArrays[j][i];		  
	  File << endl ;
        }            	     
	    
      File << endl;      
      File.close();
      
      for ( int i = 0 ; i < 9 ; i++ ) 
	{
	  delete [] PseudoPotentialArrays[i];
	}
      delete [] PseudoPotentialArrays;  
      
    }
    
  delete[] OutputFile;
    
  return 0;
}