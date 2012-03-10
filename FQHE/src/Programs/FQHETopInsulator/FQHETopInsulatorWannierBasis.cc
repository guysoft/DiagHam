#include "HilbertSpace/BosonOnSquareLatticeWannierSpace.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/ParticleOnSquareLatticeWannierInterface.h"

#include "Vector/ComplexVector.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETopInsulatorWannierBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption ('\n', "state", "input file in a Wannier basis representation");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ky", "initial momentum along the y direction", 0);

  (*SystemGroup) += new SingleIntegerOption  ('\n', "target-ky", "target momentum on torus along the y direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "no-autodetect", "do not autdetect system parameter from state file name");

  (*SystemGroup) += new SingleStringOption ('\n', "torus", "a vector on the torus, with which to evaluate the overlap");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name to save the resulting vector, instead of replacement by torus conventions");
  
  (*MiscGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (not available for the non-periodic momentum space or the case with spin)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSiteX = Manager.GetInteger("nbr-sitex"); 
  int NbrSiteY = Manager.GetInteger("nbr-sitey"); 
  int TotalKy = Manager.GetInteger("ky");
  int TargetKy = Manager.GetInteger("target-ky");

  if (Manager.GetString("state") == 0)
    {
      cout << "an input state is required: use --state"<<endl;
      return -1;
    }
  if  (Manager.GetBoolean("no-autodetect") == false)
    {
      char * Filename = Manager.GetString("state");
      char * StrWannier = strstr(Filename, "Wannier");
      if (StrWannier == 0)
	StrWannier = strstr(Filename, "wannier");
      if (StrWannier==0)
	{
	  cout << "this does not appear to be a state in a wannier basis - use --no-autodetect to override"<<endl;
	  return -1;
	}
      double Mass = 0.0;
      int TmpKx=0; // not used for Wannier...
      bool Statistics = false;
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(Manager.GetString("state"),
							      NbrParticles, NbrSiteX, NbrSiteY, TmpKx, TotalKy, Mass, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
	  return -1;
	}
    }
 
  AbstractQHEParticle* Space;
  ParticleOnSquareLatticeWannierInterface *SpacePtr;
  
  if (Manager.GetBoolean("boson") == false)
    {
      cout << "Currently, no fermionic Wannier states are implemented!"<<endl;
      return -1;
    }
  else
    {
      if (Manager.GetInteger("nbr-subbands") == 1)
	{
	  Space = new BosonOnSquareLatticeWannierSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKy);
	  // two-step typecast to accommodate double inheritance properties
	  SpacePtr = (ParticleOnSquareLatticeWannierInterface*)((BosonOnSquareLatticeWannierSpace*)Space);
	  if (Manager.GetString("save-hilbert") != 0)
	    {
	      Space->WriteHilbertSpace(Manager.GetString("save-hilbert"));
	      return 0;
	    }
	}
      else
	{
	  cout << "Wannier states not yet implemented for multiple subbands." << endl;
	  return -1;
	}
    }

  ParticleOnTorus *TargetSpace = NULL;

  if (Manager.GetBoolean("boson"))
    TargetSpace = new BosonOnTorusShort(NbrParticles,NbrSiteX*NbrSiteY,TargetKy);
  else
    TargetSpace = new FermionOnTorus(NbrParticles,NbrSiteX*NbrSiteY,TargetKy);

  cout << "Dimension of Wannier space: "<<Space->GetHilbertSpaceDimension()<<endl;
  cout << "Dimension of torus space:   "<<TargetSpace->GetHilbertSpaceDimension()<<endl;
  
  int NbrTorusComponents = 0;
  double TorusProjection = 0.0;
  double WeightOtherComponents = 0.0;
  ComplexVector InputState;
  ComplexVector TorusState(TargetSpace->GetHilbertSpaceDimension(), true);
  if (InputState.ReadVector(Manager.GetString("state")) == false)
    {
      cout << "error while reading " << Manager.GetString("state") << endl;
      return -1;
    }
  if (Space->GetHilbertSpaceDimension() != InputState.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << InputState.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  int Index;
  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
    {
      //Space->PrintState(cout,i);
      //cout << " KY="<<SpacePtr->GetLinearizedMomentum(i)<<endl;
      if (SpacePtr->GetLinearizedMomentum(i)==TargetKy)
	{
	  Index = SpacePtr->ProjectToTorus(TargetSpace, i);
	  TorusState[Index] = InputState[i];
	  TorusProjection+=SqrNorm(InputState[i]);
	  ++NbrTorusComponents;
	}
      else
	WeightOtherComponents+=SqrNorm(InputState[i]);
	  
    }
  if (NbrTorusComponents!=TargetSpace->GetHilbertSpaceDimension())
    {
      cout << "Remark: Not all torus states found: "<<NbrTorusComponents<<" vs "<<TargetSpace->GetHilbertSpaceDimension()<<endl;
    }

  cout << "Weight of projection onto torus space: "<<TorusProjection<<endl;
  cout << "Weight outside of torus space:         "<<WeightOtherComponents<<endl;

  if (Manager.GetString("torus")!=NULL)
    {
      ComplexVector ReferenceState;
      if (ReferenceState.ReadVector(Manager.GetString("torus")) == false)
	{
	  cout << "error while reading " << Manager.GetString("torus") << endl;
	  return -1;
	}
      if (TargetSpace->GetHilbertSpaceDimension() != ReferenceState.GetVectorDimension())
	{
	  cout << "dimension mismatch between the state (" << ReferenceState.GetVectorDimension() << ") and the Hilbert space (" << TargetSpace->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}
      cout << "Overlap with torus eigenstate:         "<<ReferenceState*TorusState<<endl;
    }
  
  char *Output = Manager.GetString("output-file");
  if (Output==NULL)
    {      
      Output = AddSegmentInFileName(Manager.GetString("state"), "to_torus_", "annier_", false);
      if (Output == NULL)
	Output = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "torus_vec");
      if (Output == NULL)
	Output = AddExtensionToFileName(Manager.GetString("state"), "torus_vec");
    }
  
  TorusState.WriteVector(Output);
  
  delete Output;
  delete Space;
  delete TargetSpace;

  
  return 0;
}

