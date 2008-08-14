#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"

#include "Operator/ParticleOnLatticeOneBodyOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "MainTask/QHEOnLatticeMainTask.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "GeneralTools/FilenameTools.h"

#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

#include <bitset>
using std::bitset;


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;



class GutzwillerWaveFunction
{
public:
  // constructor
  // nbrParticles = particles in the condensate (should match space)
  // condensateWF = condensate wavefunction, in notation as single-particle hilbert-space
  // space = target-space of many-body state
  GutzwillerWaveFunction(int nbrParticles, ComplexVector& condensateWF, ParticleOnLattice *space);

  // destructor
  ~GutzwillerWaveFunction();

  // get Many-Body state
  // return = resultingState
  ComplexVector & GetGutzwillerWaveFunction();    

protected:

  // main recursion to calculate State (\sum \psi_i a^\dagger_i)^N |state>
  // exponent = power N remaining to be applied
  // state = state to be acted upon
  // prefactor = previous coefficients applied to state
  // in last stage of recursion, writes to this->TargetVector
  void Product (int exponent, unsigned long state, Complex prefactor);

  // target state, internal use
  ComplexVector TargetVector;

  // target condensate wavefunction
  ComplexVector CondensateState;

  // number of particles in many-body state
  int NbrParticles;

  // target Hilbert space
  ParticleOnLattice *Space;

  // number of states on lattice
  int NbrStates;

  // dimension of Space
  int Dim;
};


// constructor
// nbrParticles = particles in the condensate (should match space)
// condensateWF = condensate wavefunction, in notation as single-particle hilbert-space
// space = target-space of many-body state
GutzwillerWaveFunction::GutzwillerWaveFunction(int nbrParticles, ComplexVector &condensateWF, ParticleOnLattice *space)
{
  this->NbrParticles = nbrParticles;
  this->Space = space;
  this->CondensateState = ComplexVector(condensateWF,true);
  // believe this is the way it should be done to get CondensateState[q] associated with quantum no. q.
//   cout << "initial condensate: "<<CondensateState<<endl;
  this->CondensateState.ReverseVector();
//   cout << "reversed condensate: "<<CondensateState<<endl;
  this->Dim = Space->GetHilbertSpaceDimension();
  this->NbrStates = CondensateState.GetVectorDimension();
}

// destructor
GutzwillerWaveFunction::~GutzwillerWaveFunction()
{
}

// get Many-Body state
// return = resultingState
ComplexVector & GutzwillerWaveFunction::GetGutzwillerWaveFunction()
{
  this->TargetVector.Resize(Dim);
  this->TargetVector.ClearVector();
  // call main recursion
  this->Product(NbrParticles, 0x0l, 1.0);
  this->TargetVector/=this->TargetVector.Norm();
//   cout <<"Test norm: "<<TargetVector.Norm()<<endl;
  return this->TargetVector;
}

 // main recursion to calculate State (\sum \psi_i a^\dagger_i)^N |state>
// exponent = power N remaining to be applied
// state = state to be acted upon
// prefactor = previous coefficients applied to state
// in last stage of recursion, writes to this->TargetVector
void GutzwillerWaveFunction::Product (int exponent, unsigned long state, Complex prefactor)
{
  int Index;
  unsigned long ResultingState;
//   for (int i=NbrParticles; i>exponent; --i) cout <<"  ";
//   bitset<20> b = state;
//   cout << "Called: P("<<exponent<<", "<<b<<", prefactor: "<<prefactor<<endl;
  if (exponent>1)
    {
      for (int q=0; q<this->NbrStates; ++q)
	{
	  ResultingState = Space->Ad(state,q);
	  if (ResultingState!=0x0l)
	    Product(exponent-1, ResultingState, prefactor*CondensateState[q]);
	}
    }
  else
    {      
      for (int q=0; q<this->NbrStates; ++q)
	{
	  ResultingState = Space->Ad(state,q);
	  if (ResultingState!=0x0l)
	    {	      
	      if ((Index=Space->CarefulFindStateIndex(ResultingState,-1))<Dim)
		{
// 		  b=ResultingState;
// 		  for (int i=0; i<NbrParticles; ++i) cout <<"  ";
// 		  cout << "Finishing: with "<<b<<", index="<<Index<<", Prefactor: "<<prefactor*CondensateState[q]<<endl;
		  TargetVector[Index]+= prefactor*CondensateState[q];
		}
	    }
	}
    }
}


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeCondensateState" , "0.01");  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "condensate", "filename of vector describing condensate WF");

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles of target many-body state", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice (grabbed from condensate)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice (grabbed from condensate)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice (grabbed from condensate)", 0);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons (~Gutzwiller projection)");
  (*SystemGroup) += new BooleanOption('n',"no-hard-core","Do not use Hilbert-space of hard-core bosons (overriding detection from filename)");

  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
    
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  unsigned long MemorySpace = 9ul << 20;

  if (Manager.GetString("condensate")==0)
    {
      cout << "A vector for the condensate is required" << endl;
      exit(1);
    }

  double Interaction=0.0;
  int TmpI=-1;
  bool Statistics=false;
  bool HardCore=false;
  int NbrBosonsCondensate=0;

  if (FQHEOnLatticeFindSystemInfoFromVectorFileName(Manager.GetString("condensate"), NbrBosonsCondensate, Lx, Ly, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
    {
      cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
      exit(1);
    }

  HardCore=(HardCore||Manager.GetBoolean("hard-core"));

  if (Manager.GetBoolean("no-hard-core"))
    HardCore=false;

  if (NbrBosons==0) NbrBosons=NbrBosonsCondensate;

  cout << "Calculating Gutzwiller state for "<<NbrBosons<<" bosons on a lattice with geometry"<<endl;
  cout << "Lx="<<Lx<<", Ly="<<Ly<<", N_phi="<<NbrFluxQuanta;
  if (HardCore)
    cout <<" and hardcore interactions."<<endl;
  else
    cout << endl;

  ComplexVector CondensateState;

  if (CondensateState.ReadVector(Manager.GetString("condensate"))==false)
    {
      cout<<"Could not read condensate state!"<<endl;
      exit(1);
    }

  if (CondensateState.GetVectorDimension()!=Lx*Ly)
    {
      cout<<"Number of sites in condensate does not match the lattice Lx="<<Lx<<", Ly="<<Ly<<"!"<<endl;
      exit(1);
    }

  ParticleOnLattice* Space;
  if (HardCore)
    Space =new HardCoreBosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
  else Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);

  GutzwillerWaveFunction GutzwillerState(NbrBosons, CondensateState, Space);

  ComplexVector GutzwillerStateVector = GutzwillerState.GetGutzwillerWaveFunction();
  char extension[20];
  sprintf(extension,"N_%d.gw",NbrBosons);
  char *TmpFileName = AddExtensionToFileName(Manager.GetString("condensate"), extension);

  GutzwillerStateVector.WriteVector(TmpFileName);

  delete [] TmpFileName;
}
