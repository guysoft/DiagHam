#include "HilbertSpace/BosonOnLattice.h"

#include "Operator/ParticleOnLatticeOneBodyOperator.h"
#include "Operator/ParticleOnLatticeTranslationOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "MainTask/QHEOnLatticeMainTask.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeDensityMatrix" , "0.01");  
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "filenames of state vectors to be processed");

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 8);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 5);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 1);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice (-1=all)", -1);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  int NbrSites = Lx*Ly;  
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;

  int NbrVectors;
  char** VectorFiles = Manager.GetStrings("states",NbrVectors);

  if (NbrVectors==0)
    {
      cout << "At least one vector file is required!"<<endl;
      exit(1);
    }

  ParticleOnLattice* Space=new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

  ParticleOnLatticeOneBodyOperator *DensityOperator= new ParticleOnLatticeOneBodyOperator(Space);
  ParticleOnLatticeTranslationOperator *TranslationOperator= new ParticleOnLatticeTranslationOperator(Space);

  int VectorDimension = Space->GetHilbertSpaceDimension();
  ComplexVector *Vectors = new ComplexVector[NbrVectors];
  for (int i=0; i<NbrVectors; ++i)
    {
      Vectors[i].Resize(VectorDimension);
      Vectors[i].ReadVector(VectorFiles[i]);
      if (Vectors[i].GetVectorDimension()!=VectorDimension)
	{
	  cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of Hilbert-space!"<<endl;
	  exit(1);
	}
      //cout << "Vector "<<i<<":"<<endl<<Vectors[i]<<endl;
    }

  int DensityMatrixDimension = NbrSites*NbrVectors;
  HermitianMatrix Rho(DensityMatrixDimension);  

  Complex Tmp;
  int CreationIndex, AnnihilationIndex, TotalIndexI, TotalIndexJ;  
  for (int CreationX=0; CreationX<Lx; ++CreationX)
    for (int CreationY=0; CreationY<Ly; ++CreationY)
      {
	CreationIndex = Space->EncodeQuantumNumber(CreationX, CreationY, 0, Tmp);	
	for (int AnnihilationX=0; AnnihilationX<Lx; ++AnnihilationX)
	  for (int AnnihilationY=0; AnnihilationY<Ly; ++AnnihilationY)
	    {
	      AnnihilationIndex = Space->EncodeQuantumNumber(AnnihilationX, AnnihilationY, 0, Tmp);
	      
	      DensityOperator->SetCreationAnnihilationIndex(CreationIndex,AnnihilationIndex);
	      // calculate possible matrix elements in subspace of vectors
	      for (int numVector=0; numVector<NbrVectors; ++numVector)
		for (int numVector2=0; numVector2<NbrVectors; ++numVector2)
		  {
		    TotalIndexI = CreationIndex+numVector*NbrSites;
		    TotalIndexJ = AnnihilationIndex+numVector2*NbrSites;
		    if (TotalIndexI<=TotalIndexJ)
		      {
			Tmp=DensityOperator->MatrixElement(Vectors[numVector], Vectors[numVector2]);
			Rho.SetMatrixElement(TotalIndexI,TotalIndexJ,Tmp);
		      }
		  }
		  
	    }
      }
  // cout << "Matrix="<<endl<<Rho<<endl;
  // calculate eigenvalues & vectors of Rho
  RealDiagonalMatrix M;
  Rho.Diagonalize(M, 1e-10, 250);
  for (int i=0; i<DensityMatrixDimension; ++i)
    cout << "EV["<<i<<"] = " << M[i] << endl;

  

  RealVector TmpState(VectorDimension);
  for (int i=0; i<NbrVectors; ++i)
    {
      TranslationOperator->SetTranslationComponents(1,0);
      VectorOperatorMultiplyOperation Operation (TranslationOperator, &(Vectors[i]), &TmpState);      
      Operation.ApplyOperation(Architecture.GetArchitecture());      
      Complex Result1 = TmpState * Vectors[i];
      if (fabs(Norm(Result1)-1.0)>1e-10)
	{
	  cout << "State "<<VectorFiles[i]<< " is not a momentum K_x eigenstate (norm "<<Norm(Result1)<<")"<<endl;
	}
      else
	{
	  cout << "Momentum K_x of "<<VectorFiles[i]<<" = "<<Arg(Result1);
	  TranslationOperator->SetTranslationComponents(3,0);
	  VectorOperatorMultiplyOperation Operation (TranslationOperator, &(Vectors[i]), &TmpState);      
	  Operation.ApplyOperation(Architecture.GetArchitecture());      
	  Complex Result3 = TmpState * Vectors[i];
	  cout << " [ check: "<<Arg(Result3)/3<<" ]"<<endl;	  
	}
      TranslationOperator->SetTranslationComponents(0,1);
      VectorOperatorMultiplyOperation Operation2 (TranslationOperator, &(Vectors[i]), &TmpState);      
      Operation.ApplyOperation(Architecture.GetArchitecture());      
      Complex Result2 = TmpState * Vectors[i];
      if (fabs(Norm(Result1)-1.0)>1e-10)
	{
	  cout << "State "<<VectorFiles[i]<< " is not a momentum K_y eigenstate (norm "<<Norm(Result2)<<")"<<endl;
	}
      else
	{
	  cout << "Momentum K_x of "<<VectorFiles[i]<<" = "<<Arg(Result2);
	  TranslationOperator->SetTranslationComponents(0,3);
	  VectorOperatorMultiplyOperation Operation2 (TranslationOperator, &(Vectors[i]), &TmpState);      
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  Complex Result4 = TmpState * Vectors[i];
	  cout << " [ check: "<<Arg(Result4)/3<<" ]"<<endl;
	}      
    }

  delete DensityOperator;
  delete TranslationOperator;
}
