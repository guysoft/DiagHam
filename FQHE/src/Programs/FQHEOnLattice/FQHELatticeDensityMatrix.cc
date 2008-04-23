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
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice", 0);
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
      cout << "Vector "<<i<<":"<<endl;
      for (int j=0; j<VectorDimension; ++j)
	{
	  cout<<Vectors[i][j]<<" ( ";
	  Space->PrintState(cout,j);
	  cout << " )"<<endl;
	}
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

  

  ComplexVector TmpState(VectorDimension);
  ComplexMatrix XTranslationMatrix(NbrVectors, NbrVectors);
  ComplexMatrix YTranslationMatrix(NbrVectors, NbrVectors);  
  for (int i=0; i<NbrVectors; ++i)
    for (int j=0; j<NbrVectors; ++j)
      {
	TranslationOperator->SetTranslationComponents(1,0);
	VectorOperatorMultiplyOperation Operation (TranslationOperator, &(Vectors[i]), &TmpState);      
	Operation.ApplyOperation(Architecture.GetArchitecture());      
	Tmp = TmpState * Vectors[j];
	XTranslationMatrix.SetMatrixElement(i,j,Tmp);

	TranslationOperator->SetTranslationComponents(0,1);
	VectorOperatorMultiplyOperation Operation2 (TranslationOperator, &(Vectors[i]), &TmpState);      
	Operation2.ApplyOperation(Architecture.GetArchitecture());      
	Tmp = TmpState * Vectors[j];
	YTranslationMatrix.SetMatrixElement(i,j,Tmp);
      }

  cout << "XTranslationMatrix="<<endl<<XTranslationMatrix<<endl;
  cout << "YTranslationMatrix="<<endl<<YTranslationMatrix<<endl;
  
  for (int i=0; i<NbrVectors; ++i)
    {
      for (int dX=1; dX<=Lx; dX++)
	{
	  TranslationOperator->SetTranslationComponents(dX,0);
	  VectorOperatorMultiplyOperation Operation (TranslationOperator, &(Vectors[i]), &TmpState);      
	  Operation.ApplyOperation(Architecture.GetArchitecture());      
	  Complex Result1 = TmpState * Vectors[i];
	  if (fabs(Norm(Result1)-1.0)>1e-10)
	    {
	      cout << "dX= "<<dX<<": State "<<VectorFiles[i]<< " is not a momentum K_x eigenstate (norm "<<Norm(Result1)<<")"<<endl;
	    }
	  else
	    {
	      cout << "Momentum K_x from translation "<<dX<<" of "<<VectorFiles[i]<<" = "<<Arg(Result1)/2.0/M_PI*Lx<<"/"<<Lx<<endl;
	    }	  
	}
      for (int dY=1; dY<=Ly; dY++)
	{
	  TranslationOperator->SetTranslationComponents(0,dY);
	  VectorOperatorMultiplyOperation Operation2 (TranslationOperator, &(Vectors[i]), &TmpState);      
	  Operation2.ApplyOperation(Architecture.GetArchitecture());      
	  Complex Result2 = TmpState * Vectors[i];
	  if (fabs(Norm(Result2)-1.0)>1e-10)
	    {
	      cout << "dY= "<<dY<<": State "<<VectorFiles[i]<< " is not a momentum K_y eigenstate (norm "<<Norm(Result2)<<")"<<endl;
	    }
	  else
	    {
	      cout << "Momentum K_y from translation "<<dY<<" of "<<VectorFiles[i]<<" = "<<Arg(Result2)/2.0/M_PI*Ly<<"/"<<Ly<<endl;
	    }	     
	}
    }
  
  delete DensityOperator;
  delete TranslationOperator;
}
