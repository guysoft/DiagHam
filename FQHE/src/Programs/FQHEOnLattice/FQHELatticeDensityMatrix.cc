#include "HilbertSpace/BosonOnLattice.h"

#include "Operator/ParticleOnLatticeOneBodyOperator.h"
#include "Operator/ParticleOnLatticeTranslationOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "MainTask/QHEOnLatticeMainTask.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MathTools/IntegerAlgebraTools.h"

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
      // cout << "Vector "<<i<<":"<<endl;
//       for (int j=0; j<VectorDimension; ++j)
// 	{
// 	  cout<<Vectors[i][j]<<" ( ";
// 	  Space->PrintState(cout,j);
// 	  cout << " )"<<endl;
// 	}
      
//       int dX, dY;
//       for (int d=1;d<VectorDimension;++d)
// 	{
// 	  if (fabs(Norm(Vectors[i][d])-Norm(Vectors[i][0]))<1e-10)
// 	    {
// 	      cout << "Potential translation invariance 0->"<<d<<endl;	      
// 	      if (Space->IsTranslation(0,d,dX,dY))
// 		cout << "States have same particle numbers"<<endl;
// 	    }
// 	}
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
  ComplexVector TmpState2(VectorDimension);
  ComplexMatrix TrRep(VectorDimension, VectorDimension);
  // testing unitarity of translation operator matrix:
  
  TranslationOperator->SetTranslationComponents(1,0);
  for (int i=0; i<VectorDimension; ++i)
    {
      TmpState2.ClearVector();
      TmpState2.Re(i)=1.0;
      VectorOperatorMultiplyOperation Operation (TranslationOperator, &TmpState2, &TmpState);      
      Operation.ApplyOperation(Architecture.GetArchitecture());      
      for (int j=0; j<VectorDimension; ++j)
	TrRep.SetMatrixElement(j,i,TmpState[j]);
    }

  cout << "Representation of T_x"<<endl<<TrRep<<endl;

  
  ComplexMatrix XTranslationMatrix(NbrVectors, NbrVectors);
  ComplexMatrix YTranslationMatrix(NbrVectors, NbrVectors);
  ComplexVector TmpState3(VectorDimension);

  int Degeneracy=1;
  int n1=1, n2=1;
  int FluxModulo = FindGCD(NbrFluxQuanta, Lx*Ly);
  int r=NbrFluxQuanta/FluxModulo;
  int t=Lx*Ly/FluxModulo;

  while ((((Ly*n1)%t)!=0) && (n1<Lx)) ++n1;
  while ((((Lx*n2)%t)!=0) && (n2<Ly)) ++n2;

  if ((Lx%n1)!=0)
    cout << "Extending range of n1 to Lx"<<endl;
  if ((Ly%n2)!=0)
    cout << "Extending range of n2 to Ly"<<endl;

  if (((n1*n2*NbrFluxQuanta)%t) != 0)
    {
      cout << "Cannot resolve translations: Brillouin zone trivial?"<<endl;
      n1=Lx;
      n2=Ly;
    }

  while ((r*NbrBosons*n1*n2*Degeneracy)%t != 0) ++Degeneracy;
  
  cout << "N_phi = "<<r<<"/"<<t<<endl;
  cout << "n1="<<n1<<", n2="<<n2<<", global degeneracy: "<<Degeneracy<<endl;
  
  for (int i=0; i<NbrVectors; ++i)
    {
      // cout<<"Vector["<<i<<"]="<<endl;
//       for (int k=0; k<VectorDimension; ++k)
// 	cout << Vectors[i].Re(k) << "+"<<Vectors[i].Im(k)<<"i"<<endl;      
      TranslationOperator->SetTranslationComponents(n1,0);
      VectorOperatorMultiplyOperation Operation (TranslationOperator, &(Vectors[i]), &TmpState);      
      Operation.ApplyOperation(Architecture.GetArchitecture());
      
      // cout<<"T_x*Vector["<<i<<"]="<<endl;
//       for (int k=0; k<VectorDimension; ++k)
// 	cout << TmpState.Re(k) << "+"<<TmpState.Im(k)<<"i"<<endl;
      
      for (int j=0; j<NbrVectors; ++j)
	{
	  Tmp = Vectors[j] * TmpState;
	  XTranslationMatrix.SetMatrixElement(i,j,Tmp);
	}
	
      TranslationOperator->SetTranslationComponents(0,n2);
      VectorOperatorMultiplyOperation Operation2 (TranslationOperator, &(Vectors[i]), &TmpState);      
      Operation2.ApplyOperation(Architecture.GetArchitecture());
      for (int j=0; j<NbrVectors; ++j)
	{
	  Tmp = Vectors[j] * TmpState;
	  YTranslationMatrix.SetMatrixElement(i,j,Tmp);
	}
      
    }

  cout << "XTranslationMatrix="<<endl<<XTranslationMatrix<<endl;
  cout << "YTranslationMatrix="<<endl<<YTranslationMatrix<<endl;

  ComplexMatrix XT2((Matrix&)XTranslationMatrix);
  
  
  ComplexMatrix EVecX(NbrVectors, NbrVectors);
  ComplexMatrix EVecY(NbrVectors, NbrVectors);
  ComplexDiagonalMatrix EValX(NbrVectors, NbrVectors);
  ComplexDiagonalMatrix EValY(NbrVectors, NbrVectors);
  
  XTranslationMatrix.Diagonalize(EValX);
  cout << "X Eigenvalues="<<endl<<EValX<<endl;
  
  XTranslationMatrix.Diagonalize(EValX,EVecX);
  YTranslationMatrix.Diagonalize(EValY,EVecY);

  cout << "X Eigenvalues="<<endl<<EValX<<endl;
  cout << "X Eigenvectors="<<endl<<EVecX<<endl;
  cout << "Y Eigenvalues="<<endl<<EValY<<endl;
  cout << "Y Eigenvectors="<<endl<<EVecY<<endl;
  
  delete DensityOperator;
  delete TranslationOperator;
}
