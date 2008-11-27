#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/FermionOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"
#include "Hamiltonian/ParticleOnLatticeDeltaHamiltonian.h"
#include "Operator/ParticleOnLatticeTranslationOperator.h"

#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "MainTask/QHEOnLatticeMainTask.h"

#include "MathTools/IntegerAlgebraTools.h"
#include "MathTools/JacobiThetaFunction.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// some global variables to avoid too many parameters in function calls
int NbrBosons;
int Lx;
int Ly;
ArchitectureManager Architecture;

// values = eigenvalues
// vectors = complex matrix with many vectors to be analyzed and overwritten with basis diagonal in k
// NbrFlux = number of flux quanta to be considered
// start, end = indices of range of eigenvalues to be considered
ComplexMatrix& DiagonalizeMomentaInSubspace(RealDiagonalMatrix &values, ComplexMatrix &vectors, ParticleOnLatticeTranslationOperator *translationOperator, int NbrFlux, int start, int end);

void GetTranslationMatrix(ParticleOnLatticeTranslationOperator *Operator, int NbrVectors,
			  ComplexVector *Vectors, ComplexMatrix &MatrixRepresentation,
			  ComplexVector &TmpState, ArchitectureManager &Architecture);


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeCompositeFermions" , "0.01");  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");  

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);  
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 3);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 6);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 9);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice ", 6);
  (*SystemGroup) += new SingleIntegerOption  ('f', "flux-per-CF", "number of flux attached to each boson (allowed values: +/-1)", 1);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");
  
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  NbrBosons = Manager.GetInteger("nbr-particles");
  Lx = Manager.GetInteger("lx");
  Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  int CFFlux = Manager.GetInteger("flux-per-CF");
  bool HardCore = Manager.GetBoolean("hard-core");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;

  char* OutputName;
  char interactionStr[20]="";
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      OutputName = new char [256];      
      if (HardCore)
	sprintf(interactionStr,"_hardcore");      
      sprintf (OutputName, "bosons_lattice_CF_n_%d_x_%d_y_%d%s_q_%d_p_%d", NbrBosons, Lx, Ly, interactionStr, NbrFluxQuanta, CFFlux);
    }
  char *TmpC = new char[strlen(OutputName)+20];

  cout << "* Full Hilbert-space: N="<<NbrBosons<<" bosons in "<<Lx<<" x "<<Ly<<" cells at N_phi="<<NbrFluxQuanta<<endl;
  ParticleOnLattice* Space;
  if (HardCore)
    Space =new HardCoreBosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
  else Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  
  int AttachedFlux = CFFlux * NbrBosons;  
  ParticleOnLatticeTranslationOperator *TranslationOperator;
  
  // constructing 1P states:
  
  cout << "* CF states contribute "<<NbrFluxQuanta-AttachedFlux<<" flux"<<endl;
  // space in which CF's live (statistics doesn't matter as we consider single particle physics!)
  //BosonOnLattice *CFSpace = new BosonOnLattice(/*NbrParticles*/ 1, Lx, Ly, NbrFluxQuanta-AttachedFlux, MemorySpace);
  FermionOnLattice *CFSpace = new FermionOnLattice(/*NbrParticles*/ 1, Lx, Ly, NbrFluxQuanta-AttachedFlux, MemorySpace);
  TranslationOperator = new ParticleOnLatticeTranslationOperator(CFSpace);
  
  // corresponding Hamiltonians
  AbstractQHEOnLatticeHamiltonian* CFHamiltonian = new ParticleOnLatticeDeltaHamiltonian(CFSpace, /*NbrParticles*/ 1, Lx, Ly, NbrFluxQuanta-AttachedFlux, /* U */ 0.0 , /*ReverseHopping*/ false, /* Delta */ 0.0, /* Random */ 0.0, Architecture.GetArchitecture(), 0, NULL);    
  
  HermitianMatrix HCF(CFHamiltonian->GetHilbertSpaceDimension(), true);
  ComplexMatrix CFEigenVecs(CFHamiltonian->GetHilbertSpaceDimension(), CFHamiltonian->GetHilbertSpaceDimension());
  CFHamiltonian->GetHamiltonian(HCF);
  RealDiagonalMatrix CFEigenVals(CFHamiltonian->GetHilbertSpaceDimension());
  HCF.Diagonalize(CFEigenVals, CFEigenVecs, /* error */ 1e-10 , /* maxIter */ 250);
  ParticleOnLatticeTranslationOperator *CFTranslationOperator= new ParticleOnLatticeTranslationOperator(CFSpace);
  DiagonalizeMomentaInSubspace(CFEigenVals, CFEigenVecs,  CFTranslationOperator, NbrFluxQuanta-AttachedFlux, 0, NbrBosons);
  delete CFTranslationOperator;

  for (int i=0; i<NbrBosons; ++i)
    cout << "E_CF["<<i<<"]="<<CFEigenVals[i]<<" norm of EVec: "<<CFEigenVecs[i].Norm()<<endl;
  cout << "E_other["<<NbrBosons<<"]="<<CFEigenVals[NbrBosons]<<endl;

  for (int i=0; i<NbrBosons; ++i)
    {
      sprintf(TmpC,"%s.CF.%d.vec",OutputName,i);
      CFEigenVecs[i].WriteVector(TmpC);
    }


  cout << "* LLL states for Jastrow-factor contribute "<<AttachedFlux<<" flux"<<endl;  
  
  
  // calculate states required to build Jastrow factor:
  // BosonOnLattice *JastrowSpace = new BosonOnLattice(/*NbrParticles*/ 1, Lx, Ly, AttachedFlux, MemorySpace);
  FermionOnLattice *JastrowSpace = new FermionOnLattice(/*NbrParticles*/ 1, Lx, Ly, AttachedFlux, MemorySpace);

  AbstractQHEOnLatticeHamiltonian* JastrowHamiltonian = new ParticleOnLatticeDeltaHamiltonian(JastrowSpace, /*NbrParticles*/ 1, Lx, Ly, AttachedFlux, /* U */ 0.0 , /*ReverseHopping*/ false, /* Delta */ 0.0, /* Random */ 0.0, Architecture.GetArchitecture(), 0, NULL);
  delete TranslationOperator;
  TranslationOperator = new ParticleOnLatticeTranslationOperator(JastrowSpace);
  
  HermitianMatrix HJastrow(JastrowHamiltonian->GetHilbertSpaceDimension(), true);
  ComplexMatrix JastrowEigenVecs(JastrowHamiltonian->GetHilbertSpaceDimension(),
				 JastrowHamiltonian->GetHilbertSpaceDimension());
  JastrowHamiltonian->GetHamiltonian(HJastrow);
  RealDiagonalMatrix JastrowEigenVals(JastrowHamiltonian->GetHilbertSpaceDimension());
  HJastrow.Diagonalize(JastrowEigenVals, JastrowEigenVecs, /* error */ 1e-10 , /* maxIter */ 250);
  ParticleOnLatticeTranslationOperator *JastrowTranslationOperator= new ParticleOnLatticeTranslationOperator(JastrowSpace);
  DiagonalizeMomentaInSubspace(JastrowEigenVals, JastrowEigenVecs,  JastrowTranslationOperator, AttachedFlux, 0, NbrBosons);
  delete JastrowTranslationOperator;

  for (int i=0; i<NbrBosons; ++i)
    cout << "E_Jastrow["<<i<<"]="<<JastrowEigenVals[i]<<" norm of EVec: "<<JastrowEigenVecs[i].Norm()<<endl;
  cout << "E_other["<<NbrBosons<<"]="<<JastrowEigenVals[NbrBosons]<<endl;

  for (int i=0; i<NbrBosons; ++i)
    {
      sprintf(TmpC,"%s.Jastrow.%d.vec",OutputName,i);
      JastrowEigenVecs[i].WriteVector(TmpC);      
    }

  // for calculation of analytic Laughlin state:
  JacobiThetaFunction ThetaRel(0.5,0.5,Complex(0.0,((double)Ly)/Lx));
  JacobiThetaFunction ThetaCM1(0.0/2.0+(NbrFluxQuanta-2.0)/4.0,1.0/2.0-(NbrFluxQuanta-2.0)/2.0,Complex(0.0,(2.0*(double)Ly)/Lx));
  JacobiThetaFunction ThetaCM2(NbrFluxQuanta/4.0,(2.0-NbrFluxQuanta)/2.0,Complex(0.0,(2.0*(double)Ly)/Lx));

  ComplexVector AnalyticJastrow(Space->GetHilbertSpaceDimension(), true);
  ComplexVector AnalyticRelative(Space->GetHilbertSpaceDimension(), true);
  ComplexVector AnalyticCM1(Space->GetHilbertSpaceDimension(), true);
  ComplexVector AnalyticCM2(Space->GetHilbertSpaceDimension(), true);
  ComplexVector Analytic1(Space->GetHilbertSpaceDimension(), true);
  ComplexVector Analytic2(Space->GetHilbertSpaceDimension(), true);

  int *PosX = new int[NbrBosons];
  int *PosY = new int[NbrBosons];
  int Subl, SumX, SumY, SumSqY;
  Complex FRel;  
  
  // cycle through all configurations of the Hilbert-space and calculate the corresponding Slater determinants

  int *QuantumNumbers = new int[NbrBosons];
  double Multiplicity;

#ifdef __LAPACK__
  ComplexLapackDeterminant SlaterCF(NbrBosons);
  ComplexLapackDeterminant SlaterJastrow(NbrBosons);
#else
  ComplexMatrix SlaterCF(NbrBosons, NbrBosons);
  ComplexMatrix SlaterJastrow(NbrBosons, NbrBosons);
#endif
  Complex Value;
  Complex Tmp;

  ComplexVector CFState(Space->GetHilbertSpaceDimension(), true);
  ComplexVector JastrowState(Space->GetHilbertSpaceDimension(), true);
  
  for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
    {
      Space->ListQuantumNumbers(i,QuantumNumbers,Multiplicity);
      //cout << "Q=["<<QuantumNumbers[0];
      //for (int k=1; k<NbrBosons; ++k) cout<<", "<<QuantumNumbers[k];
      //cout << "]"<<endl;
      if (Multiplicity<1.01) // ignore state with multiple occupation
	{	  
	  for (int p = 0; p < NbrBosons; ++p)
	    {
// 	      ComplexVector TmpEigenFctCF = CFEigenVecs[p];
// 	      ComplexVector TmpEigenFctJastrow = JastrowEigenVecs[p];
	      for (int q = 0; q < NbrBosons; ++q)
		{
		  // need to consider proper ordering of matrix elements
		  // in Hilbert-space, largest quantum number q corresponds to position 0!
		  SlaterCF.SetMatrixElement(p,q,CFEigenVecs[p][CFHamiltonian->GetHilbertSpaceDimension()-1-QuantumNumbers[q]]);
		  SlaterJastrow.SetMatrixElement(p,q,JastrowEigenVecs[p][JastrowHamiltonian->GetHilbertSpaceDimension()-1-QuantumNumbers[q]]);
		}	      
	    }
	  CFState[i] = SlaterCF.Determinant();
	  JastrowState[i] = SlaterJastrow.Determinant();
	  
	  // analytic states:
	  SumX=0; SumY=0, SumSqY=0;
	  for (int q = 0; q < NbrBosons; ++q)
	    {
	      Space->DecodeQuantumNumber(QuantumNumbers[q],PosX[q],PosY[q],Subl);
	      SumX+=PosX[q];
	      SumY+=PosY[q];
	      SumSqY+=PosY[q]*PosY[q];
	    }

	  FRel=1.0;
	  for (int bi = 1; bi < NbrBosons; ++bi)
	    for (int bj = 0; bj < bi; ++bj)
	      FRel*=ThetaRel.GetValue(Complex( ((double)(PosX[bi]-PosX[bj]))/Lx,((double)(PosY[bi]-PosY[bj]))/Lx));
	  
	  AnalyticJastrow[i] = FRel * exp(-0.5*SumSqY);
	  AnalyticRelative[i] = FRel*FRel * exp(-0.5*SumSqY);
	  AnalyticCM1[i] = ThetaCM1.GetValue(Complex((2.0*(double)SumX)/Lx,(2.0*(double)SumY)/Lx));
	  AnalyticCM2[i] = ThetaCM2.GetValue(Complex((2.0*(double)SumX)/Lx,(2.0*(double)SumY)/Lx));
 	  ThetaCM1.PrintValue(cout,Complex((2.0*(double)SumX)/Lx,(2.0*(double)SumY)/Lx))<<endl;
	  ThetaCM2.PrintValue(cout,Complex((2.0*(double)SumX)/Lx,(2.0*(double)SumY)/Lx))<<endl;
	  Analytic1[i] = AnalyticRelative[i] * AnalyticCM1[i];
	  Analytic2[i] = AnalyticRelative[i] * AnalyticCM2[i];
	}
    }
  CFState /= CFState.Norm();  
  sprintf(TmpC,"%s.CF.vec",OutputName);
  CFState.WriteVector(TmpC);

  JastrowState /= JastrowState.Norm();
  sprintf(TmpC,"%s.J.vec",OutputName);
  JastrowState.WriteVector(TmpC);
  
  sprintf(TmpC,"%s.vec",OutputName);
  for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
    CFState[i] *= JastrowState[i];
  CFState /= CFState.Norm();
  CFState.WriteVector(TmpC);

  AnalyticJastrow /= AnalyticJastrow.Norm();
  sprintf(TmpC,"%s.AJ.vec",OutputName);
  AnalyticJastrow.WriteVector(TmpC);

  AnalyticRelative /= AnalyticRelative.Norm();
  sprintf(TmpC,"%s.ACF.vec",OutputName);
  AnalyticRelative.WriteVector(TmpC);

  AnalyticCM1 /= AnalyticCM1.Norm();
  sprintf(TmpC,"%s.CM1.vec",OutputName);
  AnalyticCM1.WriteVector(TmpC);

  AnalyticCM2 /= AnalyticCM2.Norm();
  sprintf(TmpC,"%s.CM2.vec",OutputName);
  AnalyticCM2.WriteVector(TmpC);

  Analytic1 /= Analytic1.Norm();
  sprintf(TmpC,"%s.A1.vec",OutputName);
  Analytic1.WriteVector(TmpC);

  Analytic2 /= Analytic2.Norm();
  sprintf(TmpC,"%s.A2.vec",OutputName);
  Analytic2.WriteVector(TmpC);

  
  delete [] TmpC;
  delete [] QuantumNumbers;
  
  delete CFHamiltonian;
  delete CFSpace;  

  delete JastrowHamiltonian;
  delete JastrowSpace;
  
  return 0;
}



// values = eigenvalues
// vectors = complex matrix with many vectors to be analyzed and overwritten with basis diagonal in k
// NbrFlux = number of flux quanta to be considered
// start, end = indices of range of eigenvalues of the same eigenvalue
// based on code from FQHELatticeDensityMatrix
ComplexMatrix& DiagonalizeMomentaAtEnergy(ComplexMatrix &vectors, ParticleOnLatticeTranslationOperator *translationOperator, int nbrFlux, int start, int end)
{
  int NbrVectors = end - start;
  ComplexMatrix XTranslationMatrix(NbrVectors, NbrVectors);
  ComplexMatrix YTranslationMatrix(NbrVectors, NbrVectors);
  ComplexVector TmpState(vectors.GetNbrRow());

  ComplexVector TmpState3(vectors.GetNbrRow());

  ComplexVector *Vectors = new ComplexVector[NbrVectors];
  for (int i=start, j=0; i<end; ++i,++j)
    Vectors[j]=ComplexVector(vectors[i], true);

  int Degeneracy=1;
  int n1=1, n2=1;
  int FluxModulo = FindGCD(nbrFlux, Lx*Ly);
  int r=nbrFlux/FluxModulo;
  int t=Lx*Ly/FluxModulo;

  while ((((Ly*n1)%t)!=0) && (n1<Lx)) ++n1;
  while ((((Lx*n2)%t)!=0) && (n2<Ly)) ++n2;

  while ((r*NbrBosons*n1*n2*Degeneracy)%t != 0) ++Degeneracy;
  
  cout << "N_phi = "<<r<<"/"<<t<<endl;
  cout << "n1="<<n1<<", n2="<<n2<<", global degeneracy: "<<Degeneracy<<endl;

  ComplexMatrix EVecX(NbrVectors, NbrVectors);
  ComplexMatrix EVecY(NbrVectors, NbrVectors);  
  ComplexDiagonalMatrix EValX(NbrVectors, NbrVectors);
  ComplexDiagonalMatrix EValY(NbrVectors, NbrVectors);
  
  translationOperator->SetTranslationComponents(n1,0);
  GetTranslationMatrix(translationOperator, NbrVectors, Vectors, XTranslationMatrix, TmpState, Architecture);

  translationOperator->SetTranslationComponents(0,n2);
  GetTranslationMatrix(translationOperator, NbrVectors, Vectors, YTranslationMatrix, TmpState, Architecture);

  ComplexMatrix EVecXY(NbrVectors, NbrVectors);


  // form linear superposition of Tx and Ty to diagonalize:
  ComplexMatrix Z((Matrix&)XTranslationMatrix);
  Z*=log(91.0); // scale with some random number > 1
  Z+=YTranslationMatrix;
  Z.Diagonalize(EValX,EVecXY);
  ComplexMatrix QH=EVecXY.GetAdjoint();

  for (int i=0;i<NbrVectors;++i)
    {
      vectors[start+i].ClearVector();
      for (int j=0; j<NbrVectors;++j)
	vectors[start+i].AddLinearCombination(Conj(EVecXY[i][j]),Vectors[j]);      
    }
    
  delete [] Vectors;
  return vectors;
}


// values = eigenvalues
// vectors = complex matrix with many vectors to be analyzed and overwritten with basis diagonal in k
// nbrFlux = number of flux quanta to be considered
// start, end = indices of range of eigenvalues to be considered
ComplexMatrix& DiagonalizeMomentaInSubspace(RealDiagonalMatrix &values, ComplexMatrix &vectors, ParticleOnLatticeTranslationOperator *translationOperator, int nbrFlux, int start, int end)

{
  // make sure we get entire multiplets:
  while ((start > 0)&&(abs(values[start]-values[start-1])<1e-12)) --start;
  while ((end < values.GetNbrRow()-1)&&(fabs(values[end]-values[end-1])<1e-12))
    {
      ++end;
    }
  
  int startSeg=start, endSeg=start+1;

  while (endSeg <= end)
    {
      while ((endSeg<end)&&(fabs(values[endSeg]-values[endSeg-1])<1e-12))
	++endSeg;      
      if (endSeg-startSeg > 1)
	{
	  cout << "Multiplet ["<<startSeg<<", "<<endSeg-1<<"]"<<endl;
	  DiagonalizeMomentaAtEnergy(vectors, translationOperator, nbrFlux, startSeg, endSeg);
	}
      startSeg=endSeg;
      endSeg=startSeg+1;
    }
  return vectors;
}


void GetTranslationMatrix(ParticleOnLatticeTranslationOperator *Operator, int NbrVectors,
			  ComplexVector *Vectors, ComplexMatrix &MatrixRepresentation,
			  ComplexVector &TmpState, ArchitectureManager &Architecture)
{
  Complex Tmp;
  for (int i=0; i<NbrVectors; ++i)
    {
      VectorOperatorMultiplyOperation Operation (Operator, &(Vectors[i]), &TmpState);      
      Operation.ApplyOperation(Architecture.GetArchitecture());           
      for (int j=0; j<NbrVectors; ++j)
	{
	  Tmp = Vectors[j] * TmpState;
	  MatrixRepresentation.SetMatrixElement(i,j,Tmp);
	}
    }
}
