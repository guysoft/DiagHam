#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"

#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"


#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include <iomanip>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
//new namespaces (added by ba340)
using namespace std;
using std::setw;

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FCIHofstadterCorrelation" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PlotOptionGroup = new OptionGroup ("plot options");  
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += PlotOptionGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the vector file describing the state whose density has to be plotted");
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density instead of density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "k-space", "compute the density/correlation in momentum space", false);

  (*PlotOptionGroup) += new SingleStringOption ('\n', "output", "output file name (default output name replace the .vec extension of the input file with .rho or .rhorho)", 0);
  (*PlotOptionGroup) += new SingleIntegerOption ('\n', "nbr-samplesx", "number of samples along the x direction", 100, true, 10);
  (*PlotOptionGroup) += new SingleIntegerOption ('\n', "nbr-samplesy", "number of samples along the y direction", 100, true, 10);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCICorrelation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("state") == 0)
    {
      cout << "FCICorrelation requires an input state" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state")) == false)
    {
      cout << "can't find vector file " << Manager.GetString("state") << endl;
      return -1;      
    }

  int NbrParticles = 0;
  int NbrCellX = 0;
  int NbrCellY = 0;
  int MomentumX = 0;
  int MomentumY = 0;
  int UnitCellX=0;
  int UnitCellY=0;     
  int FluxPerCell=0;   
  char Axis='y';         
  double Interaction=0; 
  double GammaX=0;
  double GammaY=0;
  bool EmbeddingFlag=false;
  bool Hardcore=false;	
  int NbrState=0; 
  int NbrSamplesX = Manager.GetInteger("nbr-samplesx");
  int NbrSamplesY = Manager.GetInteger("nbr-samplesy");
  bool DensityFlag = Manager.GetBoolean("density");
  bool Statistics = true;

  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName_Hofstadter(Manager.GetString("state"), NbrParticles, NbrCellX, NbrCellY, Interaction, FluxPerCell, NbrState, Statistics, Hardcore, EmbeddingFlag, Axis, GammaX, GammaY, MomentumX, MomentumY, UnitCellX, UnitCellY) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
      return -1;
    }

  cout << setw(20) << left << "Statistics" << setw(20) << left << Statistics << endl;
  cout << setw(20) << left << "UnitCellX" << setw(20) << left << UnitCellX << endl;
  cout << setw(20) << left << "UnitCellY" << setw(20) << left << UnitCellY << endl;
  cout << setw(20) << left << "FluxPerCell" << setw(20) << left << FluxPerCell << endl;
  cout << setw(20) << left << "Axis" << setw(20) << left << Axis << endl;
  cout << setw(20) << left << "NbrParticles" << setw(20) << left << NbrParticles << endl;
  cout << setw(20) << left << "NbrCellX" << setw(20) << left << NbrCellX << endl;
  cout << setw(20) << left << "NbrCellY" << setw(20) << left << NbrCellY << endl;
  cout << setw(20) << left << "Hardcore" << setw(20) << left << Hardcore << endl;
  cout << setw(20) << left << "Interaction" << setw(20) << left << Interaction << endl;
  cout << setw(20) << left << "GammaX" << setw(20) << left << GammaX << endl;
  cout << setw(20) << left << "GammaY" << setw(20) << left << GammaY << endl;
  cout << setw(20) << left << "EmbeddingFlag" << setw(20) << left << EmbeddingFlag << endl;
  cout << setw(20) << left << "MomentumX" << setw(20) << left << MomentumX << endl;
  cout << setw(20) << left << "MomentumY" << setw(20) << left << MomentumY << endl;
  cout << setw(20) << left << "NbrState" << setw(20) << left << NbrState << endl;

  ParticleOnSphere* Space = 0;
  if (Statistics == true)
    Space = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrCellX, NbrCellY, MomentumX, MomentumY);
  else
    Space = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrCellX, NbrCellY, MomentumX, MomentumY);
  ComplexVector ComplexState;
  if (ComplexState.ReadVector (Manager.GetString("state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state") << endl;
      return -1;      
    }
  Complex* PrecalculatedValues_rho = 0;
  Complex* PrecalculatedValues_rhorho = 0;
  int* PrecalculatedIndices = 0;
  int NbrPrecalculatedValues = 0;
  if (DensityFlag == false)
    {
      for (int kx1 =0; kx1 < NbrCellX; ++kx1)
	for (int kx2 =0; kx2 < NbrCellX; ++kx2)
	  for (int kx3 =0; kx3 < NbrCellX; ++kx3)
	    for (int kx4 =0; kx4 < NbrCellX; ++kx4)
	      {
		if (((kx1 + kx2 - kx3 - kx4) % NbrCellX) == 0)
		  {
		    for (int ky1 = 0; ky1 < NbrCellY; ++ky1)
		      for (int ky2 = 0; ky2 < NbrCellY; ++ky2)
			for (int ky3 = 0; ky3 < NbrCellY; ++ky3)
			  for (int ky4 = 0; ky4 < NbrCellY; ++ky4)
			    {
			      if (((ky1 + ky2 - ky3 - ky4) % NbrCellY) == 0)
				{
				  ++NbrPrecalculatedValues;
				}
			    }
		  }
	      }
      PrecalculatedValues_rhorho = new Complex [NbrPrecalculatedValues];
      PrecalculatedIndices = new int [4 * NbrPrecalculatedValues];
      NbrPrecalculatedValues = 0; 
      for (int kx1 =0; kx1 < NbrCellX; ++kx1)
	for (int kx2 =0; kx2 < NbrCellX; ++kx2)
	  for (int kx3 =0; kx3 < NbrCellX; ++kx3)
	    for (int kx4 =0; kx4 < NbrCellX; ++kx4)
	      {
		if (((kx1 + kx2 - kx3 - kx4) % NbrCellX) == 0)
		  {
		    for (int ky1 = 0; ky1 < NbrCellY; ++ky1)
		      for (int ky2 = 0; ky2 < NbrCellY; ++ky2)
			for (int ky3 = 0; ky3 < NbrCellY; ++ky3)
			  for (int ky4 = 0; ky4 < NbrCellY; ++ky4)
			    {
			      if (((ky1 + ky2 - ky3 - ky4) % NbrCellY) == 0)
				{
				  int Index1 = (kx1 * NbrCellY) + ky1;
				  int Index2 = (kx2 * NbrCellY) + ky2;
				  int Index3 = (kx3 * NbrCellY) + ky3;
				  int Index4 = (kx4 * NbrCellY) + ky4;
				  ParticleOnSphereDensityDensityOperator Operator (Space, Index1, Index2, Index3, Index4);
				  PrecalculatedValues_rhorho[NbrPrecalculatedValues] = Operator.MatrixElement(ComplexState, ComplexState);
				  PrecalculatedIndices[(NbrPrecalculatedValues << 2)] = Index1;
				  PrecalculatedIndices[(NbrPrecalculatedValues << 2) + 1] = Index2;
				  PrecalculatedIndices[(NbrPrecalculatedValues << 2) + 2] = Index3;
				  PrecalculatedIndices[(NbrPrecalculatedValues << 2) + 3] = Index4;
				  ++NbrPrecalculatedValues;
				}
			    }
		  }
	      }
    }
  else
    {
      NbrPrecalculatedValues = NbrCellX * NbrCellY;
      PrecalculatedValues_rhorho = new Complex [NbrPrecalculatedValues];
      for (int kx =0; kx < NbrCellX; ++kx)
	for (int ky = 0; ky < NbrCellY; ++ky)
	  {
	    int Index = (kx * NbrCellY) + ky;
	    ParticleOnSphereDensityOperator Operator (Space, Index);	    
	    PrecalculatedValues_rho[Index] = Operator.MatrixElement(ComplexState, ComplexState);
	  }
    }
   
  PrecalculatedValues_rho = new Complex [NbrPrecalculatedValues];
  for (int kx =0; kx < NbrCellX; ++kx)
    for (int ky = 0; ky < NbrCellY; ++ky)
      {
	int Index = (kx * NbrCellY) + ky;
	ParticleOnSphereDensityOperator Operator (Space, Index);	    
	PrecalculatedValues_rho[Index] = Operator.MatrixElement(ComplexState, ComplexState);
      }
      
  delete Space;
  ofstream File;
  File.precision(14);
  double XStep = ((double) NbrCellX) / ((double) NbrSamplesX);
  double YStep = ((double) NbrCellY) / ((double) NbrSamplesY);
  RealVector Position(2, true);
  if (Manager.GetString("output") != 0)
    File.open(Manager.GetString("output"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = 0;
      if (DensityFlag == false)
	{
	  TmpFileName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "rhorho");
	}
      else
	{
	  TmpFileName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "rho");
	}
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("state") << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }

  if (Manager.GetBoolean("k-space") == true)
    {
      if (DensityFlag == true)
	{
	  File << "# kx ky n(kx,ky)" << endl;
	  for (int kx =0; kx < NbrCellX; ++kx)
	    for (int ky = 0; ky < NbrCellY; ++ky)
	      {
		int Index = (kx * NbrCellY) + ky;
		File << kx << " " << ky << " " << PrecalculatedValues_rhorho[Index].Re << endl;
	      }
	  File.close();
	}
      return 0;
    }

  // use tight binding model to provide function basis
  TightBindingModelHofstadterSquare tightBindingModel(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, GammaX, GammaY, Architecture.GetArchitecture(), true, EmbeddingFlag);
  //
  int TotalNbrMomenta = NbrCellX * NbrCellY;
  int NbrSublattices = UnitCellX * UnitCellY;
  Complex* Coefficients = new Complex[TotalNbrMomenta];
  Complex* Coefficients2 = new Complex[TotalNbrMomenta];
  for (int i = 0; i < TotalNbrMomenta; ++i){
    int alpha=0;
    tightBindingModel.GetFunctionValue(Position, Coefficients[i], i, alpha, 0); 
  }
    
  //routine to calculate the two-particle correlation function of the form <psi|n_i n_0|psi> - <psi|n_i|psi><psi|n_0|psi>, with <psi|n_0|psi>(<psi|n_0|psi>-1) at the origin
  //
  double Normalisation = (1.0/(double)(NbrCellX*NbrCellY)); //normalisation factor for the <c^+ c> term i.e. 1/N_c
  double Normalisation2 = (1.0/(double)(NbrCellX*NbrCellY*NbrCellX*NbrCellY)); //normalisation factor for the <c^+ c^+ c c> term i.e. 1/N_c^2

  int subX, subY;
  for (int Rjx = 0; Rjx < NbrCellX; ++Rjx) //loop over MUCs
    {
      for (int Rjy = 0; Rjy < NbrCellY; ++Rjy)
	{
	  for (int alphaJ=0; alphaJ < NbrSublattices; ++alphaJ) // sublattice for r_j
	    {
              tightBindingModel.DecodeSublatticeIndex(alphaJ, subX, subY); 
	      
              Position[0] = Rjx*UnitCellX + subX; //position expressed in terms of MUC + sublattice positions
	      Position[1] = Rjy*UnitCellY + subY;
              
              Complex TmpValue = 0.0; //<psi|n_i n_0|psi>
	      Complex TmpValue_aux = 0.0; //<psi|n_i|psi><psi|n_0|psi>
	      Complex TmpValue_aux_zero = 0.0; //<psi|n_0|psi>
	      
	      if (DensityFlag == false)
		{
		  for (int i = 0; i < TotalNbrMomenta; ++i)
		    {
		      tightBindingModel.GetFunctionValue(Position, Coefficients2[i], i, alphaJ, 0);
		    }
		  for (int i = 0; i < NbrPrecalculatedValues; ++i) 
		    {
		      TmpValue += Normalisation2*(PrecalculatedValues_rhorho[i] 
						  * Conj(Coefficients[PrecalculatedIndices[(i << 2)]])
						  * Coefficients[PrecalculatedIndices[(i << 2) + 2]] 
						  * Conj(Coefficients2[PrecalculatedIndices[(i << 2) + 1]])
						  * Coefficients2[PrecalculatedIndices[(i << 2) + 3]]); //calculate <psi|n_i n_0|psi>
		      
		      TmpValue_aux += Normalisation2*(PrecalculatedValues_rho[i]*PrecalculatedValues_rho[0]
						      * Conj(Coefficients[PrecalculatedIndices[(i << 2)]])
						      * Coefficients[PrecalculatedIndices[(i << 2) + 2]] 
						      * Conj(Coefficients2[PrecalculatedIndices[(i << 2)]])
						      * Coefficients2[PrecalculatedIndices[(i << 2) + 2]]); //calculate <psi|n_i|psi><psi|n_0|psi>
		    }
		  if (Position[0]==0 && Position[1]==0)
		    {
		      TmpValue_aux_zero = Normalisation*(PrecalculatedValues_rho[0]
							 * Conj(Coefficients[PrecalculatedIndices[(0 << 2)]])
							 * Coefficients2[PrecalculatedIndices[(0 << 2) + 3]]);

		      TmpValue=TmpValue_aux_zero*(TmpValue_aux_zero - 1); //the origin is a special case, where we take <psi|n_0|psi>(<psi|n_0|psi>-1) to compensate for unwanted on-site contribution
		    }
		  else
		    {
		      TmpValue -= TmpValue_aux; //subtract <psi|n_i|psi><psi|n_0|psi>
		    }
		}
	      else
		{	      
		  // this branch is untested
		  cout << "Attention - untested branch"<<endl;
		  for (int i = 0; i < NbrPrecalculatedValues; ++i) 
		    {
		      Complex TmpValue2;
		      tightBindingModel.GetFunctionValue(Position, TmpValue2, i, alphaJ, 0);
		      TmpValue += PrecalculatedValues_rho[i] * SqrNorm(TmpValue2);
		    }
		}
              cout << Position[0] << " " << Position[1] << " " << TmpValue.Re << endl;
	      File << Position[0] << " " << Position[1] << " " << TmpValue.Re << endl;
            }
	  Position[1] += YStep;
	}
      Position[0] += XStep;
    }
  File << endl;
  File.close();
  delete[] Coefficients;
  delete[] Coefficients2;
  return 0;
}
