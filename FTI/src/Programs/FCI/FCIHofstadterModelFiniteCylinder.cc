#include "Options/Options.h"


#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"


#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"


#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterTriangularQuarter.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterFiniteCylinder.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FCIHofstadterModelFiniteCylinder" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbrsitex", "number of unit cells along the x direction", 5);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbrsitey", "number of sites in the y direction", 7);
  (*SystemGroup) += new SingleIntegerOption  ('q', "total-flux", "number of flux quanta per unit cell", 1);
  

  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive onsite(boson) or NN (fermion) potential strength", 1.0);

  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "use the real space representation when considering the system with all bandswithout the translations");
  (*SystemGroup) += new BooleanOption  ('\n', "synthetic-dimension", "use synthetic dimension coupling");

  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption  ('\n',"testhermitian-error", "precision of the hermeticity test",0);

  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSiteX = Manager.GetInteger("nbrsitex"); 
  int NbrSiteY = Manager.GetInteger("nbrsitey");


  int Flux = Manager.GetInteger("total-flux");


  char Axis ='y';

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [16];
  sprintf (StatisticPrefix, "fermions");
  
/*  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }
  */
  
  char* FilePrefix = new char [512];
  int lenFilePrefix=0;


if (Manager.GetBoolean("synthetic-dimension") == false)
      lenFilePrefix += sprintf (FilePrefix, "%s_realspace_hofstadter_q_%d", StatisticPrefix,Flux);
else
      lenFilePrefix += sprintf (FilePrefix, "%s_realspace_synth_hofstadter_q_%d", StatisticPrefix, Flux);

  // common naming options:
  lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_n_%d_x_%d_y_%d", NbrParticles, NbrSiteX, NbrSiteY);

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    {
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
      delete [] FilePrefix;
       FilePrefix = RemoveExtensionFromFileName(EigenvalueOutputFile,".dat");
      if (FilePrefix == 0)
	strcpy(FilePrefix, EigenvalueOutputFile);
    }
  else
    {
      lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_u_%g",Manager.GetDouble("u-potential"));
      sprintf (EigenvalueOutputFile,"%s.dat",FilePrefix);
    }
  
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
	ExportOneBody = true;

      Abstract2DTightBindingModel *TightBindingModel;
      TightBindingModel = new TightBindingModelHofstadterFiniteCylinder(NbrSiteX, NbrSiteY, Flux,Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);
      
      
      //TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);

    HermitianMatrix TmpHam (TightBindingModel->GetRealSpaceTightBindingHamiltonian());
    RealDiagonalMatrix TmpHam2(TmpHam.GetNbrRow());
    TmpHam.LapackDiagonalize(TmpHam2);
    for (int i = 0; i < TmpHam.GetNbrRow(); ++i)
	{
	  cout << i << " : " << TmpHam2[i] << endl;
	}

      for (int n=0; n<TightBindingModel->GetNbrBands()-1; ++n)
	{
	  double BandSpread = TightBindingModel->ComputeBandSpread(n);
	  double DirectBandGap = TightBindingModel->ComputeDirectBandGap(n);
	  cout << "Spread("<<n<<") = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;
	}
      double BandSpread = TightBindingModel->ComputeBandSpread(TightBindingModel->GetNbrBands()-1);
      cout << "Spread("<<TightBindingModel->GetNbrBands()-1<<") = " << BandSpread << endl;

      if (Manager.GetBoolean("singleparticle-chernnumber") == true)
	{
	  cout << "Chern number("<<TightBindingModel->GetNbrBands()-1<<") = " << TightBindingModel->ComputeChernNumber(TightBindingModel->GetNbrBands()-1) << endl;
	}
      if (ExportOneBody == true)
	{
	  char* BandStructureOutputFile = new char [512];
	  if (Manager.GetString("export-onebodyname") != 0)
	    strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	  else
	    sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
	  if (Manager.GetBoolean("export-onebody") == true)
	    {
	      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
	    }
	  else
	    {
	      TightBindingModel->WriteBandStructureASCII(BandStructureOutputFile);
	    }
	  delete[] BandStructureOutputFile;
	}
      delete TightBindingModel;
      return 0;
    }


  int MinKx = 0;
  int MaxKx = NbrSiteX - 1;

  if(Manager.GetBoolean("no-translation") == true)
    {  
      MaxKx = 0;
    }
  double * ChemicalPotential= new double[NbrSiteX* NbrSiteY];
  for(int i = 0 ; i <NbrSiteX* NbrSiteY ; i++)
   ChemicalPotential[i] =0.0 ;
  TightBindingModel2DAtomicLimitLattice * TightBindingModel1 = new  TightBindingModel2DAtomicLimitLattice(NbrSiteX, 1,NbrSiteY, ChemicalPotential,0,0, Architecture.GetArchitecture());

  TightBindingModelHofstadterFiniteCylinder  * TightBindingModel  = new TightBindingModelHofstadterFiniteCylinder(NbrSiteX, NbrSiteY, Flux,Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());


  HermitianMatrix TightBindingMatrix = TightBindingModel1->GetRealSpaceTightBindingHamiltonian();  
  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
	  cout << "(kx=" << i << ") : " << endl;
	  ParticleOnSphere* Space = 0;
	  AbstractQHEHamiltonian* Hamiltonian = 0;
    
        RealSymmetricMatrix DensityDensityInteraction(TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), true);
	double UPotential = Manager.GetDouble("u-potential");
        if(Manager.GetBoolean("synthetic-dimension") ==false)
        {
  		for (int x = 0; x <  NbrSiteX; ++x)
		{
		  int y = 0;
		  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x+1, y), UPotential);
         	  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y+1), UPotential);
		  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-1, y), UPotential);
		  for (y=1; y <  NbrSiteY-1; y++)
		    {
			  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x+1, y), UPotential);
			  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y+1), UPotential);
			  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-1, y), UPotential);
			  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y-1), UPotential);
	  	    }
                  y = NbrSiteY - 1;
		  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x+1, y), UPotential);
		  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-1, y), UPotential);
		  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y-1), UPotential);
		}
	}
	else
	{
		for (int x = 0; x <  NbrSiteX; ++x)
		{
		  for (int y = 0; y <  NbrSiteY; ++y)
		    {
                       for (int yint = 0; yint <  NbrSiteY; ++yint)
		    {
	        	  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x,  yint), UPotential);
	  	    }
}
		}
	
	}

cout <<"interaction Matrix " <<endl<<DensityDensityInteraction <<endl;
         if(Manager.GetBoolean("no-translation") == true)
{
   Space = new FermionOnLatticeRealSpace(NbrParticles,TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());
   cout << "Fermion in real space not yet supported"<<endl;
     exit(1);
}
   else
{
 Space = new FermionOnLatticeRealSpaceAnd1DTranslation(NbrParticles,NbrSiteX*NbrSiteY,i,NbrSiteX);

 if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
 	Memory = Architecture.GetArchitecture()->GetLocalMemory();
 Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	

 Hamiltonian = new ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(),  i,  NbrSiteX, TightBindingMatrix, DensityDensityInteraction, Architecture.GetArchitecture(), Memory);
}

  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
  double Shift = 0.0;
  Hamiltonian->ShiftHamiltonian(Shift);


 if (Manager.GetString("energy-expectation") != 0 )
	{
	  char* StateFileName = Manager.GetString("energy-expectation");
	  if (IsFile(StateFileName) == false)
	    {
	      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
	      return -1;           
	    }
	  ComplexVector State;
	  if (State.ReadVector(StateFileName) == false)
	    {
	      cout << "error while reading " << StateFileName << endl;
	      return -1;
	    }
	  if (State.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	    {
	      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
	      return -1;
	    }
	  ComplexVector TmpState(Space->GetHilbertSpaceDimension());
	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  Complex EnergyValue = State * TmpState;
	  cout << "< Energy > = " << (EnergyValue.Re - Shift) << " " << EnergyValue.Im << endl;
	  return 0; 
	}

	  
	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d", i);
	  char* EigenstateOutputFile = new char [512];

	  sprintf (EigenstateOutputFile,"%s_kx_%d",FilePrefix, i);
	    
	  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	  FirstRunFlag = false;
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  cout << "------------------------------------" << endl;
	  delete Hamiltonian;
	  delete Space;
	  delete[] EigenstateOutputFile;
	  delete[] ContentPrefix;
    }
  delete TightBindingModel;
  return 0;
}

