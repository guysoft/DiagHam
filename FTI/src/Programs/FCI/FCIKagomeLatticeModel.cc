#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU3SpinMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandFourBodyHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandFiveBodyHamiltonian.h"

#include "Hamiltonian/ParticleOnLatticeKagomeLatticeTwoBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeKagomeLatticeThreeBandHamiltonian.h"

#include "Hamiltonian/ExplicitHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelKagomeLattice.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

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


// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, double t1, double t2, double lambda1, double lambda2);


int main(int argc, char** argv)
{
  OptionManager Manager ("FCIKagomeLatticeModel" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);  
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");

  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive next nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "five-body", "use a five body interaction instead of a two body interaction");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "real part of the nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "real part of the next nearest neighbor hopping amplitude", -0.3);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l1", "imaginary part of the nearest neighbor hopping amplitude", 0.28);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l2", "imaginary part of the next nearest neighbor hopping amplitude", 0.2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-index", "index of the band that has to be partially filled, should be 0 (lower band), 1 or 2 (upper band)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "two-bands", "use the two lowest energy bands", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-bands", "use the full three band model", 0);
  (*SystemGroup) += new BooleanOption ('\n', "project-threebands", "project the hamiltonian from the thre band model to the single band model");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice chemical potential on A site", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
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
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIKagomeLatticeModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }

  char* FilePrefix = new char [256];

  if ((Manager.GetBoolean("three-bands") == false) && (Manager.GetBoolean("two-bands") == false))
    {
      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
	{ 
	  sprintf (FilePrefix, "%s_singleband_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
	}
      else
	{
	  if (Manager.GetBoolean("three-body") == true)
	    {
	      if (Manager.GetBoolean("flat-band") == true)
		sprintf (FilePrefix, "%s_singleband_flat_threebody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
	      else
		sprintf (FilePrefix, "%s_singleband_threebody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
	    }
	  else
	    {
	      if (Manager.GetBoolean("four-body") == true)
		{
		  if (Manager.GetBoolean("flat-band") == true)
		    sprintf (FilePrefix, "%s_singleband_flat_fourbody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
		  else
		    sprintf (FilePrefix, "%s_singleband_fourbody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
		}
	      else
		{
		  if (Manager.GetBoolean("flat-band") == true)
		    sprintf (FilePrefix, "%s_singleband_flat_fivebody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
		  else
		    sprintf (FilePrefix, "%s_singleband_fivebody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
		}
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("three-bands") == false)
	{
	}
      else
	{
	  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
	    { 
	      sprintf (FilePrefix, "%s_threeband_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
	    }
	  else
	    {
	      if (Manager.GetBoolean("three-body") == true)
		{
		  if (Manager.GetBoolean("flat-band") == true)
		    sprintf (FilePrefix, "%s_threeband_flat_threebody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
		  else
		    sprintf (FilePrefix, "%s_threeband_threebody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
		}
	      else
		{
		  if (Manager.GetBoolean("four-body") == true)
		    {
		      if (Manager.GetBoolean("flat-band") == true)
			sprintf (FilePrefix, "%s_threeband_flat_fourbody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
		      else
			sprintf (FilePrefix, "%s_threeband_fourbody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
		    }
		  else
		    {
		      if (Manager.GetBoolean("flat-band") == true)
			sprintf (FilePrefix, "%s_threeband_flat_fivebody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
		      else
			sprintf (FilePrefix, "%s_threeband_fivebody_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
		    }
		}
	    }
	}
    }

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    {
      if ((Manager.GetBoolean("flat-band") == true) && (Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
	{
	  if (Manager.GetDouble("v-potential") == 0.0)
	    {
	      if (Manager.GetDouble("mu-s") == 0.0)
		sprintf (EigenvalueOutputFile, "%s_t1_%g_t2_%g_l1_%g_l2_%g_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	      else
		sprintf (EigenvalueOutputFile, "%s_t1_%g_t2_%g_l1_%g_l2_%g_gx_%g_gy_%g_mus_%g.dat",FilePrefix, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	    }
	  else
	    {
	      if (Manager.GetDouble("mu-s") == 0.0)
		sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_t1_%g_t2_%g_l1_%g_l2_%g_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	      else
		sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_t1_%g_t2_%g_l1_%g_l2_%g_gx_%g_gy_%g_mus_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	    }
	}
      else
	{
	  if (Manager.GetDouble("v-potential") == 0.0)
	    {
	      if (Manager.GetDouble("mu-s") == 0.0)
		sprintf (EigenvalueOutputFile, "%s_u_%g_t1_%g_t2_%g_l1_%g_l2_%g_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	      else
		sprintf (EigenvalueOutputFile, "%s_u_%g_t1_%g_t2_%g_l1_%g_l2_%g_gx_%g_gy_%g_mus_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	    }
	  else
	    {
	      if (Manager.GetDouble("mu-s") == 0.0)
		sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_t1_%g_t2_%g_l1_%g_l2_%g_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	      else
		sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_t1_%g_t2_%g_l1_%g_l2_%g_gx_%g_gy_%g_mus_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));

	    }
	}
    }

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
	ExportOneBody = true;
      TightBindingModelKagomeLattice TightBindingModel(NbrSitesX, NbrSitesY,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
					   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), ExportOneBody);
      if (Manager.GetBoolean("singleparticle-chernnumber") == true)
	{
	  cout << "Chern number = " << TightBindingModel.ComputeChernNumber(0) << endl;
	}
      TightBindingModel.WriteAsciiSpectrum(EigenvalueOutputFile);
      double BandSpread = TightBindingModel.ComputeBandSpread(0);
      double DirectBandGap = TightBindingModel.ComputeDirectBandGap(0);
      cout << "Spread = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	{
	  char* BandStructureOutputFile = new char [512];
	  if (Manager.GetString("export-onebodyname") != 0)
	    strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	  else
	    sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
	  if (Manager.GetBoolean("export-onebody") == true)
	    {
	      TightBindingModel.WriteBandStructure(BandStructureOutputFile);
	    }
	  else
	    {
	      TightBindingModel.WriteBandStructureASCII(BandStructureOutputFile);
	    }
	  delete[] BandStructureOutputFile;
	}	  
      return 0;

    }

  int MinKx = 0;
  int MaxKx = NbrSitesX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrSitesY - 1;
  if (Manager.GetInteger("only-ky") >= 0)
    {						
      MinKy = Manager.GetInteger("only-ky");
      MaxKy = MinKy;
    }
  
  TightBindingModelKagomeLattice TightBindingModel(NbrSitesX, NbrSitesY,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
						   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  cout << "(kx=" << i << ",ky=" << j << ") : " << endl;

	  ParticleOnSphere* Space = 0;
	  AbstractHamiltonian* Hamiltonian = 0;
	  if (Manager.GetBoolean("three-bands") == false)
	    {
	      if (Manager.GetBoolean("boson") == false)
		{
		  if ((NbrSitesX * NbrSitesY) <= 63)
		    {
		      Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		    }
		  else
		    {
		  Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		    }
		}
	      else
		{
		  Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
		{ 
		  Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), &TightBindingModel, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		}
	      else
		{ 
		  if (Manager.GetBoolean("three-body") == true)
		    {
		      // use unit three-body interaction by default
		      Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian
			(Space, NbrParticles, NbrSitesX, NbrSitesY, /* three-body */ 1.0, Manager.GetDouble("u-potential"),
			 Manager.GetDouble("v-potential"),  &TightBindingModel, 			 Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		    }
		  else
		    {
		      if (Manager.GetBoolean("four-body") == true)
			{
			  Hamiltonian = 0;
			  // 		      Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
			  // 												    Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"),
			  // 												    Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
			  // 												    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			}
		      else
			{
			  Hamiltonian = 0;
			  // 		      Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandFiveBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
			  // 												    Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"),
			  // 												    Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
			  // 												    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			}
		      
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("boson") == false)
		{
		  if ((NbrSitesX * NbrSitesY) <= 20)
		    {
		      Space = new FermionOnSquareLatticeWithSU3SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		    }
		  else
		    {
//		      Space = new FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		    }
		}
	      else
		{
		  Space = new BosonOnSquareLatticeWithSU3SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
//	      return 0;
	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
		{ 
		  Hamiltonian = new ParticleOnLatticeKagomeLatticeThreeBandHamiltonian((ParticleOnSphereWithSU3Spin*) Space, NbrParticles, NbrSitesX, NbrSitesY,
										       Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"),
										       Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
										       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		}
	      else
		{
		  Hamiltonian = 0;
		}
	      if (Manager.GetBoolean("project-threebands") == true)
		{
		  ParticleOnSphere* TargetSpace = 0;
		  ComplexMatrix* OneBodyBasis = ComputeSingleParticleTransformationMatrices(NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"));
		  ComplexMatrix SU3U1TransformationMatrix;
		  ComplexMatrix TransformedHRep;
		  ComplexMatrix HRep (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension());
 		  Hamiltonian->GetHamiltonian(HRep);
		  if (Manager.GetBoolean("boson") == false)
		    {
//		      ComplexMatrix NBodyTransformationMatrix = ((FermionOnSquareLatticeWithSU3SpinMomentumSpace*) Space)->TransformationMatrixOneBodyBasis(OneBodyBasis);
//		      TransformedHRep = HRep.Conjugate(NBodyTransformationMatrix);
//		      TargetSpace = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
//		      SU3U1TransformationMatrix = ((FermionOnSquareLatticeWithSU3SpinMomentumSpace*) Space)->TransformationMatrixSU3ToU1((FermionOnSquareLatticeMomentumSpace*) TargetSpace);
		    }
		  else
		    {
		      ComplexMatrix NBodyTransformationMatrix = ((BosonOnSquareLatticeWithSU3SpinMomentumSpace*) Space)->TransformationMatrixOneBodyBasis(OneBodyBasis);
		      TransformedHRep = HRep.Conjugate(NBodyTransformationMatrix);
		      TargetSpace = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		      SU3U1TransformationMatrix = ((BosonOnSquareLatticeWithSU3SpinMomentumSpace*) Space)->TransformationMatrixSU3ToU1((BosonOnSquareLatticeMomentumSpace*) TargetSpace, 2);
		    }

 		  ComplexMatrix TransformedHRep2 = TransformedHRep.InvConjugate(SU3U1TransformationMatrix);
		  if (Manager.GetDouble("u-potential") != 0.0)
		    TransformedHRep2 /= Manager.GetDouble("u-potential");

 		  RealDiagonalMatrix TmpDiag;
 		  HermitianMatrix HRep2(TransformedHRep2);
 		  delete Hamiltonian;
 		  delete Space;
 		  delete[] OneBodyBasis;
 		  Hamiltonian = new ExplicitHamiltonian(TargetSpace, &HRep2);
 		  Space = TargetSpace;
		}
	    }
	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d %d", i, j);
	  char* EigenstateOutputFile;
	  if (Manager.GetString("eigenstate-file") != 0)
	    {
	      EigenstateOutputFile = new char [512];
	      sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), i, j);
	    }
	  else
	    {
	      char* TmpExtention = new char [512];
	      sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
	      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	    }
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
    }
  return 0;
}


// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites

ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, double t1, double t2, double lambda1, double lambda2)
{
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[nbrSitesX * nbrSitesY];
  double KxFactor = 2.0 * M_PI / ((double) nbrSitesX);
  double KyFactor = 2.0 * M_PI / ((double) nbrSitesY);
  double NNHopping = t1;
  double NNSpinOrbit = lambda1;
  double NextNNHopping = t2;
  double NextNNSpinOrbit = lambda2;
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  int Index = (kx * nbrSitesY) + ky;
	  double KX = 0.5 * KxFactor * ((double) kx);
	  double KY = 0.5 * KyFactor * ((double) ky);
	  Complex HAB (-2.0 * NNHopping, -2.0 * NNSpinOrbit);
	  HAB *= cos (KX);
	  Complex HAC(-2.0 * NNHopping, 2.0 * NNSpinOrbit);
	  HAC *= cos (KY);
	  Complex HBC(-2.0 * NNHopping, -2.0 * NNSpinOrbit);
	  HBC *= cos(KX - KY);
	  
	  Complex HAB2 (-2.0 * NextNNHopping, 2.0 * NextNNSpinOrbit);
	  HAB2 *= cos (KX - 2.0 * KY);
	  Complex HAC2 (-2.0 * NextNNHopping, -2.0 * NextNNSpinOrbit);
	  HAC2 *= cos (2.0 * KX - KY);
	  Complex HBC2 (-2.0 * NextNNHopping, 2.0 * NextNNSpinOrbit);
	  HBC2 *= cos (KX + KY);
	  
	  HAB += HAB2;
	  HAC += HAC2;
	  HBC += HBC2;
	  
	  HermitianMatrix TmpOneBodyHamiltonian(3, true);	  
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HAB);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HAC);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HBC);
	  RealDiagonalMatrix TmpDiag;	  
	  OneBodyBasis[Index] = ComplexMatrix (3, 3);
	  OneBodyBasis[Index].SetToIdentity();
#ifdef __LAPACK__
	  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, OneBodyBasis[Index]);
#else
	  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, OneBodyBasis[Index]);
#endif   
	}
    }
  return OneBodyBasis;
}
