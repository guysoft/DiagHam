#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;

// evaluate Hilbert space dimension for fermions within a single band
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a two band
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a two band and spin conserved basis
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrSpinUp, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions on a cubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension
long FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx = 0, int currentTotalKy = 0, int currentTotalKz = 0);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETopInsulatorGetDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('z', "nbr-sitez", "number of sites along the z direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "no-inversion", "do not assume inversion symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "spin-conserved", "assume that the spin is conserved in the two band model");
  (*SystemGroup) += new BooleanOption  ('\n', "3d", "consider a 3d model instead of a 2d model");
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_topinsulator_nbrsubbands_n_nbrparticles_x_nbrsitex_y_nbrsitey.dim");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETopInsulatorGetDimension -h" << endl;
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
  int NbrSitesZ = Manager.GetInteger("nbr-sitez"); 
  
  if (Manager.GetInteger("nbr-subbands") == 1)
    {
      for (int kx = 0; kx < NbrSitesX; ++kx)
	{
	  for (int ky = 0; ky < NbrSitesY; ++ky)
	    {
	      if ((Manager.GetBoolean("no-inversion") == true) || ((kx <= ((NbrSitesX - kx) % NbrSitesX)) && (ky <= ((NbrSitesY - ky) % NbrSitesY))))
		{
		  long Dimension = FermionSingleBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
		  cout << "(kx=" << kx << ",ky=" << ky << ") : " << Dimension << endl;
		}
	    }
	}
    }
  if (Manager.GetInteger("nbr-subbands") == 2)
    {
      for (int kx = 0; kx < NbrSitesX; ++kx)
	{
	  for (int ky = 0; ky < NbrSitesY; ++ky)
	    {
	      if (Manager.GetBoolean("3d") == false)
		{
		  if ((Manager.GetBoolean("no-inversion") == true) || ((kx <= ((NbrSitesX - kx) % NbrSitesX)) && (ky <= ((NbrSitesY - ky) % NbrSitesY))))
		    {
		      if (Manager.GetBoolean("spin-conserved") == false)
			{
			  long Dimension = FermionTwoBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
			  cout << "(kx=" << kx << ",ky=" << ky << ") : " << Dimension << endl;
			}
		      else
			{
			  long TotalDimension = 0l;
			  for (int NbrSpinUp = 0; NbrSpinUp <= NbrParticles; ++NbrSpinUp)
			    {
			      long Dimension = FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1, NbrSpinUp);
			      TotalDimension += Dimension;
			      cout << "(kx=" << kx << ",ky=" << ky << ") 2Sz=" << ((2 * NbrSpinUp) - NbrParticles) << " : " << Dimension << endl;			  
			    }
			  //		      cout << "(kx=" << kx << ",ky=" << ky << ") : " << TotalDimension << endl;
			}
		    }
		}
	      else
		{
		  long TotalDimension = 0l;
		  for (int kz = 0; kz < NbrSitesZ; ++kz)
		    {
		      long Dimension = FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, kz, NbrSitesX, NbrSitesY, NbrSitesZ, NbrSitesX - 1, NbrSitesY - 1, NbrSitesZ - 1);
		      TotalDimension += Dimension;
		      cout << "(kx=" << kx << ",ky=" << ky << ",kz=" << kz << ") : " << Dimension << endl;
		    }
		  cout << TotalDimension << endl;
		}
	    }
	}
    }
      
  return 0;
}

// evaluate Hilbert space dimension for fermions within a single band
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// return value = Hilbert space dimension

long FermionSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    ++Count;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		++Count;
	    }
	}
      return Count;
    }
  Count += FermionSingleBandEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy);
  Count += FermionSingleBandEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}

// evaluate Hilbert space dimension for fermions within a two band
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// return value = Hilbert space dimension

long FermionTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count += 2l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count += 2l;
	    }
	}
      return Count;
    }
  Count += FermionTwoBandEvaluateHilbertSpaceDimension(nbrParticles - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
  Count += (2 * FermionTwoBandEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  Count += FermionTwoBandEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}

// evaluate Hilbert space dimension for fermions within a two band and spin conserved basis
//
// nbrParticles = number of nbrParticles
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int nbrSpinUp, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = nbrSiteY - 1;
      currentKx--;
    }
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrParticles))
    return 0l;

  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
	    Count++;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum))
		Count++;
	    }
	}
      return Count;
    }
  Count += FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(nbrParticles - 2, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrSpinUp - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy));
  Count += FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrSpinUp, currentTotalKx + currentKx, currentTotalKy + currentKy);
  Count += FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrSpinUp - 1, currentTotalKx + currentKx, currentTotalKy + currentKy);
  Count += FermionTwoBandWithSpinEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, nbrSiteX, nbrSiteY, currentKx, currentKy - 1, nbrSpinUp, currentTotalKx, currentTotalKy);
  return Count;
}


// evaluate Hilbert space dimension for fermions on a cubic lattice within two bands
//
// nbrParticles = number of nbrParticles
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// kxMomentum = total momentum along x
// kyMomentum = total momentum along y
// kzMomentum = total momentum along z
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// nbrSiteZ = number of sites along z
// return value = Hilbert space dimension

long FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int kxMomentum, int kyMomentum, int kzMomentum, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz)
{
  if (currentKz < 0)
    {
      currentKz = nbrSiteZ - 1;
      currentKy--;
      if (currentKy < 0)
	{
	  currentKy = nbrSiteY - 1;
	  currentKx--;
	}
    }
  if (nbrParticles == 0)
    {
      if (((currentTotalKx % nbrSiteX) == kxMomentum) && ((currentTotalKy % nbrSiteY) == kyMomentum)
	  && ((currentTotalKz % nbrSiteZ) == kzMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrParticles == 1)
    {
      for (int k = currentKz; k >= 0; --k)
	{
	  if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((currentKy + currentTotalKy) % nbrSiteY) == kyMomentum) && 
	      (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
	    Count += 2l;
	}
      for (int j = currentKy - 1; j >= 0; --j)
	{
	  for (int k = nbrSiteZ - 1; k >= 0; --k)
	    {
	      if ((((currentKx + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum)
		  && (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
		Count += 2l;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = nbrSiteY - 1; j >= 0; --j)
	    {
	      for (int k = nbrSiteZ - 1; k >= 0; --k)
		{
		  if ((((i + currentTotalKx) % nbrSiteX) == kxMomentum) && (((j + currentTotalKy) % nbrSiteY) == kyMomentum)
		      && (((k + currentTotalKz) % nbrSiteZ) == kzMomentum))
		    Count += 2l;
		}
	    }
	}
      return Count;
    }
  Count += FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles - 2, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz));
  Count += (2 * FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles - 1, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz));
  Count += FermionCubicLatticeTwoBandEvaluateHilbertSpaceDimension(nbrParticles, kxMomentum, kyMomentum, kzMomentum, nbrSiteX, nbrSiteY, nbrSiteZ, currentKx, currentKy, currentKz - 1, currentTotalKx, currentTotalKy, currentTotalKz);
  return Count;
}
