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
// KxMomentum = total momentum along x
// KyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionSingleBandEvaluateHilbertSpaceDimension(int nbrParticles, int KxMomentum, int KyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);

// evaluate Hilbert space dimension for fermions within a two band
//
// nbrParticles = number of nbrParticles
// KxMomentum = total momentum along x
// KyMomentum = total momentum along y
// nbrSiteX = number of sites along x
// nbrSiteY = number of sites along y
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension
long FermionTwoBandEvaluateHilbertSpaceDimension(int nbrParticles, int KxMomentum, int KyMomentum, int nbrSiteX, int nbrSiteY, int currentKx, int currentKy, int currentTotalKx = 0, int currentTotalKy = 0);


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
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
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
  
  if (Manager.GetInteger("nbr-subbands") == 1)
    {
      for (int kx = 0; kx < NbrSitesX; ++kx)
	{
	  for (int ky = 0; ky < NbrSitesY; ++ky)
	    {
	      long Dimension = FermionSingleBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
	      cout << "(kx=" << kx << ",ky=" << ky << ") : " << Dimension << endl;
	    }
	}
    }
  if (Manager.GetInteger("nbr-subbands") == 2)
    {
      for (int kx = 0; kx < NbrSitesX; ++kx)
	{
	  for (int ky = 0; ky < NbrSitesY; ++ky)
	    {
	      long Dimension = FermionTwoBandEvaluateHilbertSpaceDimension(NbrParticles, kx, ky, NbrSitesX, NbrSitesY, NbrSitesX - 1, NbrSitesY - 1);
	      cout << "(kx=" << kx << ",ky=" << ky << ") : " << Dimension << endl;
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
// KxMomentum = total momentum along x
// KyMomentum = total momentum along y
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
// KxMomentum = total momentum along x
// KyMomentum = total momentum along y
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
