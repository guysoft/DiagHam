////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of FQHE wave function manager                    //
//                                                                            //
//                        last modification : 18/01/2005                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "Options/Options.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/ExtendedHalperinWavefunction.h"
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/LaughlinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/MooreReadOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereTwoQuasiholeWaveFunction.h"
#include "Tools/FQHEWaveFunction/UnprojectedJainCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PairedCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/LaughlinOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/MooreReadOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/PairedCFOnSphereWithSpinWaveFunction.h"
#include "Tools/FQHEWaveFunction/TwoThirdSingletState.h"
#include "Tools/FQHEWaveFunction/HundRuleCFStates.h"
#include "Tools/FQHEWaveFunction/HundRuleBilayerSinglet.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/NumericalAnalysis/AntisymmetrizedComplexFunction.h"

#include <iostream>


using std::endl;
using std::cout;

// constructor
//
// geometry = id of the geometry to use

QHEWaveFunctionManager::QHEWaveFunctionManager(int geometry)
{
  this->GeometryID = geometry;
  this->Options = 0;
}

// destructor
//

QHEWaveFunctionManager::~QHEWaveFunctionManager()
{
}

// add an option group containing all options related to the wave functions
//
// manager = pointer to the option manager

void QHEWaveFunctionManager::AddOptionGroup(OptionManager* manager)
{
  this->Options = manager;
  OptionGroup* WaveFunctionGroup  = new OptionGroup ("analytical wave function options");
  (*(this->Options)) += WaveFunctionGroup;
  (*WaveFunctionGroup) += new SingleStringOption  ('\n', "test-wavefunction", "name of the test wave fuction",0);
  (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "nbr-flux", "number of quantum flux attached to each particle (laughlin and *cf only)", 1, true, 1);
  (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "cluster-size", "number of particles per cluster (read only)", 3, true, 2);
  (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "nbr-level", "number of pseudo Landau levels (filledcf only)", 1, true, 1);
  (*WaveFunctionGroup) += new SingleStringOption  ('\n', "cf-file", "name of the file describing the composite fermion state (genericcf only)");
  if (this->GeometryID & QHEWaveFunctionManager::SphereGeometry)
    {
      // PairedCFOptions:
      (*WaveFunctionGroup) += new SingleDoubleOption  ('\n', "MR-coeff", "coefficient for Moore-Read contribution (pairedcf only)",1.0);
      (*WaveFunctionGroup) += new MultipleDoubleOption  ('\n', "pair-coeff", "sequence of pairing coefficients (pairedcf only)",'+');
      (*WaveFunctionGroup) += new BooleanOption  ('\n', "pair-compatibility", "adopt old conventions for normalisation (pairedcf only)");
      (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "Jz-Value", "Total angular momentum Jz (hund only)", 0);
    }
  else if (this->GeometryID & QHEWaveFunctionManager::SphereWithSpinGeometry)
    {
      // PairedCF(CB)Options:
      (*WaveFunctionGroup) += new SingleDoubleOption  ('\n', "bosons", "coefficient for boson contribution (pairedcfcb only)",1.0);
      (*WaveFunctionGroup) += new MultipleDoubleOption  ('\n', "pair-coeff", "sequence of pairing coefficients (pairedcf* only)",'+');
      (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "pair-wave", "choose pairing channel s,+p,... (pairedcf* only)",1);
      (*WaveFunctionGroup) += new BooleanOption  ('\n', "pair-compatibility", "adopt old conventions for normalisation (pairedcf* only)");
      (*WaveFunctionGroup) += new BooleanOption  ('\n', "jastrow-inside", "move jastrow factors inside determinant (1s only)");
      (*WaveFunctionGroup) += new MultipleIntegerOption  ('\n', "XHC", "coefficients (k,l,q,r) of extended Halperin wavefunction  (XH only)",',' ,',', "2,2,-2,0");
      (*WaveFunctionGroup) += new BooleanOption  ('\n', "antisymmetrize", "antisymmetrize wavefunction (XH only)");
    }
}

// get list of all available wave functions
// 
// str = reference on the output stream

ostream& QHEWaveFunctionManager::ShowAvalaibleWaveFunctions (ostream& str)
{
  str << "list of avalaible wave functions:" << endl;
  if (this->GeometryID == QHEWaveFunctionManager::SphereGeometry)
    {
      str << "  * laughlin : Laughlin wave function" << endl;
      str << "  * pfaffian : pfaffian state wave function" << endl;
      str << "  * pfaffian2qh : pfaffian state wave function with 2 quasiholes" << endl;
      str << "  * read : Read-Rezayi state wave function" << endl;
      str << "  * filledcf : composite fermions wave function (only with filled pseudo Landau levels)" << endl;
      str << "  * genericcf : generic composite fermions wave function" << endl;            
      str << "  * unprojectedcf : generic unprojected composite fermions wave function" << endl;
      str << "  * pairedcf : paired composite fermion wave function at flux 2N-3" << endl;
      str << "  * hund : composite fermion state, if half filled highest shell using Hund's rule" << endl;
    }
  else
    if (this->GeometryID == QHEWaveFunctionManager::DiskGeometry)
      {
	str << "  * laughlin : Laughlin wave function" << endl;
	str << "  * pfaffian : pfaffian state wave function" << endl;
	str << "  * read : Read-Rezayi state wave function" << endl;	
      }
    else
      if (this->GeometryID == QHEWaveFunctionManager::SphereWithSpinGeometry)
	{	  
	  str << "  * 111 : 111-state" << endl;
	  str << "  * HR : Haldane-Rezayi state" << endl;
	  str << "  * XH : Extended Halperin wavefunctions (z-z)^k (w-w)^k (z-w)^l det((z-w)^q) per((z-w)^r)" << endl;
	  str << "  * 1s : Spin-singlet state at filling one" << endl;
	  str << "  * 2-3s : Spin-singlet state at filling two-thirds" << endl;
	  str << "  * pairedcf : paired composite fermion wave function at flux 2N_1-1" << endl;	  
	  str << "  * pairedcfcb : paired composite fermion wave function at flux 2N_1-1 with CB component" << endl;
	  str << "  * hund : singlet state with each layer formed according to Hund's rule" << endl;
	}
  return str;
}
  
// get the wave function corresponding to the option constraints
//
// return value = pointer to the wave function (null if an error occurs)

Abstract1DComplexFunction* QHEWaveFunctionManager::GetWaveFunction()
{
  if ((*(this->Options))["test-wavefunction"] == 0)
    {
      return 0;
    }
  if (this->Options->GetString("test-wavefunction") == 0)
    {
      return 0;
    }
  if (this->GeometryID == QHEWaveFunctionManager::SphereGeometry)
    {
      if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "laughlin") == 0)
	{
	  return new LaughlinOnSphereWaveFunction(((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger(), 
						  ((SingleIntegerOption*) (*(this->Options))["nbr-flux"])->GetInteger() + 1);
	}
      if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pfaffian") == 0)
	{
	  return new PfaffianOnSphereWaveFunction(((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger());
	}
      if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pfaffian2qh") == 0)
	{
	  return new PfaffianOnSphereTwoQuasiholeWaveFunction(((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger(), 0.0, 0.0, M_PI, 0.0);
	}
      if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "read") == 0)
	{
	  return new MooreReadOnSphereWaveFunction(((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger(), 
						   ((SingleIntegerOption*) (*(this->Options))["cluster-size"])->GetInteger());
	}
      if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "filledcf") == 0)
	{
	  return new JainCFFilledLevelOnSphereWaveFunction(((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger(), 
							   ((SingleIntegerOption*) (*(this->Options))["nbr-level"])->GetInteger(),
							   ((SingleIntegerOption*) (*(this->Options))["nbr-flux"])->GetInteger());
	}
      if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "genericcf") == 0) && ((*(this->Options))["cf-file"] != 0))
	{	  
	  return new JainCFOnSphereWaveFunction(((SingleStringOption*) (*(this->Options))["cf-file"])->GetString());
	}
      if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "unprojectedcf") == 0) && ((*(this->Options))["cf-file"] != 0))
	{	  
	  return new UnprojectedJainCFOnSphereWaveFunction(((SingleStringOption*) (*(this->Options))["cf-file"])->GetString());
	}
      if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pairedcf") == 0))
	{
	  int N= ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
	  double *Coefficients = ((MultipleDoubleOption*) (*(this->Options))["pair-coeff"])->GetDoubles();
	  int LL;
	  if (Coefficients==NULL)
	    {
	      Coefficients = new double[1];
	      Coefficients[0]=0.0;
	      LL=1;
	    }
	  else
	    LL = ((MultipleDoubleOption*) (*(this->Options))["pair-coeff"])->GetLength();
	  double MR =((SingleDoubleOption*) (*(this->Options))["MR-coeff"])->GetDouble();
	  bool conventions = ((BooleanOption*) (*(this->Options))["pair-compatibility"])->GetBoolean();
	  PairedCFOnSphereWaveFunction* rst = new PairedCFOnSphereWaveFunction(N, LL, -1, MR, Coefficients, conventions, 2);
	  rst->AdaptAverageMCNorm();
	  delete [] Coefficients;
	  return rst;
	}
      if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "hund") == 0))
	{
	  int N = this->Options->GetInteger("nbr-particles");
	  int JastrowP = this->Options->GetInteger("nbr-flux")/2;
	  int effectiveFlux = this->Options->GetInteger("lzmax")-2*JastrowP*(N-1);
	  if (JastrowP==0)
	    {
	      cout << "To obtain CF's, at least two flux need to be attached. Try:  --nbr-flux 2"<<endl;
	      exit(1);
	    }
	  cout << "Effective flux for CF's: "<<effectiveFlux<<endl;
	  HundRuleCFStates* rst = new HundRuleCFStates(N, effectiveFlux, JastrowP);
	  rst->SelectMValue(this->Options->GetInteger("Jz-Value"));
	  rst->AdaptAverageMCNorm();
	  return rst;
	}
      return 0;
    }
  else
    if (this->GeometryID == QHEWaveFunctionManager::DiskGeometry)
      {
	if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "laughlin") == 0)
	  {
	    return new LaughlinOnDiskWaveFunction(((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger(), 
						  ((SingleIntegerOption*) (*(this->Options))["nbr-flux"])->GetInteger() + 1);
	  }
	if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pfaffian") == 0)
	  {
	    return new PfaffianOnDiskWaveFunction(((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger());
	  }
	if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "read") == 0)
	  {
	    return new MooreReadOnDiskWaveFunction(((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger(), 
						   ((SingleIntegerOption*) (*(this->Options))["cluster-size"])->GetInteger());
	  }
	return 0;
      }
    else
      if (this->GeometryID == QHEWaveFunctionManager::SphereWithSpinGeometry)
	{
	  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pairedcf") == 0))
	    {
	      int N= ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
	      int Sz= Options->GetInteger("SzTotal");
	      if ((N&1)||(Sz!=0))
		{
		  cout << "Paired CF bilayer states require equal population of both (pseudo-)spin species and even N!" << endl;
		  exit(1);
		}
	      double *Coefficients = ((MultipleDoubleOption*) (*(this->Options))["pair-coeff"])->GetDoubles();
	      int LL;
	      if (Coefficients==NULL)
		{
		  Coefficients = new double[1];
		  Coefficients[0]=0.0;
		  LL=1;
		}
	      else
		LL = ((MultipleDoubleOption*) (*(this->Options))["pair-coeff"])->GetLength();
	      bool conventions = ((BooleanOption*) (*(this->Options))["pair-compatibility"])->GetBoolean();
	      int pairWave = this->Options->GetInteger("pair-wave");
	      PairedCFOnSphereWithSpinWaveFunction* rst = new PairedCFOnSphereWithSpinWaveFunction(N, LL, pairWave, false, 0.0, Coefficients, conventions, 2);
	      rst->AdaptAverageMCNorm();
	      delete [] Coefficients;
	      return rst;
	    }
	  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pairedcfcb") == 0))
	    {
	      int N= ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
	      int Sz= Options->GetInteger("SzTotal");
	      if ((N&1)||(Sz!=0))
		{
		  cout << "Paired CF-CB bilayer states require equal population of both (pseudo-)spin species and even N!" << endl;
		  exit(1);
		}
	      double *Coefficients = ((MultipleDoubleOption*) (*(this->Options))["pair-coeff"])->GetDoubles();
	      int LL;
	      if (Coefficients==NULL)
		{
		  Coefficients = new double[1];
		  Coefficients[0]=0.0;
		  LL=1;
		}
	      else
		LL = ((MultipleDoubleOption*) (*(this->Options))["pair-coeff"])->GetLength();
	      double BC = ((SingleDoubleOption*) (*(this->Options))["bosons"])->GetDouble();
	      bool conventions = ((BooleanOption*) (*(this->Options))["pair-compatibility"])->GetBoolean();
	      PairedCFOnSphereWithSpinWaveFunction* rst = new PairedCFOnSphereWithSpinWaveFunction(N, LL, 1, true, BC, Coefficients, conventions, 2);
	      rst->AdaptAverageMCNorm();
	      delete [] Coefficients;
	      return rst;
	    }
	  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "hund") == 0))
	    {
	      int N= Options->GetInteger("nbr-particles");
	      int n=N/2;
	      HundRuleBilayerSinglet* rst = new HundRuleBilayerSinglet(n);
	      rst->AdaptAverageMCNorm();
	      return rst;	      
	    }
	  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "111-old") == 0))
	    {
	      int N= ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
	      int Sz= Options->GetInteger("SzTotal");
	      if ((N&1)||(Sz!=0))
		{
		  cout << "For now, the implementation of the 111 state requires Sz=0 and even N!" << endl;
		  exit(1);
		}
	      double* Coefficients = new double[1];
	      Coefficients[0]=0.0;
	      bool conventions = ((BooleanOption*) (*(this->Options))["pair-compatibility"])->GetBoolean();
	      PairedCFOnSphereWithSpinWaveFunction* rst = new PairedCFOnSphereWithSpinWaveFunction(N, 1, 1, true, 1.0, Coefficients, conventions, 2);
	      rst->AdaptAverageMCNorm();
	      delete[] Coefficients;
	      return rst;
	    }
	  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "111") == 0))
	    {
	      int N= ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
	      ExtendedHalperinWavefunction* rst = new ExtendedHalperinWavefunction(N, 1, 1, 0);
	      rst->AdaptAverageMCNorm();
	      return rst;
	    }
	  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "HR") == 0))
	    {
	      int N= ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
	      ExtendedHalperinWavefunction* rst = new ExtendedHalperinWavefunction(N, 2, 2, -2);
	      rst->AdaptAverageMCNorm();
	      return rst;
	    }
	  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "XH") == 0))
	    {
	      int N= ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
	      int length;
	      int *Params = Options->GetIntegers("XHC",length);
	      int K=1, M=1, P=0, Q=0, R=1, S=1;
	      if (length>0) K = Params[0];
	      if (length>1) M = Params[1];
	      if (length>2) P = Params[2];
	      if (length>3) Q = Params[3];
	      if (length>4) R = Params[4];
	      if (length>5) S = Params[5];
	      ExtendedHalperinWavefunction* rst = new ExtendedHalperinWavefunction(N, K, M, P, Q, R, S);
	      rst->AdaptAverageMCNorm();
	      if (Options->GetBoolean("antisymmetrize"))
		{
		  AntisymmetrizedComplexFunction* rst2 = new AntisymmetrizedComplexFunction(rst,N,2);		  
		  return rst2;
		}
	      else
		return rst;
	    }
	  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "1s") == 0))
	    {
	      int N= ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
	      bool inside = Options->GetBoolean("jastrow-inside");
	      ExtendedHalperinWavefunction* rst = new ExtendedHalperinWavefunction(N, 0, 2, -2, 0, 1, 1, inside);
	      rst->AdaptAverageMCNorm();
	      return rst;
	    }
	  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "2-3s") == 0))
	    {
	      int N= ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();	      
	      TwoThirdSingletState* rst = new TwoThirdSingletState(N);
	      rst->AdaptAverageMCNorm();
	      return rst;
	    }
	}
  return 0;
}

char* QHEWaveFunctionManager::GetDescription()
{
  if ((*(this->Options))["test-wavefunction"] == 0)
    {
      return 0;
    }
  char * buffer = new char[1000];
  sprintf(buffer,"%s N=%d",this->Options->GetString("test-wavefunction"), this->Options->GetInteger("nbr-particles"));
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pairedcf") == 0))
    {
      double *Coefficients = ((MultipleDoubleOption*) (*(this->Options))["pair-coeff"])->GetDoubles();
      if (Coefficients==NULL)
	{
	  if (this->GeometryID & QHEWaveFunctionManager::SphereWithSpinGeometry)
	    sprintf(buffer,"%s, B: %g, c: 0",buffer, this->Options->GetDouble("bosons"));
	  else
	    sprintf(buffer,"%s, MR: %g, c: 0",buffer, this->Options->GetDouble("MR-coeff"));
	}
      else
	{
	  if (this->GeometryID & QHEWaveFunctionManager::SphereWithSpinGeometry)
	    sprintf(buffer,"%s, B: %g, c: %g",buffer, this->Options->GetDouble("bosons"), Coefficients[0]);
	  else
	    sprintf(buffer,"%s, MR: %g, c: %g",buffer, this->Options->GetDouble("MR-coeff"), Coefficients[0]);
	  int LL = ((MultipleDoubleOption*) (*(this->Options))["pair-coeff"])->GetLength();
	  for (int i=1; i<LL; ++i)
	    sprintf(buffer,"%s+%g",buffer, Coefficients[i]);
	  if(this->Options->GetBoolean("pair-compatibility"))
	    sprintf(buffer,"%s (c)",buffer);
	}
    }
  char *rst = new char[strlen(buffer)+1];
  strcpy(rst,buffer);
  delete [] buffer;
  return rst;
}


int QHEWaveFunctionManager::GetWaveFunctionType()
{
  if ((*(this->Options))["test-wavefunction"] == 0)
    return QHEWaveFunctionManager::InvalidWaveFunction;
  if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "laughlin") == 0)
    return QHEWaveFunctionManager::Laughlin;
  if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pfaffian") == 0)
    return QHEWaveFunctionManager::Pfaffian;
  if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pfaffian2qh") == 0)
    return QHEWaveFunctionManager::Pfaffian2QH;
  if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "read") == 0)
    return QHEWaveFunctionManager::ReadRezayi;
  if (strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "filledcf") == 0)
    return QHEWaveFunctionManager::FilledCF;
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "genericcf") == 0) && ((*(this->Options))["cf-file"] != 0))
    return QHEWaveFunctionManager::GenericCF;
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "unprojectedcf") == 0) && ((*(this->Options))["cf-file"] != 0))
    return QHEWaveFunctionManager::UnprojectedCF;
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pairedcf") == 0))
    return QHEWaveFunctionManager::PairedCF;  
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "pairedcfcb") == 0))
    return QHEWaveFunctionManager::PairedCFCB;
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "111") == 0))
    return QHEWaveFunctionManager::OneOneOne;
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "HR") == 0))
    return QHEWaveFunctionManager::HaldaneRezayi;
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "XH") == 0))
    return QHEWaveFunctionManager::ExtendedHalperin;
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "1s") == 0))
    return QHEWaveFunctionManager::OneS;
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "2-3s") == 0))
    return QHEWaveFunctionManager::TwoThirdsS;
  if ((strcmp (((SingleStringOption*) (*(this->Options))["test-wavefunction"])->GetString(), "hund") == 0))
    return QHEWaveFunctionManager::HundRuleSinglet;
  return QHEWaveFunctionManager::InvalidWaveFunction;
}
