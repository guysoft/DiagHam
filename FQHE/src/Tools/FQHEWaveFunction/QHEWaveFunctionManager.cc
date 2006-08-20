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
#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/LaughlinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/MooreReadOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereTwoQuasiholeWaveFunction.h"
#include "Tools/FQHEWaveFunction/UnprojectedJainCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/LaughlinOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/MooreReadOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnDiskWaveFunction.h"


#include <iostream>


using std::endl;


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
  (*WaveFunctionGroup) += new SingleStringOption  ('\n', "test-wavefunction", "name of the test wave fuction");
  (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "nbr-flux", "number of quantum flux attached to each particle (laughlin and *cf only)", 1, true, 1);
  (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "cluster-size", "number of particles per cluster (read only)", 3, true, 2);
  (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "nbr-level", "number of pseudo Landau levels (filledcf only)", 1, true, 1);
  (*WaveFunctionGroup) += new SingleStringOption  ('\n', "cf-file", "name of the file describing the composite fermion state (genericcf only)");
  
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
    }
  else
    if (this->GeometryID == QHEWaveFunctionManager::DiskGeometry)
      {
	str << "  * laughlin : Laughlin wave function" << endl;
	str << "  * pfaffian : pfaffian state wave function" << endl;
	str << "  * read : Read-Rezayi state wave function" << endl;	
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
  return 0;
}
