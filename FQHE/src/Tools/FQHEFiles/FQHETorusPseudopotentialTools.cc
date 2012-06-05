////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          set of functions used to handle torus pseudopotential files       //
//                                                                            //
//                        last modification : 17/04/2012                      //
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
#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"
#include "GeneralTools/ConfigurationParser.h"


#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>

using std::cout;
using std::endl;
using std::string;


// get pseudopototentials for particles on torus from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = reference on the number of pseudopotentials
// pseudoPotentials = reference on the array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// return value = true if no error occured

bool FQHETorusGetPseudopotentials (char* fileName, int& nbrPseudoPotentials, double*& pseudoPotentials)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', pseudoPotentials, nbrPseudoPotentials) == false)
    {
      cout << "Pseudopotentials has a wrong value in " << fileName << endl;
      return false;
    }
  return true;
}


// get pseudopototentials for particles on torus from file
// 
// fileName = name of the file that contains the pseudopotantial description
// landauLevel = index of Coulomb Landau-level
// nbrPseudoPotentials = reference on the number of pseudopotentials
// pseudoPotentials = reference on the array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
// interactionName = naming convention read from definition, or generated from LL index
// return value = true if no error occured

bool FQHETorusGetPseudopotentials (char* fileName, bool haveCoulomb, int &landauLevel, int& nbrPseudoPotentials, double*& pseudoPotentials, char*& interactionName)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }

  if (InteractionDefinition["CoulombLandauLevel"] != NULL)
    {
      landauLevel = atoi(InteractionDefinition["CoulombLandauLevel"]);
      haveCoulomb=true;
    }
  else
    {
      haveCoulomb=false;
    }
  if (InteractionDefinition["Name"] == NULL)
    {
      if ((InteractionDefinition["CoulombLandauLevel"] != NULL) && (InteractionDefinition["Pseudopotentials"] == NULL))
	{
	  interactionName = new char[18];
	  if (landauLevel>=0)
	    sprintf(interactionName,"coulomb_l_%d",landauLevel);
	  else
	    sprintf(interactionName,"graphene_l_%d",-landauLevel);
	}
      else
	{
	  cout << "Attention, using unnamed interaction! Please include a line 'Name = ...'" << endl;
	  interactionName = new char[10];
	  sprintf(interactionName,"unnamed");
	}
    }
  else
    {
      interactionName = new char[strlen(InteractionDefinition["Name"])+1];
      strcpy(interactionName, InteractionDefinition["Name"]);
    }
  InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', pseudoPotentials, nbrPseudoPotentials);
  
  return true;
}


// get pseudopototentials for particles on torus with SU(2) spin from file
// 
// fileName = name of the file that contains the pseudopotantial description
// nbrPseudoPotentials = number of pseudopotentials per interaction type
// pseudoPotentials = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// return value = true if no error occured

bool FQHETorusSU2GetPseudopotentials (char* fileName, int* nbrPseudoPotentials, double** pseudoPotentials)
{
  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(fileName) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      return false;
    }
  int TmpNbrPseudoPotentials;
  double* TmpPseudoPotentials;
  bool Flag = false;
  if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      for (int i = 0; i < 3; ++i)
	{
	  nbrPseudoPotentials[i] = TmpNbrPseudoPotentials;
	  pseudoPotentials[i] = new double[nbrPseudoPotentials[i]];
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    pseudoPotentials[i][j] = TmpPseudoPotentials[j];
	}
    }
  else
    if (InteractionDefinition["Pseudopotentials"] != 0)
      {
	cout << "Pseudopotentials has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpUp", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[0] = TmpNbrPseudoPotentials;
      pseudoPotentials[0] = new double[nbrPseudoPotentials[0]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[0][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpUp"] != 0)
      {
	cout << "PseudopotentialsUpUp has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[1] = TmpNbrPseudoPotentials;
      pseudoPotentials[1] = new double[nbrPseudoPotentials[1]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[1][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsDownDown"] != 0)
      {
	cout << "PseudopotentialsDownDown has a wrong value in " << fileName << endl;
	return false;
      }
  if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
    {
      Flag = true;
      nbrPseudoPotentials[2] = TmpNbrPseudoPotentials;
      pseudoPotentials[2] = new double[nbrPseudoPotentials[2]];
      for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	pseudoPotentials[2][j] = TmpPseudoPotentials[j];
    }
  else
    if (InteractionDefinition["PseudopotentialsUpDown"] != 0)
      {
	cout << "PseudopotentialsUpDown has a wrong value in " << fileName << endl;
	return false;
      }
  return true;
}


