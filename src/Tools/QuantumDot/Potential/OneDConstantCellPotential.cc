////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//            class of potential in one direction with constant cells         //
//                                                                            //
//                        last modification : 02/11/2004                      //
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


#include "Tools/QuantumDot/Potential/OneDConstantCellPotential.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;


// destructor
//

OneDConstantCellPotential::~OneDConstantCellPotential()
{
}

// save the diagram of atoms in a file
//
// fileName = name of the file to stock the diagram

void OneDConstantCellPotential::SaveDiagram(char* fileName)
{
}

// load the diagram of atoms from a file
//
// fileName = name of the file in which the diagram is stocked

void OneDConstantCellPotential::LoadDiagram(char* fileName)
{
}

// save the potential in a file
//
// fileName = name of the file to stock the potential

void OneDConstantCellPotential::SavePotential(char* fileName)
{
  ofstream file(fileName);
  if (!file.is_open())
    cout << "Error when open the file: " << fileName << " to write potential." << endl;
  double tmp = 0.0;
  for (int i = 0; i < this->Number; ++i)
    {
      tmp = this->GetPotential(i);
      file << tmp << '\t';
    }
  file.close();
}

// load the potential from a file
//
// fileName = name of the file in which the potential is stocked

void OneDConstantCellPotential::LoadPotential(char* fileName)
{
  ifstream file(fileName);
  if (!file.is_open())
    cout << "Error when open the file: " << fileName << " to load potential." << endl;
  double tmp = 0.0;
  for (int i = 0; i < this->Number; ++i)
    {
      file >> tmp;
      this->SetPotential(i, tmp);
    }
  file.close();
}


