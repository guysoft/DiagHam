////////////////////////////////////////////////////////////////////////////////
//                                                                            // 
//                     Copyright (C) 2004 Duc-Phuong Nguyen                   //
//                                                                            //
//    class of dot potential in three directions with constant cylinders      //
//                                                                            //
//                        last modification : 04/22/2004                      //
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


#include "Tools/QuantumDot/Potential/QuantumDotThreeDConstantCylinderPotential.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;


// constructor from data
//
// belowHeight = height of cylinder below the wetting layer
// wettingWidth = width of the wetting layer
// nbrCylinderDot = number of cylinders in the dot
// dotHeight = height of the dot
// baseRadius = radius of the dot base
// topRadius = radius of the dot top
// aboveHeight = height of cylinder above the dot

QuantumDotThreeDConstantCylinderPotential::QuantumDotThreeDConstantCylinderPotential(double belowHeight, double wettingWidth, int nbrCylinderDot, double dotHeight, double baseRadius, double topRadius, double aboveHeight)
{
  this->NbrCylinderDot = nbrCylinderDot;
  this->NumberZ = this->NbrCylinderDot + 3;
  this->CylinderHeight = new double[this->NumberZ];
  this->CylinderRadius = new double[this->NumberZ];
  this->PotentialValue =  new double[this->NumberZ];

  this->BelowHeight = belowHeight;
  this->CylinderHeight[0] = this->BelowHeight;
  this->CylinderRadius[0] = -1.0;

  this->WettingWidth = wettingWidth;
  this->CylinderHeight[1] = this->WettingWidth;
  this->CylinderRadius[1] = -1.0;

  this->DotHeight = dotHeight;
  this->BaseRadius = baseRadius;
  this->TopRadius = topRadius;
  for (int k = 2; k < (2 + this->NbrCylinderDot); ++k)
    {
      this->CylinderHeight[k] = this->DotHeight / double(this->NbrCylinderDot);
      this->CylinderRadius[k] = this->BaseRadius - double(k - 2) * (this->BaseRadius - this->TopRadius) / double(this->NbrCylinderDot - 1);
    }
  
  this->AboveHeight = aboveHeight;
  this->CylinderHeight[this->NumberZ - 1] = this->AboveHeight;
  this->CylinderRadius[this->NumberZ - 1] = -1.0;

}

// destructor
//

QuantumDotThreeDConstantCylinderPotential::~QuantumDotThreeDConstantCylinderPotential()
{
}

// construct the potential profile from potential values
//

void QuantumDotThreeDConstantCylinderPotential::ConstructPotential(double dotPotential)
{
  this->PotentialValue[0] = 0.0;
  this->PotentialValue[1] = dotPotential;
  for (int k = 0; k < this->NbrCylinderDot; ++k)
    this->PotentialValue[k + 2] = dotPotential;
  this->PotentialValue[this->NumberZ] = 0.0;
}
