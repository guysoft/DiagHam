#include "Tools/QuantumDot/Potential/AbstractPotential.h"
#include "Tools/QuantumDot/Potential/OneDConstantCellPotential.h"
#include "Tools/QuantumDot/Potential/TwoDConstantCellPotential.h"
#include "Tools/QuantumDot/Potential/ThreeDConstantCellPotential.h"
#include "Tools/QuantumDot/Potential/HardBoxPyramidQuantumDotThreeDConstantCellPotential.h"

#include <iostream>
#include <fstream>
#include <math.h>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

int main(int argc, char** argv)
{  

  HardBoxPyramidQuantumDotThreeDConstantCellPotential* potential = new HardBoxPyramidQuantumDotThreeDConstantCellPotential(50, 50, 30, 6, 10, 4, 20, 5);  
  potential->LoadDiagram("Diagram.txt");
  potential->ConstructPotential(0.1, 0.295, -0.001, 0.0245, 0.0245, -0.001, 0.9, 2.64, false, "Parameter.txt");
  potential->SavePotentialWithConstantField("PyramidDotPotential.txt"); 
  //Potential->SaveDiagram("Diagrambis.txt");


  return 0;
}
