#include "Tools/QuantumDot/Potential/Potential.h"
#include "Tools/QuantumDot/Potential/TwoDPotential.h"
#include "Tools/QuantumDot/Potential/ThreeDPotential.h"

#include <iostream>
#include <fstream>
#include <math.h>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

double Weight_E[] = {0.128917, 0.119204, 0.103178, 0.0802035, 0.0568751, 0.0372012, 0.0217138, 0.0122897, 0.00653256, 0.00326567, 0.00153112, 0.000692242};
double Weight_H[] = {3.70785e-08, 3.38324e-07, 2.70859e-06, 1.90038e-05, 0.000116127, 0.000612045, 0.0030628, 0.0113668, 0.0345713, 0.0835946, 0.153434, 0.210873};

int main(int argc, char** argv)
{
  /*
  TwoDPotential a;
  a.ReadTwoDPotential("TestPotential1.txt", 3, 3, 2);
  a.PrintDimension(cout).flush();
  a.PrintPotential(cout).flush();
  */
  /*
  TwoDPotential a(7, 5, 3);  
  double * b= new double[7]; 
  b[0] = 0.0; b[1] = 0.3; b[2] = 0.7;  
  a.UniformWell(0.175, 0.9, b, true);
  //a.PrintDimension(cout).flush();
  a.PrintDiagram(cout).flush();
  a.PrintPotential(cout).flush();
  */
/*
  TwoDPotential c(7, 5, 12);
  ifstream fichier("Potentiel_3D.txt");
  c.ReadDiagram(fichier);
  //c.ArbitraryDistribution(0.175);
  c.UniformWell(0.175, 1.8, Weight_E, true);
  //c.PrintDiagram(cout).flush();
  c.PrintPotential(cout).flush();
  cout << c(5, 3) << '\n';
  c.SetPotentialElement(5, 3, 0.01);
  cout << c(5, 3) << '\n';
  return 0;
  */
/*
  ThreeDPotential d(3, 5, 7, 2, 2);
  d.UniformWell(0.5, -1.0, 1.0, 0.2, 2.0, true);
  d.PrintDiagram(cout).flush();
  d.PrintDimension(cout).flush();
  d.PrintPotential(cout).flush();
*/

  //ThreeDPotential(int M, int N, int H, int under, int above)
  //ThreeDPotential e(11, 11, 7, 1, 1);
  //void ArbitraryPyramidDot(double proportion, int Rb, int Rt, int w, double down_field, double wetting_field, double dot_field, double up_field, double offset, double c);
  //e.ArbitraryPyramidDot(0.9, 3, 1, 2, -0.001, 0.025, 0.025, -0.001, 1.8, 2.64, true, 0);
  /*
  e.PrintDiagram(cout).flush();
  e.PrintDimension(cout).flush(); 
  e.PrintPotential(cout).flush(); 
  */
  //bool SaveBmpPicture(int under, int above, int sizeX, int sizeY, PicRGB& InN, PicRGB& GaN, PicRGB& background, int NbrX, char* fileName);


  return 0;
}
