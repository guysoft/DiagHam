#ifndef SEGREGATIONTHREEDPOTENTIAL_H
#define SEGREGATIONTHREEDPOTENTIAL_H

#include "Potential/ThreeDPotential.h"

#include "config.h"

#include <iostream>
#include <fstream>

class SegregationThreeDPotential : public ThreeDPotential
{
 public:

  // proportion1 : probability without InN neighborhood
  // proportion2 : probability with InN neighborhood
  void SegregationWell(double down_field, double field, double up_field, double cell, double proportion1, double proportion2, double offset, bool scratch, char* fileName);

  // m, n, h : three coordinations of the considered cell
  // return true if there is any InN neighborhood
  bool Neighborhood(int m, int n, int h);

  // m, n, h : three coordinations of the considered cell
  // return true if there is any InN neighborhood
  bool NeighborhoodBis(int m, int n, int h);

 private:

};

#endif
