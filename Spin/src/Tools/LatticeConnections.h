////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
//                                                                            //
//                                                                            //
//  class that defines sites and nearest neighbor relationships on a lattice  //
//                                                                            //
//                        last modification : 05/26/2009                      //
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


#ifndef LATTICECONNECTIONS_H
#define LATTICECONNECTIONS_H

#include "config.h"
#include "Options/OptionManager.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Matrix/RealMatrix.h"

#include <iostream>


using std::ostream;


class LatticeConnections
{
 protected:  
  
  // total number of sites
  int NbrSites;

  // number of unit cells
  int NbrCells;
  
  // number of sites per unit cell
  int NbrSitesPerCell;
  
  // dimension
  int Dimension;

  // Periodic repetitions in the 1,2,...,Dim direction
  int *PeriodicRep;

  // Descriptor for lattice
  char *Descriptor;

  // Matrix containing lattice vectors
  RealMatrix LatticeVectors;

  // Matrix containing coordinates of  sublattices 
  RealMatrix SubLatticeVectors;

  // neigborship relations of cells
  int **NeighborCells;
  // and number thereof
  int NbrNeighborCells;

  // compact coding of neighborship relations  
  int **Neighbors;
  // number of neighbors for spin i
  int *NbrNeighbors;

  // reduced encoding of the pairs which are interacting
  // such that Partner[i]>i
  int **Partners;
  // number of partners for spin i
  int *NbrPartners;

 public:

  // generate the object using options from Option Manager
  //
  LatticeConnections();

  // destructor
  //
  ~LatticeConnections();

  // get cell coordinates given the number of the unit cell
  // nbrCell = cell to be looked up
  // cellCoordinates = resulting coordinates, has to be reserved prior to call
  void GetCellCoordinates(int nbrCell, int *cellCoordinates);

  // get cell coordinates given the number of the unit cell
  // nbrSite = cell to be looked up
  // cellCoordinates = resulting coordinates, has to be reserved prior to call
  // sublattice = resulting sublattice
  void GetSiteCoordinates(int nbrSite, int *cellCoordinates, int &sublattice);
  
  // get number of a site in cell nbrCell
  // nbrCell = cell to be addressed
  // sublattice = sublattice index
  inline int GetSiteNumber(int nbrCell, int sublattice);

  // get number of a site in cell nbrCell
  // cellCoordinates = coordinates of cell to be addressed
  // sublattice = sublattice index
  int GetSiteNumber(int *cellCoordinates, int sublattice);

  // request total number of sites
  //
  int GetNbrSites(){return this->NbrSites;}

  // request address of partners of site
  // nbrSite = number of site whose partners to request
  // partners = array to partner sites
  // nbrPartners = number of partners found
  void GetPartners(int nbrSite, int * &partners, int &nbrPartners);

  // get a string describing the lattice geometry
  // 
  char *GeometryString();
  
  // add an option group containing all options related to the LatticeGeometry options
  //
  // manager = pointer to the option manager
  static void AddOptionGroup(OptionManager* manager);

  // pointer to the option manager
  static OptionManager* Options;

 private:
  // simple sort algorithm
  // array = integer array to be sorted
  // length = length of array
  void ArraySort(int* array, int length);

  // periodize index within fundamental interval along direction d
  // coordinate = number to periodize
  // dimension = index of dimension
  inline int Periodize(int coordinate, int dimension);
  
};


// get number of a site in cell nbrCell
// nbrCell = cell to be addressed
// sublattice = sublattice index
int LatticeConnections::GetSiteNumber(int nbrCell, int sublattice)
{
  return (nbrCell*NbrSitesPerCell+sublattice)%NbrSites;
}

// periodize index within fundamental interval along direction d
// coordinate = number to periodize
// dimension = index of dimension
int LatticeConnections::Periodize(int coordinate, int dimension)
{
  if (coordinate<0)
    coordinate += (coordinate/PeriodicRep[dimension] + 1)*PeriodicRep[dimension];
  return coordinate%PeriodicRep[dimension];
}

#endif
