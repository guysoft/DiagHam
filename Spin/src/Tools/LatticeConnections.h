////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar M�ller                        //
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
  // array indicating shift of sites around periodic boundaries in the different dimensions (site,neighborIdx,dimension)
  int ***NeighborShift;

  
  // reduced encoding of the pairs which are interacting
  // such that Partner[i]>i
  int **Partners;
  // number of partners for spin i
  int *NbrPartners;

  // plaquettes of the lattice, with indices going counterclockwise
  // number of plaquettes
  int NbrPlaquettes;
  // number of spins on plaquette
  int *NbrPlaquetteSpins;
  // spins involved in each plaquette
  int **PlaquetteSpins;
    
 public:

  // generate the object using options from Option Manager
  //
  LatticeConnections();

  // destructor
  //
  ~LatticeConnections();

  // request address of partners of site
  // nbrSite = number of site whose partners to request
  // nbrNeighbors = number of partners found
  // Neighbors = array to partner sites
  // periodicTranslations = translations into the fundamental domain
  void GetNeighbors(int nbrSite, int &nbrNeighbors, int * &neighbors, int **&periodicTranslations);

  // get cell coordinates given the number of the unit cell
  // nbrCell = cell to be looked up
  // cellCoordinates = resulting coordinates, has to be reserved prior to call
  void GetCellCoordinates(int nbrCell, int *cellCoordinates);

  // get cell coordinates given the number of the unit cell
  // nbrSite = cell to be looked up
  // cellCoordinates = resulting coordinates, has to be reserved prior to call
  // sublattice = resulting sublattice
  void GetSiteCoordinates(int nbrSite, int *cellCoordinates, int &sublattice);

  // retrieve the position of a given site
  // cellCoordinates = resulting coordinates, has to be reserved prior to call
  // sublattice = resulting sublattice  
  RealVector GetSitePosition(int *cellCoordinates, int sublattice);
  
  // get number of a site in cell nbrCell
  // nbrCell = cell to be addressed
  // sublattice = sublattice index
  inline int GetSiteNumber(int nbrCell, int sublattice);

  // get number of a site in cell nbrCell
  // cellCoordinates = coordinates of cell to be addressed
  // sublattice = sublattice index
  int GetSiteNumber(int *cellCoordinates, int sublattice);

  // get number of a site in cell nbrCell, and return translation vector back into the simulation cell
  // cellCoordinates = coordinates of cell to be addressed
  // sublattice = sublattice index
  // translation = vector of tranlation back into simulation cell
  int GetSiteNumber(int *cellCoordinates, int sublattice, int *translation);

  // request total number of sites
  //
  int GetNbrSites(){return this->NbrSites;}

  // access lattice extension in d-th direction
  int GetNbrSubLattices(){return this->NbrSitesPerCell;}

  // access lattice extension in d-th direction
  int GetLatticeLength(int direction){return this->PeriodicRep[direction];}

  // access dimension
  int GetLatticeDimension(){return this->Dimension;}

  // request address of partners of site
  // nbrSite = number of site whose partners to request
  // partners = array to partner sites
  // nbrPartners = number of partners found
  void GetPartners(int nbrSite, int * &partners, int &nbrPartners);

  // get number of plaquettes
  int GetNbrPlaquettes(){return this->NbrPlaquettes;}

  // get spins of a given plaquette
  void GetPlaquetteSpins(int nbrPlaquette, int * &spins, int &nbrSpins);
  
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

  // simple sort algorithm
  // array = integer array to be sorted
  // length = length of array
  void ArraySort2(int* array, int length, int** array2);

  // periodize index within fundamental interval along direction d
  // coordinate = number to periodize
  // dimension = index of dimension
  inline int Periodize(int coordinate, int dimension);

  // periodize index within fundamental interval along direction d
  // coordinate = number to periodize
  // dimension = index of dimension
  // shift = translation of coordinate necessary to end up in unit cell
  inline int Periodize(int coordinate, int dimension, int &shift);  


  // find spins within a plaquette of the lattice
  // origin = number of spin where plaquette originates
  // vec = coordinates in lattice vectors of direction
  // length = number of spins in plaquette
  int *DeterminePlaquetteSpins(int origin, double *vec, int & length);
  
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

// periodize index within fundamental interval along direction d
// coordinate = number to periodize
// dimension = index of dimension
// shift = translation of coordinate necessary to end up in unit cell, in units of lattice vectors
int LatticeConnections::Periodize(int coordinate, int dimension, int &shift)
{
  int result;
  shift=0;
  //std::cout << "Raw value: "<<coordinate;
  if (coordinate<0)
    {
      shift = (coordinate/PeriodicRep[dimension] + 1)*PeriodicRep[dimension];
      coordinate += shift;
    }
  shift += (result=(coordinate%PeriodicRep[dimension])) - coordinate;
  //std::cout << ", shift="<<shift<<", result="<<result<<std::endl;
  return result;
}


#endif
