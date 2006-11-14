////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of SU(4) generalized Halperine wave function on sphere       //
//                                                                            //
//                        last modification : 12/11/2006                      //
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


#ifndef SU4HALPERINONSPHEREWAVEFUNCTION_H
#define SU4HALPERINONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"


class SU4HalperinOnSphereWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // total number of particles
  int TotalNbrParticles;

  // number of particles with spin up and isopsin plus
  int NbrSpinUpIsospinPlusParticles;
  // number of particles with spin up and isopsin minus
  int NbrSpinUpIsospinMinusParticles;
  // number of particles with spin down and isopsin plus
  int NbrSpinDownIsospinPlusParticles;
  // number of particles with spin down and isopsin minus
  int NbrSpinDownIsospinMinusParticles;

  // power of the laughlin-like part for spin up - isopsin plus
  int MupIndex;
  // power of the laughlin-like part for spin up - isopsin minus
  int MumIndex;
  // power of the laughlin-like part for spin down - isopsin plus
  int MdpIndex;
  // power of the laughlin-like part for spin down - isopsin minus
  int MdmIndex;
  // power of the intra-isospin part (i.e (z_u{pm} -z_d{pm}))
  int NIntraIsopinIndex;
  // power of the intra-spin part (i.e (z_{ud}p -z_{ud}m))
  int NIntraIsopinIndex;
  // power of the cross spin-isospin part (i.e (z_up -z_dm) and (z_um -z_dp))
  int NCrossSpinIsopinIndex;
  

 public:

  // constructor
  //
  // nbrSpinUpIsospinPlusParticles = number of particles with spin up and isopsin plus
  // nbrSpinUpIsospinMinusParticles = number of particles with spin up and isopsin minus
  // nbrSpinDownIsospinPlusParticles = number of particles with spin down and isopsin plus
  // nbrSpinDownIsospinMinusParticles = number of particles with spin down and isopsin minus
  // mupIndex = power of the laughlin-like part for spin up - isopsin plus
  // mumIndex = power of the laughlin-like part for spin up - isopsin minus
  // mdpIndex = power of the laughlin-like part for spin down - isopsin plus
  // mdmIndex = power of the laughlin-like part for spin down - isopsin minus
  // nIntraIsospinIndex = power of the intra-isospin part (i.e (z_u{pm} -z_d{pm}))
  // nIntraSpinIndex = power of the intra-spin part (i.e (z_{ud}p -z_{ud}m))
  // nCrossSpinIsopinIndex = power of the cross spin-isospin part (i.e (z_up -z_dm) and (z_um -z_dp))
  SU4HalperinOnSphereWaveFunction(int nbrSpinUpIsospinPlusParticles, int nbrSpinUpIsospinMinusParticles, 
				  int nbrSpinDownIsospinPlusParticles, int nbrSpinDownIsospinMinusParticles,
				  int mupIndex, int mumIndex, int mdpIndex, int mdmIndex
				  int nIntraIsopinIndex, int nIntraSpinIndex, int nCrossSpinIsopinIndex);

  // copy constructor
  //
  // function = reference on the wave function to copy
  SU4HalperinOnSphereWaveFunction(const SU4HalperinOnSphereWaveFunction& function);

  // destructor
  //
   ~HalperinOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point (the first 2*nbrSpinUpIsopinPlusParticles coordinates correspond to the position of the spin up - isposin plus particles, 
  //                                     the following 2*nbrSpinDownIsopinPlusParticles coordinates correspond to the position of spin down - isposin plus particles,
  //                                     the other coordinates obey to the same scheme for isopin minus)
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

};

#endif
