////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                           class author: Yang-Le Wu                         //
//                                                                            //
//      class of MPS matrix for quasihole in the Read-Rezayi k=3 state        //
//                                                                            //
//                        last modification : 07/12/2013                      //
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


#ifndef FQHEMPSREADREZAYI3QUASIHOLEMATRIX_H
#define FQHEMPSREADREZAYI3QUASIHOLEMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3Matrix.h"
#include "Vector/LongRationalVector.h"


class LongRationalMatrix;
class BosonOnDiskShort;

class FQHEMPSReadRezayi3QuasiholeMatrix : public FQHEMPSReadRezayi3Matrix
{

protected:

public:

    // constructor 
    //
    // pLevel = |P| level truncation
    // nbrBMatrices = number of B matrices to compute (max occupation per orbital)
    // cftDirectory = path to the directory where all the pure CFT matrices are stored
    // useRational = use arbitrary precision numbers for all the CFT calculations
    // trimChargeIndices = trim the charge indices
    // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
    // kappa = cylinder aspect ratio
    // architecture = architecture to use for precalculation
    FQHEMPSReadRezayi3QuasiholeMatrix(int pLevel, int nbrBMatrices, char* cftDirectory = 0,
            bool useRational = true, bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
            AbstractArchitecture* architecture = 0);

    // destructor
    //
    ~FQHEMPSReadRezayi3QuasiholeMatrix();

    // create the B matrices for the sigma field
    //
    // cftDirectory = an optional path to the directory where all the CFT matrices are stored
    // architecture = architecture to use for precalculation
    void CreateBMatrices(char* cftDirectory, AbstractArchitecture* architecture);

    // dummy
    virtual bool SaveMatrices(char* fileName);
};


#endif

