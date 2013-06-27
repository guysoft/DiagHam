////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of MPS matrix for quasihole in the clustered (k=2,r) states      //
//                                                                            //
//                        last modification : 25/06/2013                      //
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


#ifndef FQHEMPSCLUSTERED2RQUASIHOLEMATRIX_H
#define FQHEMPSCLUSTERED2RQUASIHOLEMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RQuasiholeSectorMatrix.h"
#include "Vector/LongRationalVector.h"


class LongRationalMatrix;
class BosonOnDiskShort;

class FQHEMPSClustered2RQuasiholeMatrix : public FQHEMPSClustered2RQuasiholeSectorMatrix
{

protected:

public:

    // default constructor
    //
    FQHEMPSClustered2RQuasiholeMatrix();

    // constructor
    //
    // rindex = r index (i.e. clustered (k=2,r) states)
    // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)
    // pLevel = |P| level truncation
    // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
    // useRational = use arbitrary precision numbers for all the CFT calculations
    // trimChargeIndices = trim the charge indices
    // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
    // kappa = cylinder aspect ratio
    // architecture = architecture to use for precalculation
    FQHEMPSClustered2RQuasiholeMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool useRational = true,
            bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0,
            AbstractArchitecture* architecture = 0);

    // constructor
    //
    // rindex = r index (i.e. clustered (k=2,r) states)
    // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)
    // pLevel = |P| level truncation
    // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
    // useRational = use arbitrary precision numbers for all the CFT calculations
    // trimChargeIndices = trim the charge indices
    // cftDirectory = path to the directory where all the pure CFT matrices are stored
    // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
    // kappa = cylinder aspect ratio
    // architecture = architecture to use for precalculation
    FQHEMPSClustered2RQuasiholeMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool useRational = true,
            bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0,
            AbstractArchitecture* architecture = 0);

    // constructor from a file describing the state
    //
    // pLevel = |P| level truncation
    // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
    // fileName = name of the file that contains the state description
    // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
    // kappa = cylinder aspect ratio
    // architecture = architecture to use for precalculation
    FQHEMPSClustered2RQuasiholeMatrix(int pLevel, int nbrBMatrices, char* fileName, bool cylinderFlag = false, double kappa = 1.0,
            AbstractArchitecture* architecture = 0);

    // constructor from stored B matrices
    //
    // rindex = r index (i.e. clustered (k=2,r) states)
    // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)
    // pLevel = |P| level truncation
    // fileName = name of the file that contains the B matrices
    // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
    // kappa = cylinder aspect ratio
    FQHEMPSClustered2RQuasiholeMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag = false, double kappa = 1.0);

    // destructor
    //
    ~FQHEMPSClustered2RQuasiholeMatrix();

    // create the B matrices for the sigma field
    //
    // cftDirectory = an optional path to the directory where all the CFT matrices are stored
    // architecture = architecture to use for precalculation
    void CreateBMatrices(char* cftDirectory, AbstractArchitecture* architecture);

    // dummy
    virtual bool SaveMatrices(char* fileName);
};


#endif
