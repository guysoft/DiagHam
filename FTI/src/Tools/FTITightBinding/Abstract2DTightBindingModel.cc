////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 2D tight binding model                  //
//                                                                            //
//                        last modification : 13/10/2012                      //
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


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/Endian.h"

#include <fstream>
#include <iostream>
#include <sys/time.h>

using std::ofstream;
using std::endl;
using std::cout;
using std::max;
using std::min;

#ifdef __FFTW__
#include <fftw3.h>
#endif

#ifdef HAVE_LAPACK
// solve general linear system
extern "C" void FORTRAN_NAME(dgesv)(const int* n, const int* nrhs, const double* a, const int* lda, const int* ipiv, const double* b, const int* ldb, const int* info );
#endif


// default constructor
//

Abstract2DTightBindingModel::Abstract2DTightBindingModel()
{
  this->EmbeddingY = RealVector();
  this->TwistAngle = M_PI / 2;
  this->Curvature = NULL;
  this->Chern = NULL;
  this->LLLGammaX = NULL;
  this->LLLGammaY = NULL;
}

// destructor
//

Abstract2DTightBindingModel::~Abstract2DTightBindingModel()
{
  delete[] this->Curvature;
  delete[] this->Chern;
  delete[] this->LLLGammaX;
  delete[] this->LLLGammaY;
}

// write an header that describes the tight binding model
// 
// output = reference on the output stream
// return value  = reference on the output stream

ofstream& Abstract2DTightBindingModel::WriteHeader(ofstream& output)
{
  int Dimension = 2;
  int HeaderSize = (((this->NbrBands + 2) * Dimension + 1) * sizeof(double)) + ((Dimension + 1) * sizeof(int));
  WriteLittleEndian(output, HeaderSize);
  WriteLittleEndian(output, Dimension);
  WriteLittleEndian(output, this->NbrSiteX);
  WriteLittleEndian(output, this->KxFactor);
  WriteLittleEndian(output, this->GammaX);
  if (this->EmbeddingX.GetVectorDimension() != this->NbrBands)
  {
      double Tmp = 0.0;
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, Tmp);
  }
  else
  {
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, this->EmbeddingX[i]);
  }
  WriteLittleEndian(output, this->NbrSiteY);
  WriteLittleEndian(output, this->KyFactor);
  WriteLittleEndian(output, this->GammaY);
  if (this->EmbeddingY.GetVectorDimension() != this->NbrBands)
  {
      double Tmp = 0.0;
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, Tmp);
  }
  else
  {
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, this->EmbeddingY[i]);
  }
  WriteLittleEndian(output, this->TwistAngle);
  return output; 
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool Abstract2DTightBindingModel::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(kx, ky);
	  File << kx << " " << ky; 
	  for (int i = 0; i < this->NbrBands; ++i)
	    File << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex];
	  File << endl;
	}
    }
  File.close();
  return true;
}

// write the full band structure information in an ASCII file
//
// fileName = name of the output file 
// return value = true if no error occured  

bool Abstract2DTightBindingModel::WriteBandStructureASCII(char* fileName)
{
  ofstream File;
  File.open(fileName);
  File.precision(14);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  for (int i = 0; i < this->NbrBands; ++i)
    for (int j = 0; j < this->NbrBands; ++j)
      File <<  "    U_{" << i << ", " << j << "}";
  File << endl;
  Complex Tmp;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(kx, ky);
	  File << kx << " " << ky; 
	  for (int i = 0; i < this->NbrBands; ++i)
	    File << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex];
	  for (int i = 0; i < this->NbrBands; ++i)
	    for (int j = 0; j < this->NbrBands; ++j)
	      {
		this->GetOneBodyMatrix(LinearizedMomentumIndex).GetMatrixElement(i, j, Tmp);
		File <<  "    " << Tmp;
	      }
	  File << endl;
	}
    }

  File.close();
  return true;
}

// compute the exponentiated, unitary Abelian connection
//
// kx = momentum along x
// ky = momentum along y
// qx = momentum transfer along x
// qy = momentum transfer along y
// band = band index
// return value = < u(k) | u(k+q) >
Complex Abstract2DTightBindingModel::GetAbelianConnection(int kx, int ky, int qx, int qy, int band)
{
    // [band][orbital] = Conj(<orbital|band>)
    int k1 = this->GetLinearizedMomentumIndexSafe(kx, ky);
    int k2 = this->GetLinearizedMomentumIndexSafe(kx + qx, ky + qy);
    ComplexVector bra = this->GetOneBodyMatrix(k1)[band];
    ComplexVector ket = this->GetOneBodyMatrix(k2)[band];

    Complex inner = 0.0;
    if (this->EmbeddingX.GetVectorDimension() != this->NbrBands)
        inner = ket * bra;
    else
    {
        for (int i = 0; i < this->NbrBands; ++i)
            inner += Phase(- 2.0 * M_PI * (qx * this->EmbeddingX[i] / this->NbrSiteX + qy * this->EmbeddingY[i] / this->NbrSiteY)) * bra[i] * Conj(ket[i]);
    }
    double n = Norm(inner);
    if (n < 1e-13)
    {
        cout << "Cannot make connection unitary due to orthogonality!" << endl;
        exit(1);
    }
    return inner / n;
}

// compute the unitary Abelian Wilson loop
//
// ky = momentum along y
// band = band index
// return value = value of the Wilson loop

Complex Abstract2DTightBindingModel::GetAbelianWilsonLoopX(int ky, int band)
{
    Complex prod = 1.0;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        prod *= this->GetAbelianConnection(kx, ky, 1, 0, band);
    return prod;
}

// compute the unitary Abelian Wilson loop
//
// kx = momentum along x
// band = band index
// return value = value of the Wilson loop

Complex Abstract2DTightBindingModel::GetAbelianWilsonLoopY(int kx, int band)
{
    Complex prod = 1.0;
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
        prod *= this->GetAbelianConnection(kx, ky, 0, 1, band);
    return prod;
}

// compute the stream function for the part of Berry connections that accounts for curvature fluctuations
//
// band = band index
// phi = reference to the vector storing the stream function over linearized BZ
// return = 0 if succeed, otherwise fail

int Abstract2DTightBindingModel::ComputeStreamFunction(int band, RealVector& phi, double& vx, double& vy)
{
    if (this->Curvature == NULL)
        this->ComputeCurvature();


#ifdef __FFTW__
    fftw_complex* source = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->NbrStatePerBand);
    fftw_complex* fsource = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->NbrStatePerBand);

    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int k = kx * this->NbrSiteY + ky; // need to enforce row-major layout
            source[k][0] = this->Curvature[band][ky][kx] - ((double) this->Chern[band]) / this->NbrStatePerBand;
            source[k][1] = 0.0;
        }
    fftw_plan plan_source = fftw_plan_dft_2d(this->NbrSiteX, this->NbrSiteY, source, fsource, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_source);
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            fsource[kx * this->NbrSiteY + ky][0] /= this->NbrStatePerBand;
            fsource[kx * this->NbrSiteY + ky][1] /= this->NbrStatePerBand;
        }


    fftw_complex* solution = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->NbrStatePerBand);
    fftw_complex* fsolution = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->NbrStatePerBand);

    fsolution[0][0] = 0.0;
    fsolution[0][1] = 0.0;
    for (int nx = 0; nx < this->NbrSiteX; ++nx)
    {
        for (int ny = 0; ny < this->NbrSiteY; ++ny)
        {
            int n = nx * this->NbrSiteY + ny;
            if (n == 0)
                continue;

            double nsq = - 4.0;
            nsq += 2 * cos(2 * M_PI * ((double) nx) / this->NbrSiteX);
            nsq += 2 * cos(2 * M_PI * ((double) ny) / this->NbrSiteY);

            fsolution[n][0] = fsource[n][0] / nsq;
            fsolution[n][1] = fsource[n][1] / nsq;
        }
    }

    fftw_plan plan_solution = fftw_plan_dft_2d(this->NbrSiteX, this->NbrSiteY, fsolution, solution, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_solution);

    double maxerror = 0.0;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        int kxm = (kx - 1 + this->NbrSiteX) % this->NbrSiteX;
        int kxp = (kx + 1) % this->NbrSiteX;
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int kym = (ky - 1 + this->NbrSiteY) % this->NbrSiteY;
            int kyp = (ky + 1) % this->NbrSiteY;

            int k = kx * this->NbrSiteY + ky;
            int ku = kx * this->NbrSiteY + kyp;
            int kd = kx * this->NbrSiteY + kym;
            int kl = kxm * this->NbrSiteY + ky;
            int kr = kxp * this->NbrSiteY + ky;

            if (fabs(solution[k][1]) > maxerror)
                maxerror = fabs(solution[k][1]);
            double error = fabs(solution[kl][0] + solution[kr][0] + solution[ku][0] + solution[kd][0] - 4 * solution[k][0] - source[k][0]);
            if (error > maxerror)
                maxerror = error;
        }
    }
    cout << "FFTW-Poisson error: " << maxerror << endl;

    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
            phi[this->GetLinearizedMomentumIndexSafe(kx, ky)] = solution[kx * this->NbrSiteY + ky][0];

    fftw_destroy_plan(plan_source);
    fftw_destroy_plan(plan_solution);
    fftw_free(source);
    fftw_free(fsource);
    fftw_free(solution);
    fftw_free(fsolution);

#else
#ifdef HAVE_LAPACK
    // curvature fluctuations
    double* raw = new double[this->NbrStatePerBand];
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
            raw[this->GetLinearizedMomentumIndexSafe(kx, ky)] = this->Curvature[band][ky][kx] - ((double) this->Chern[band]) / this->NbrStatePerBand;

    // laplacian with PBC
    double* laplacian = new double[this->NbrStatePerBand * this->NbrStatePerBand];
    for (int k = 0; k < this->NbrStatePerBand; ++k)
    {
        for (int l = 0; l < this->NbrStatePerBand; ++l)
        {
            laplacian[k + l * this->NbrStatePerBand] = 0.0;
        }
        laplacian[k + k * this->NbrStatePerBand] = -4.0;
    }
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int k = this->GetLinearizedMomentumIndexSafe(kx, ky);
            laplacian[k + this->GetLinearizedMomentumIndexSafe(kx, ky + 1) * this->NbrStatePerBand] = 1.0;
            laplacian[k + this->GetLinearizedMomentumIndexSafe(kx, ky - 1) * this->NbrStatePerBand] = 1.0;
            laplacian[k + this->GetLinearizedMomentumIndexSafe(kx + 1, ky) * this->NbrStatePerBand] = 1.0;
            laplacian[k + this->GetLinearizedMomentumIndexSafe(kx - 1, ky) * this->NbrStatePerBand] = 1.0;
        }
    }

    // solve poisson: \nabla^2 phi = curvature - 1/Nf
    int n = this->NbrStatePerBand;
    int nrhs = 1;
    int* ipiv = new int[this->NbrStatePerBand];
    int info = 42;
    FORTRAN_NAME(dgesv)(&n, &nrhs, laplacian, &n, ipiv, raw, &n, &info);
    delete[] ipiv;
    if (info)
    {
        cout << "Lapack DGESV failed with exit code " << info << endl;
        return -1;
    }
    for (int k = 0; k < this->NbrStatePerBand; ++k)
        phi[k] = raw[k];
    delete[] laplacian;
    delete[] raw;
#else
    cout << "Need FFTW or LAPACK to solve the Poisson equation for stream function." << endl;
    return -1;
#endif
#endif

    // overall shift in Ax and Ay, to reset Wx(0) and Wy(0) to 1
    vx = 0.0;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        vx += phi[this->GetLinearizedMomentumIndexSafe(kx, 0)] - phi[this->GetLinearizedMomentumIndexSafe(kx, -1)];
    vx /= this->NbrSiteX;

    vy = 0.0;
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
        vy -= phi[this->GetLinearizedMomentumIndexSafe(0, ky)] - phi[this->GetLinearizedMomentumIndexSafe(-1, ky)];
    vy /= this->NbrSiteY;

    return 0;
}

// build the gauge transform such that gauge(k) * |k⟩_lat is in the "Γ"-shaped parallel-transport gauge
// 
// band = band index
// gauge = reference to the gauge transform

void Abstract2DTightBindingModel::BuildParallelTransportGauge(int band, ComplexMatrix& gauge)
{
    if (this->Curvature == NULL)
        this->ComputeCurvature();

    gauge[0][0] = 1.0;

    // p-t along ky
    // we don't evenly spread Wy(0) over the bonds, in order to match with the quasi-periodic LLL gauge
    // Wy(0) is thus hidden in the bond crossing periodic boundary
    for (int ky = 1; ky < this->NbrSiteY; ++ky)
        gauge[ky][0] = gauge[ky - 1][0] / this->GetAbelianConnection(0, ky - 1, 0, 1, band);

    // p-t with rotation along kx
    // need to fit the phase of the LLL connection Ay = e^{- i 2 Pi (ky + gammay) / Nf}
    double* theta = new double[this->NbrSiteY];
    theta[0] = Arg(this->GetAbelianWilsonLoopX(0, band)); // need to be consistent with GetLLLGammaY
    if (this->Chern[band] > 0)
    {
        for (int ky = 1; ky < this->NbrSiteY; ++ky)
        {
            theta[ky] = Arg(this->GetAbelianWilsonLoopX(ky, band));
            while (theta[ky] <= theta[ky - 1])
                theta[ky] += 2 * M_PI;
            while (theta[ky] > theta[ky - 1])
                theta[ky] -= 2 * M_PI;
        }
    }
    else
    {
        for (int ky = 1; ky < this->NbrSiteY; ++ky)
        {
            theta[ky] = Arg(this->GetAbelianWilsonLoopX(ky, band));
            while (theta[ky] >= theta[ky - 1])
                theta[ky] -= 2 * M_PI;
            while (theta[ky] < theta[ky - 1])
                theta[ky] += 2 * M_PI;
        }
    }

    cout << "Theta(ky) / (2 Pi)" << endl;
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
        cout << "  " << theta[ky] / (2 * M_PI)<< endl;
    double total = fabs(theta[this->NbrSiteY - 1] - theta[0]) / (2 * M_PI);
    if ((total > abs(this->Chern[band])) || (total < (abs(this->Chern[band]) - 1)))
        cout << "Bad branch choice!" << endl;
    cout << endl;

    for (int ky = 0; ky < this->NbrSiteY; ++ky)
    {
        Complex lambda = Phase(theta[ky] / this->NbrSiteX);
        for (int kx = 1; kx < this->NbrSiteX; ++kx)
            gauge[ky][kx] = gauge[ky][kx - 1] * lambda / this->GetAbelianConnection(kx - 1, ky, 1, 0, band);
    }
    delete[] theta;
}

// build the gauge transform such that gauge(k) * |k⟩_lat is in the generalized π/2-rotated Landau gauge
//
// band = band index
// gauge = reference to the gauge transform
// return = 0 if succeed, otherwise fail

int Abstract2DTightBindingModel::BuildGeneralizedLandauGauge(int band, ComplexMatrix& gauge)
{
    double vx, vy;
    // matrices are initialized by (nrow, ncol) and accessed by [col][row]
    RealMatrix ax(this->NbrSiteX, this->NbrSiteY, true);
    RealMatrix ay(this->NbrSiteX, this->NbrSiteY, true);
    RealVector phi(this->NbrSiteX * this->NbrSiteY, true);

    // handle curvature fluctuations
    if (this->Curvature == NULL)
        this->ComputeCurvature();
    if (this->ComputeStreamFunction(band, phi, vx, vy) != 0)
        return -1;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            ax[ky][kx] += vx + phi[this->GetLinearizedMomentumIndexSafe(kx, ky - 1)] - phi[this->GetLinearizedMomentumIndexSafe(kx, ky)];
            ay[ky][kx] += vy + phi[this->GetLinearizedMomentumIndexSafe(kx, ky)] - phi[this->GetLinearizedMomentumIndexSafe(kx - 1, ky)];
        }
    }

    // twisting angles needed to handle Wx(0) and Wy(0)
    double gammaX = this->GetLLLGammaX(band);
    double gammaY = this->GetLLLGammaY(band);

    // handle curvature average
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            ax[ky][kx] += - this->Chern[band] * (ky + gammaY) / this->NbrStatePerBand;
            if (ky == this->NbrSiteY - 1)
                ay[ky][kx] += this->Chern[band] * (kx + gammaX) / this->NbrSiteX;
        }
    }

    // construct gauge transform to the target gauge given by A?Phase, gauge(k) * |k⟩_lat
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
            gauge[ky][kx] = 1.0;
    for (int kx = 1; kx < this->NbrSiteX; ++kx)
        gauge[0][kx] = gauge[0][kx - 1] * Phase(2 * M_PI * ax[0][kx - 1]) / this->GetAbelianConnection(kx - 1, 0, 1, 0, band);
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 1; ky < this->NbrSiteY; ++ky)
            gauge[ky][kx] = gauge[ky - 1][kx] * Phase(2 * M_PI * ay[ky - 1][kx]) / this->GetAbelianConnection(kx, ky - 1, 0, 1, band);
    }

    // check gauge construction
    bool Fail = false;

    // check connections against plaquette Wilson loops
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            double curvature = ax[ky][kx] + ay[ky][(kx + 1) % this->NbrSiteX] - ax[(ky + 1) % this->NbrSiteY][kx] - ay[ky][kx];
            if ((kx == this->NbrSiteX - 1) && (ky == this->NbrSiteY - 1)) // trivial fix. does not affect phase
                curvature += Chern[band];

            double diff = fabs(curvature - this->Curvature[band][ky][kx]);
            if (diff > 1e-13)
            {
                Fail = true;
                cout << "Curvature Mismatch @ k = (" << kx << "," << ky << "), Error = " << diff << endl;
            }
        }
    }

    // check connections against large Wilson loops
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
    {
        double sum = 0;
        for (int kx = 0; kx < this->NbrSiteX; ++kx)
            sum += ax[ky][kx];
        double diff = Arg(Phase(2 * M_PI * sum) / this->GetAbelianWilsonLoopX(ky, band));
        if (diff > 1e-13)
        {
            Fail = true;
            cout << "Wx Mismatch @ ky = " << ky << ", Error = " << diff << endl;
        }
    }
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        double sum = 0;
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
            sum += ay[ky][kx];
        double diff = Arg(Phase(2 * M_PI * sum) / this->GetAbelianWilsonLoopY(kx, band));
        if (diff > 1e-13)
        {
            Fail = true;
            cout << "Wy Mismatch @ kx = " << kx << ", Error = " << diff << endl;
        }
    }

    // check gauge transform
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        int kxm = (kx - 1 + this->NbrSiteX) % this->NbrSiteX;
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int kym = (ky - 1 + this->NbrSiteY) % this->NbrSiteY;
            double diff = Arg((gauge[ky][kx] / gauge[kym][kx]) / (Phase(2 * M_PI * ay[kym][kx]) / this->GetAbelianConnection(kx, kym, 0, 1, band)));
            if (diff > 1e-13)
            {
                Fail = true;
                cout << "Gauge Mismatch @ k = (" << kx << "," << ky << "), along y, Error = " << diff << endl;
            }
            diff = Arg((gauge[ky][kx] / gauge[ky][kxm]) / (Phase(2 * M_PI * ax[ky][kxm]) / this->GetAbelianConnection(kxm, ky, 1, 0, band)));
            if (diff > 1e-13)
            {
                Fail = true;
                cout << "Gauge Mismatch @ k = (" << kx << "," << ky << "), along x, Error = " << diff << endl;
            }
        }
    }

    return int(Fail);
}

// compute the curvature over each plaquette in the BZ, and also Chern number
//
// band = band index
// return = 0 if succeed, otherwise fail

void Abstract2DTightBindingModel::ComputeCurvature()
{
    if (this->Curvature == NULL)
        this->Curvature = new RealMatrix[this->NbrBands];
    if (this->Chern == NULL)
        this->Chern = new int[this->NbrBands];

    for (int b = 0; b < this->NbrBands; ++b)
        this->Curvature[b].ResizeAndClean(this->NbrSiteX, this->NbrSiteY);

    for (int b = 0; b < this->NbrBands; ++b)
    {
        for (int kx = 0; kx < this->NbrSiteX; ++kx)
            for (int ky = 0; ky < this->NbrSiteY; ++ky)
                this->Curvature[b][ky][kx] = this->ComputeCurvatureSinglePlaquette(kx, ky, b);

        double sum = 0.0;
        for (int kx = 0; kx < this->NbrSiteX; ++kx)
            for (int ky = 0; ky < this->NbrSiteY; ++ky)
                sum += this->Curvature[b][ky][kx];

        if (sum == 0.0)
            this->Chern[b] = 0;
        else
        {
            sum += (sum / fabs(sum)) * 0.5;
            this->Chern[b] = int(sum);
        }

        if (fabs(this->Chern[b] - (sum - (sum / fabs(sum)) * 0.5)) > 1e-13) // shouldn't happen unless connections are seriously screwed up
            cout << "non-integer Chern number for band " << b << "?! (curvature sum = " << sum << ")" << endl;
        if (this->Chern[b] == 0)
            cout << "Zero Chern number for band " << b << "?!" << endl;
    }

    if (this->LLLGammaX == NULL)
        this->LLLGammaX = new double[this->NbrBands];
    if (this->LLLGammaY == NULL)
        this->LLLGammaY = new double[this->NbrBands];
    for (int b = 0; b < this->NbrBands; ++b)
    {
        double gammaX = Arg(this->GetAbelianWilsonLoopY(0, b)); // Arg takes value in (-π, π]
        gammaX *= this->NbrSiteX / (2 * M_PI * this->Chern[b]);
        this->LLLGammaX[b] = gammaX;

        double gammaY = Arg(this->GetAbelianWilsonLoopX(0, b)); // Arg takes value in (-π, π]
        gammaY *= - this->NbrSiteY / (2 * M_PI * this->Chern[b]);
        this->LLLGammaY[b] = gammaY;
    }
}


// write the eigenvalues of the D matrix in an ASCII file
//
// fileName = name of the ASCII file 
//nbrOccupiedBands = nbr of occupied bands
// return value = true if no error occured

bool Abstract2DTightBindingModel::WriteAsciiDMatrixEigenValues(char* fileName, int nbrOccupiedBands)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# ky" ;
  for (int i = 0; i < nbrOccupiedBands; ++i)
    File <<  "   DEigenValue_" << i << "    Theta_" << i;
  File << endl;
  
  double distancePlus;
  double distanceMoins;
  double distanceMod2PiPlus;
  double distanceMod2PiMoins;
 
  double theta1;
  double theta2;
  
  Complex** Lambda = this->ComputeDMatrixEigenvalues(nbrOccupiedBands, 0, this->NbrSiteY - 1, this->NbrSiteY); 
  double** Theta = new double*[this->NbrSiteY];
  for (int ky = 0; ky < this->NbrSiteY; ++ky)
  {
   Theta[ky] = new double[2];
   for (int i = 0; i < 2; ++i)
   {
     theta1 = atan2(Lambda[ky][nbrOccupiedBands - 2].Im,Lambda[ky][nbrOccupiedBands - 2].Re);
     theta2 = atan2(Lambda[ky][nbrOccupiedBands - 1].Im,Lambda[ky][nbrOccupiedBands - 1].Re);
     Theta[ky][0] = max(theta1, theta2);
     Theta[ky][1] = min(theta1, theta2); 
   }
  }
  
  for (int ky = 0; ky < this->NbrSiteY - 1; ++ ky)
  {
    distancePlus = abs(Theta[ky][0] - Theta[ky + 1][0]);
    distanceMod2PiPlus = abs(Theta[ky][0] - Theta[ky + 1][1] - 2*M_PI);
    distanceMoins = abs(Theta[ky][1] - Theta[ky + 1][1]);
    distanceMod2PiMoins = abs(Theta[ky][1] - Theta[ky + 1][0] + 2*M_PI);
    
    if (distanceMod2PiPlus < distancePlus)
    {
     double Tmp = Theta[ky + 1][0];
     Theta[ky + 1][0] = Theta[ky + 1][1] + 2*M_PI;
     Theta[ky + 1][1] = Tmp;
    }
    
    if (distanceMod2PiMoins < distanceMoins)
    {
     double Tmp = Theta[ky + 1][1];
     Theta[ky + 1][1] = Theta[ky + 1][0] - 2*M_PI;
     Theta[ky + 1][0] = Tmp;
    }
  }
  for (int ky = 0; ky < this->NbrSiteY; ++ky)
    {
      File << ky; 
      File << " " << Lambda[ky][0] << " " << atan2(Lambda[ky][0].Im,Lambda[ky][0].Re) << " " << Lambda[ky][1] << " "<< atan2(Lambda[ky][1].Im,Lambda[ky][1].Re) << " " << Theta[ky][0] << " " << Theta[ky][1] ;
      File << endl;
    }
   
  File.close();
  return true;
}

// compute the Chern number of a given band
//
// band = band index
// return value = Chern number

double Abstract2DTightBindingModel::ComputeChernNumber(int band)
{
  if (this->HaveOneBodyBasis() == false)
    {
      cout << "error, the tight binding model does not provide the one body basis" << endl;
      return 0.0;
    }
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  Complex TmpChernNumber = 0.0;
  Complex Tmp1[4];
  Complex Tmp2[8];
  for (long LinearizedMomentumIndex = 0l; LinearizedMomentumIndex < this->NbrStatePerBand; ++LinearizedMomentumIndex)
    {
      int Kx;
      int Ky;
      this->GetLinearizedMomentumIndex(LinearizedMomentumIndex, Kx, Ky);
      int LinearizedMomentumIndexIncX = this->GetLinearizedMomentumIndex((Kx + 1) % this->NbrSiteX, Ky);
      int LinearizedMomentumIndexDecX;
      if (Kx > 0)
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex((Kx - 1) % this->NbrSiteX, Ky);
      else
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex(this->NbrSiteX - 1, Ky);
      int LinearizedMomentumIndexIncY = this->GetLinearizedMomentumIndex(Kx, (Ky + 1) % this->NbrSiteY);
      int LinearizedMomentumIndexDecY;
      if (Ky > 0)
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, (Ky - 1) % this->NbrSiteY);
      else
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, this->NbrSiteY - 1);

      ComplexMatrix& LocalBasis = this->OneBodyBasis[LinearizedMomentumIndex];
      ComplexMatrix& LocalBasisIncX = this->OneBodyBasis[LinearizedMomentumIndexIncX];
      ComplexMatrix& LocalBasisDecX = this->OneBodyBasis[LinearizedMomentumIndexDecX];
      ComplexMatrix& LocalBasisIncY = this->OneBodyBasis[LinearizedMomentumIndexIncY];
      ComplexMatrix& LocalBasisDecY = this->OneBodyBasis[LinearizedMomentumIndexDecY];  
      Tmp1[0] = 0.0;
      Tmp1[1] = 0.0;
      Tmp1[2] = 0.0;
      Tmp1[3] = 0.0;

      Tmp2[0] = 0.0;
      Tmp2[1] = 0.0;
      Tmp2[2] = 0.0;
      Tmp2[3] = 0.0;
      Tmp2[4] = 0.0;
      Tmp2[5] = 0.0;
      Tmp2[6] = 0.0;
      Tmp2[7] = 0.0;

      for (int i = 0; i < this->NbrBands; ++i)
	{
	  Tmp1[0] += LocalBasis[band][i] * Conj(LocalBasisIncX[band][i]);
	  Tmp1[1] += LocalBasis[band][i] * Conj(LocalBasisDecX[band][i]);
	  Tmp1[2] += LocalBasis[band][i] * Conj(LocalBasisIncY[band][i]);
	  Tmp1[3] += LocalBasis[band][i] * Conj(LocalBasisDecY[band][i]);

	  Tmp2[0] += Conj(LocalBasisIncX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[1] += Conj(LocalBasisDecX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[2] += Conj(LocalBasisIncX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[3] += Conj(LocalBasisDecX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[4] += Conj(LocalBasisIncY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[5] += Conj(LocalBasisDecY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[6] += Conj(LocalBasisIncY[band][i]) * LocalBasisDecX[band][i];
	  Tmp2[7] += Conj(LocalBasisDecY[band][i]) * LocalBasisDecX[band][i];
	}

      TmpChernNumber += (Tmp1[2] * Conj(Tmp1[0]) * Tmp2[0]);
      TmpChernNumber -= (Tmp1[2] * Conj(Tmp1[1]) * Tmp2[1]);
      TmpChernNumber -= (Tmp1[3] * Conj(Tmp1[0]) * Tmp2[2]);
      TmpChernNumber += (Tmp1[3] * Conj(Tmp1[1]) * Tmp2[3]);
	  
      TmpChernNumber -= (Tmp1[0] * Conj(Tmp1[2]) * Tmp2[4]);
      TmpChernNumber += (Tmp1[0] * Conj(Tmp1[3]) * Tmp2[5]);
      TmpChernNumber += (Tmp1[1] * Conj(Tmp1[2]) * Tmp2[6]);
      TmpChernNumber -= (Tmp1[1] * Conj(Tmp1[3]) * Tmp2[7]);

    }
  TmpChernNumber /= 8.0 * M_PI;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
  cout << "Chern number computed in  " << Dt << "s" << endl;
  return TmpChernNumber.Im;
}

// compute the Berry curvature  of a given band
//
// band = band index
// fileName = name of the output file 
// return value = Chern number

double Abstract2DTightBindingModel::ComputeBerryCurvature(int band, char* fileName)
{
  if (this->HaveOneBodyBasis() == false)
    {
      cout << "error, the tight binding model does not provide the one body basis" << endl;
      return 0.0;
    }
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  Complex TmpChernNumber = 0.0;
  Complex Tmp1[4];
  Complex Tmp2[8];
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky    Berry_curvature";
  double Fluctations = 0.0;
  for (long LinearizedMomentumIndex = 0l; LinearizedMomentumIndex < this->NbrStatePerBand; ++LinearizedMomentumIndex)
    {
      int Kx;
      int Ky;
      this->GetLinearizedMomentumIndex(LinearizedMomentumIndex, Kx, Ky);
      int LinearizedMomentumIndexIncX = this->GetLinearizedMomentumIndex((Kx + 1) % this->NbrSiteX, Ky);
      int LinearizedMomentumIndexDecX;
      if (Kx > 0)
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex((Kx - 1) % this->NbrSiteX, Ky);
      else
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex(this->NbrSiteX - 1, Ky);
      int LinearizedMomentumIndexIncY = this->GetLinearizedMomentumIndex(Kx, (Ky + 1) % this->NbrSiteY);
      int LinearizedMomentumIndexDecY;
      if (Ky > 0)
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, (Ky - 1) % this->NbrSiteY);
      else
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, this->NbrSiteY - 1);

      ComplexMatrix& LocalBasis = this->OneBodyBasis[LinearizedMomentumIndex];
      ComplexMatrix& LocalBasisIncX = this->OneBodyBasis[LinearizedMomentumIndexIncX];
      ComplexMatrix& LocalBasisDecX = this->OneBodyBasis[LinearizedMomentumIndexDecX];
      ComplexMatrix& LocalBasisIncY = this->OneBodyBasis[LinearizedMomentumIndexIncY];
      ComplexMatrix& LocalBasisDecY = this->OneBodyBasis[LinearizedMomentumIndexDecY];  
      Tmp1[0] = 0.0;
      Tmp1[1] = 0.0;
      Tmp1[2] = 0.0;
      Tmp1[3] = 0.0;

      Tmp2[0] = 0.0;
      Tmp2[1] = 0.0;
      Tmp2[2] = 0.0;
      Tmp2[3] = 0.0;
      Tmp2[4] = 0.0;
      Tmp2[5] = 0.0;
      Tmp2[6] = 0.0;
      Tmp2[7] = 0.0;

      for (int i = 0; i < this->NbrBands; ++i)
	{
	  Tmp1[0] += LocalBasis[band][i] * Conj(LocalBasisIncX[band][i]);
	  Tmp1[1] += LocalBasis[band][i] * Conj(LocalBasisDecX[band][i]);
	  Tmp1[2] += LocalBasis[band][i] * Conj(LocalBasisIncY[band][i]);
	  Tmp1[3] += LocalBasis[band][i] * Conj(LocalBasisDecY[band][i]);

	  Tmp2[0] += Conj(LocalBasisIncX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[1] += Conj(LocalBasisDecX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[2] += Conj(LocalBasisIncX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[3] += Conj(LocalBasisDecX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[4] += Conj(LocalBasisIncY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[5] += Conj(LocalBasisDecY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[6] += Conj(LocalBasisIncY[band][i]) * LocalBasisDecX[band][i];
	  Tmp2[7] += Conj(LocalBasisDecY[band][i]) * LocalBasisDecX[band][i];
	}

      Complex TmpCurvature = 0.0;
      TmpCurvature += (Tmp1[2] * Conj(Tmp1[0]) * Tmp2[0]);
      TmpCurvature -= (Tmp1[2] * Conj(Tmp1[1]) * Tmp2[1]);
      TmpCurvature -= (Tmp1[3] * Conj(Tmp1[0]) * Tmp2[2]);
      TmpCurvature += (Tmp1[3] * Conj(Tmp1[1]) * Tmp2[3]);
	  
      TmpCurvature -= (Tmp1[0] * Conj(Tmp1[2]) * Tmp2[4]);
      TmpCurvature += (Tmp1[0] * Conj(Tmp1[3]) * Tmp2[5]);
      TmpCurvature += (Tmp1[1] * Conj(Tmp1[2]) * Tmp2[6]);
      TmpCurvature -= (Tmp1[1] * Conj(Tmp1[3]) * Tmp2[7]);

      TmpCurvature *= 0.25;

      Fluctations += (TmpCurvature.Im - (2.0 * M_PI / ((double) this->NbrStatePerBand))) * (TmpCurvature.Im - (2.0 * M_PI / ((double) this->NbrStatePerBand)));

      File << Kx << " " << Ky << " " << TmpCurvature.Im << endl;

      TmpChernNumber += TmpCurvature;
    }
//  Fluctations *= ((double) this->NbrStatePerBand) ;
  cout << "Berry curvature fluctuations = " << Fluctations << " " <<  sqrt(Fluctations) << " " << (TmpChernNumber.Im * TmpChernNumber.Im) << " " << sqrt(Fluctations - (TmpChernNumber.Im * TmpChernNumber.Im)) 
       << "( " << (sqrt(Fluctations - (TmpChernNumber.Im * TmpChernNumber.Im)) / (2.0 * M_PI) )<< ") in 2 pi units" << endl;
  TmpChernNumber /= 2.0 * M_PI;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
  cout << "Chern number computed in  " << Dt << "s" << endl;
//  cout << "Berry curvature fluctuations = " << sqrt ()<< endl;

  File.close();

  return TmpChernNumber.Im;
}


//compute the complex eigenvalues of the D(ky) matrix (in order to compute the Z2 invariant)
//
//bandIndex = band index (corresponds to two bands that are related by time reversal symmetry)
//nbrOccupiedBands = dimension of the D matrix
//DMatrixEigenvalues = array of eigenvalues of the D Matrix, for all values of ky
//kyMin = minimal value of ky for which the D matrix has to be diagonalized
//kyMax = maximal value of ky for which the D matrix has to be diagonalized
//nbrKy = number of ky values for which the D matrix has to be diagonalized
//return value = array of eigenvalues of the D matrix
Complex** Abstract2DTightBindingModel::ComputeDMatrixEigenvalues(int nbrOccupiedBands, int kyMin, int kyMax, int nbrKy)
{
  Complex** DMatrixEigenvalues;
  DMatrixEigenvalues = new Complex*[this->NbrSiteY];
  for (int i = 0; i < this->NbrSiteY; ++i)
    DMatrixEigenvalues[i] = new Complex[nbrOccupiedBands];
  ComplexMatrix TmpDMatrix(nbrOccupiedBands, nbrOccupiedBands, true);
  ComplexMatrix FMatrix(nbrOccupiedBands, nbrOccupiedBands, true);
  
//   ComplexMatrix Rotation(this->NbrBands, this->NbrBands, true);
//   Rotation.SetMatrixElement(0, 0, M_SQRT1_2);
//   Rotation.SetMatrixElement(0, 1, Complex(0.0, M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 0, Complex(0.0, -1.0*M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 1, -1.0*M_SQRT1_2);
//   Rotation.SetMatrixElement(2, 2, 1.0);
//   Rotation.SetMatrixElement(3, 3, 1.0);
  
//   ComplexMatrix Rotation1(this->NbrBands, this->NbrBands, true);
//   Rotation.SetMatrixElement(0, 0, M_SQRT1_2);
//   Rotation.SetMatrixElement(0, 1, Complex(0.0, M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 0, Complex(0.0, -1.0*M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 1, -1.0*M_SQRT1_2);
//   Rotation1.SetMatrixElement(2, 2, 1.0);
//   Rotation1.SetMatrixElement(3, 3, 1.0);
  
  
  for (int ky = kyMin; ky <= kyMax; ++ky)
  {
    TmpDMatrix.SetToIdentity();
    for (int i = 0; i < this->NbrSiteX; ++i)
    {
      double KX = (double) i *2.0 * M_PI / ((double) this->NbrSiteX);
      double KX1 = (double) (i + 1) *2.0 * M_PI / ((double) this->NbrSiteX);
//       cout << cos(KX) << " " << sin(KX) << endl;
//       Rotation.SetMatrixElement(0, 0, Complex(M_SQRT1_2*cos(KX), M_SQRT1_2*sin(KX)));
//       Rotation.SetMatrixElement(0, 1, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 0, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 1, Complex(-1.0*M_SQRT1_2*cos(KX), M_SQRT1_2*sin(KX)));
// 
//       
//       Rotation1.SetMatrixElement(0, 0, Complex(M_SQRT1_2*cos(KX1), M_SQRT1_2*sin(KX1)));
//       Rotation1.SetMatrixElement(0, 1, M_SQRT1_2);
//       Rotation1.SetMatrixElement(1, 0, M_SQRT1_2);
//       Rotation1.SetMatrixElement(1, 1, Complex(-1.0*M_SQRT1_2*cos(KX1), M_SQRT1_2*sin(KX1)));
//       Rotation.SetMatrixElement(0, 0, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 1, -1.0*M_SQRT1_2);
//       Rotation.SetMatrixElement(0, 1, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 0, M_SQRT1_2);
 
	int LinearizedMomentumIndex1 = this->GetLinearizedMomentumIndex(i, ky);
	int LinearizedMomentumIndex2 = this->GetLinearizedMomentumIndex((i + 1) % this->NbrSiteX, ky);
// 	cout << i << " " << LinearizedMomentumIndex1 << " " << LinearizedMomentumIndex2 << " " << endl;
	ComplexMatrix& LocalBasis = this->OneBodyBasis[LinearizedMomentumIndex1];
	ComplexMatrix& LocalBasisIncX = this->OneBodyBasis[LinearizedMomentumIndex2];
// 	ComplexMatrix LocalBasis(this->NbrBands, this->NbrBands, true);
// 	ComplexMatrix LocalBasisIncX(this->NbrBands, this->NbrBands, true);
	
// 	LocalBasis = TmpLocalBasis*Rotation;
// 	LocalBasisIncX = TmpLocalBasisIncX*Rotation1;
	
// 	LocalBasis = TmpLocalBasis;
// 	LocalBasisIncX = TmpLocalBasisIncX;
	
	for (int n = 0; n < nbrOccupiedBands; ++n)
	  {
	    for (int m = 0; m < nbrOccupiedBands; ++m)
	      {
// 		Complex Tmp = 0.0;
// 		for (int alpha = 0; alpha < this->NbrBands; ++alpha)
// 		{
// 		  Tmp += Conj(LocalBasis[n][alpha]) * LocalBasisIncX[m][alpha];
// 		}
// 		FMatrix.SetMatrixElement(n, m, Tmp);
		FMatrix.SetMatrixElement(n, m, LocalBasis[n] * LocalBasisIncX[m]);
	      }
	  }
// 	  cout << i << endl;
// 	  cout << FMatrix << endl;
// 	  ComplexMatrix TmpMatrix = TmpDMatrix;
// 	  TmpDMatrix = TmpMatrix*FMatrix;
	  TmpDMatrix.Multiply(FMatrix);
      }
    
// 	  if (ky == 0)
// 	  {
// 	   cout << "ky = 0" << endl;
// 	   cout << TmpDMatrix << endl; 
// 	  }
    ComplexDiagonalMatrix TmpDiag(nbrOccupiedBands);
#ifdef __LAPACK__
    TmpDMatrix.LapackDiagonalize(TmpDiag);
#else
    TmpDMatrix.Diagonalize(TmpDiag);
#endif
//     cout << TmpDiag << endl;
    
    for (int j = 0; j < nbrOccupiedBands; ++j)
    {
      DMatrixEigenvalues[ky][j] = TmpDiag[j] ;
//       cout << DMatrixEigenvalues[ky][j] << endl;
    }
  }
  
  return DMatrixEigenvalues;
}

//compute the Z2 topological invariant for a system with time reversal symmetry
//
//nbrOccupiedBands = number of occupied bands
//return value = Z2 invariant
int Abstract2DTightBindingModel::ComputeZ2Invariant(int nbrOccupiedBands)
{
  int z2Invariant = 0;
  double referenceLine = 0.9267;
  
  double distancePlus;
  double distanceMoins;
  double distanceMod2PiPlus;
  double distanceMod2PiMoins;
  
  int ModPiPlus = 0;
  int ModPiMoins = 0;
  
  double theta1;
  double theta2;
  
  Complex** Lambda = this->ComputeDMatrixEigenvalues(nbrOccupiedBands, 0, this->NbrSiteY - 1, this->NbrSiteY); 
  double** Theta = new double*[this->NbrSiteY];
  for (int ky = 0; ky < this->NbrSiteY; ++ky)
  {
   Theta[ky] = new double[2];
   for (int i = 0; i < 2; ++i)
   {
     theta1 = atan2(Lambda[ky][nbrOccupiedBands - 2].Im,Lambda[ky][nbrOccupiedBands - 2].Re);
     theta2 = atan2(Lambda[ky][nbrOccupiedBands - 1].Im,Lambda[ky][nbrOccupiedBands - 1].Re);
     Theta[ky][0] = max(theta1, theta2);
     Theta[ky][1] = min(theta1, theta2); 
   }
  }
  
  for (int ky = 0; ky < this->NbrSiteY  - 1; ++ ky)
  {
    distancePlus = abs(Theta[ky][0] - Theta[ky + 1][0]);
    distanceMod2PiPlus = abs(Theta[ky][0] - Theta[ky + 1][1] - 2*M_PI);
    distanceMoins = abs(Theta[ky][1] - Theta[ky + 1][1]);
    distanceMod2PiMoins = abs(Theta[ky][1] - Theta[ky + 1][0] + 2*M_PI);
    
    if (distanceMod2PiPlus < distancePlus)
    {
     ModPiPlus += 1;
     double Tmp = Theta[ky + 1][0];
     Theta[ky + 1][0] = Theta[ky + 1][1] + 2*M_PI;
     Theta[ky + 1][1] = Tmp;
    }
    
    if (distanceMod2PiMoins < distanceMoins)
    {
     ModPiMoins += 1;
     double Tmp = Theta[ky + 1][1];
     Theta[ky + 1][1] = Theta[ky + 1][0] - 2*M_PI;
     Theta[ky + 1][0] = Tmp;
    }
  }
//   cout << ModPiPlus << " " << ModPiMoins << endl;
  for (int ky = 1 ; ky < this->NbrSiteY/2 - 1; ++ky)
  {
//     cout << ky << " " << Theta[ky][0] << " " << Theta[ky][1] << endl;
    for (int i = 0; i <= ModPiPlus; ++i)
    {
//       cout << Theta[ky + 1][0] - (referenceLine + i*2*M_PI) << " " << Theta[ky][0] - (referenceLine + i*2*M_PI) << " " << (Theta[ky + 1][1] - (referenceLine - i*2*M_PI)) << " " << (Theta[ky][1] - (referenceLine - i*2*M_PI)) << endl;
      if ((Theta[ky + 1][0] - (referenceLine + i*2*M_PI)) * (Theta[ky][0] - (referenceLine + i*2*M_PI)) < 0)
	z2Invariant += 1;
    }
    for (int i = 0; i <= ModPiMoins; ++i)
    {
      if ((Theta[ky + 1][1] - (referenceLine - i*2*M_PI)) * (Theta[ky][1] - (referenceLine - i*2*M_PI)) < 0)
	z2Invariant += 1;
    }
  }
  return (z2Invariant % 2); 
}