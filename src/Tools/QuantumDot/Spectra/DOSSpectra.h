#ifndef DOSSPECTRA_H
#define DOSSPECTRA_H

#include "config.h"

#include "Tools/QuantumDot/Spectra/Spectra.h"

class DOSSpectra : public Spectra
{
 public:

  // constructor from a set of energy files. Each peak is assimilated to a Lorentzian function.
  //
  // FileNumber: number of files, Files: name of files
  // StateNumber: integer array containing number of states in each file
  // Gamma: FWHM
  // Emin, Emax, dE: two energy bounds and the step in energy
  DOSSpectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE);
};

#endif
