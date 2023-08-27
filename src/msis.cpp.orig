// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
//
// initial version - A. Ridley - April 30, 2023

#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------------
// Call msis model if we enable this in the cmake file
// -----------------------------------------------------------------------------

extern "C" void init_msis(void);
extern "C" void call_msis_f(int *iYear,
                            int *iDay,
                            float *second,
                            float *gLonDeg,
                            float *gLatDeg,
                            float *altKm,
                            float *f107,
                            float *f107a,
                            float *ap,
                            float[], float[]);

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

Msis::Msis() {

#ifdef FORTRAN
  // Initialize msis (reading in the data file):
  init_msis();
  isCompiled = true;
#else
  isCompiled = false;
#endif

  // This is the order in which we will store values:
  if (isCompiled) {
    value_lookup["He"] = 0;
    value_lookup["O"] = 1;
    value_lookup["N2"] = 2;
    value_lookup["O2"] = 3;
    value_lookup["Ar"] = 4;
    value_lookup["H"] = 5;
    value_lookup["N"] = 6;
    value_lookup["N_4S"] = 6;
    value_lookup["NO"] = 7;
    value_lookup["Tn"] = 8;
    nVars = 9;
  } else {
    value_lookup["unknown"] = 0;
    nVars = 0;
  }

  didChange = true;
  nX = -1;
  nY = -1;
  nZ = -1;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

bool Msis::is_ok() {
  return isCompiled;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

bool Msis::is_valid_species(std::string cValue) {
  if (value_lookup.contains(cValue))
    return true;
  else
    return false;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

bool Msis::set_time(Times time) {
  bool didWork = true;
  std::vector<int> iCurrent = time.get_iCurrent();
  iYear = iCurrent[0];
  iDay = time.get_julian_day();
  second =
    float(iCurrent[3]) * 3600.0 +
    float(iCurrent[4]) * 60.0 +
    float(iCurrent[5]) +
    float(iCurrent[6]) / 1000.0;
  didChange = true;
  return didWork;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

bool Msis::set_f107(precision_t f107in, precision_t f107ain) {
  bool didWork = true;
  f107 = f107in;
  f107a = f107ain;
  didChange = true;
  return didWork;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

bool Msis::set_ap(precision_t apin) {
  bool didWork = true;
  ap = apin;
  didChange = true;
  return didWork;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

bool Msis::reset_interface_variable_sizes() {
  bool didWork = true;
  lonDeg.set_size(nX, nY, nZ);
  latDeg.set_size(nX, nY, nZ);
  altKm.set_size(nX, nY, nZ);
  msis_results = make_cube_vector(nX, nY, nZ, nVars);
  return didWork;
}

// -----------------------------------------------------------------------------
// longitude is in radians
// latitude is in radians
// altitude is in meters
// -----------------------------------------------------------------------------

bool Msis::set_locations(arma_vec longitude,
                         arma_vec latitude,
                         arma_vec altitude) {
  bool didWork = true;

  // Set the size of the results if we need to:
  if (longitude.n_rows != nX) {
    nX = longitude.n_rows;
    nY = 1;
    nZ = 1;
    didWork = reset_interface_variable_sizes();
  }

  lonDeg.tube(0, 0) = longitude * cRtoD;
  latDeg.tube(0, 0) = latitude * cRtoD;
  altKm.tube(0, 0) = altitude / 1000.0;
  didChange = true;

  return didWork;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

bool Msis::set_locations(arma_mat longitude,
                         arma_mat latitude,
                         arma_mat altitude) {
  bool didWork = true;

  // Set the size of the results if we need to:
  if (longitude.n_rows != nX ||
      longitude.n_cols != nY) {
    nX = longitude.n_rows;
    nY = longitude.n_cols;
    nZ = 1;
    didWork = reset_interface_variable_sizes();
  }

  lonDeg.slice(0) = longitude * cRtoD;
  latDeg.slice(0) = latitude * cRtoD;
  altKm.slice(0) = altitude / 1000.0;
  didChange = true;

  return didWork;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

bool Msis::set_locations(arma_cube longitude,
                         arma_cube latitude,
                         arma_cube altitude) {
  bool didWork = true;

  // Set the size of the results if we need to:
  if (longitude.n_rows != nX ||
      longitude.n_cols != nY ||
      longitude.n_slices != nZ) {
    nX = longitude.n_rows;
    nY = longitude.n_cols;
    nZ = longitude.n_slices;
    didWork = reset_interface_variable_sizes();
  }

  lonDeg = longitude * cRtoD;
  latDeg = latitude * cRtoD;
  altKm = altitude / 1000.0;
  didChange = true;

  return didWork;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

bool Msis::reset_results() {

  bool didWork = true;
  float density_back[10];
  float temperature_back[2];
  float gAltKm;
  float gLatDeg;
  float gLonDeg;

  for (int i = 0; i < 10; i++)
    density_back[i] = -1.0;

  temperature_back[0] = -1.0;
  temperature_back[1] = -1.0;

  int64_t iX, iY, iZ;

  for (iX = 0; iX < nX; iX++) {
    for (iY = 0; iY < nY; iY++) {
      for (iZ = 0; iZ < nZ; iZ++) {
        gLonDeg = lonDeg(iX, iY, iZ);
        gLatDeg = latDeg(iX, iY, iZ);
        gAltKm = altKm(iX, iY, iZ);

#ifdef FORTRAN
        // Call msis in fortran:
        call_msis_f(&iYear, &iDay, &second, &gLonDeg, &gLatDeg, &gAltKm,
                    &f107, &f107a, &ap, density_back, temperature_back);
#endif
        // convert from /cm3 to /m3
        msis_results[0](iX, iY, iZ) = density_back[0] * 1e6;
        msis_results[1](iX, iY, iZ) = density_back[1] * 1e6;
        msis_results[2](iX, iY, iZ) = density_back[2] * 1e6;
        msis_results[3](iX, iY, iZ) = density_back[3] * 1e6;
        msis_results[4](iX, iY, iZ) = density_back[4] * 1e6;
        msis_results[5](iX, iY, iZ) = density_back[6] * 1e6;
        msis_results[6](iX, iY, iZ) = density_back[7] * 1e6;
        msis_results[7](iX, iY, iZ) = density_back[9] * 1e6;
        msis_results[8](iX, iY, iZ) = temperature_back[1];
      } // iZ
    } // iY
  } // iX

  return didWork;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

arma_vec Msis::get_vec(std::string cVar) {
  bool didWork = true;
  int item_ = -1;

  if (didChange)
    didWork = reset_results();

  if (didWork)
    item_ = value_lookup[cVar];

  return msis_results[item_].tube(0, 0);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

arma_mat Msis::get_mat(std::string cVar) {
  bool didWork = true;
  int item_ = -1;

  if (didChange)
    didWork = reset_results();

  if (didWork)
    item_ = value_lookup[cVar];

  return msis_results[item_].slice(0);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

arma_cube Msis::get_cube(std::string cVar) {
  bool didWork = true;
  int item_ = -1;

  if (didChange)
    didWork = reset_results();

  if (didWork)
    item_ = value_lookup[cVar];

  return msis_results[item_];
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

