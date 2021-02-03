// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_EUV_H_
#define INCLUDE_EUV_H_

#include <vector>
#include <string>

#include "inputs.h"
#include "times.h"
#include "indices.h"
#include "planets.h"
#include "grid.h"
#include "ions.h"
#include "report.h"

class Euv {

public:

  int nWavelengths;
  int nLines;

  struct waveinfotype {
    std::string name;
    std::string to;
    std::string type;
    std::string units;
    std::string note;
    std::vector<float> values;
  };

  std::vector<waveinfotype> waveinfo;

  std::vector<float> wavelengths_short;
  std::vector<float> wavelengths_long;
  std::vector<float> wavelengths_energy;
  std::vector<float> wavelengths_intensity_1au;
  std::vector<float> wavelengths_intensity_top;

  std::vector<float> euvac_f74113;
  std::vector<float> euvac_afac;

  // --------------------------------------------------------------------------
  // Initialize EUV
  // --------------------------------------------------------------------------

  Euv(Inputs args, Report report);

  // -------------------------------------------------------------------------
  //
  // -------------------------------------------------------------------------

  int euvac(Times time, Indices indices, Report &report);

  // -------------------------------------------------------------------------
  //
  // -------------------------------------------------------------------------

  int scale_from_1au(Planets planet, Times time);

private:

  // --------------------------------------------------------------------------
  // Read in the EUV file that describes all of the wavelengths and
  // cross sections
  // --------------------------------------------------------------------------

  int read_file(Inputs args, Report report);

  // --------------------------------------------------------------------------
  //
  // --------------------------------------------------------------------------

  int slot_euv(std::string item,
               std::string item2,
               std::vector<float> &values,
               Report report);
};

#endif  // INCLUDE_EUV_H_
