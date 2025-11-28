// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------------
// Initialize EUV
// -----------------------------------------------------------------------------

Euv::Euv() {

  precision_t ave;

  IsOk = true;

  if (input.get_euv_douse()) {
    doUse = true;

    // Read in the EUV file:
    IsOk = read_file();

    if (IsOk) {
      // Slot the short and long wavelengths into their arrays:
      IsOk = slot_euv("Long", "", wavelengths_long);

      if (IsOk)
        IsOk = slot_euv("Short", "", wavelengths_short);

      // This means we found both long and short wavelengths:
      if (IsOk) {
        for (int iWave = 0; iWave < nWavelengths; iWave++) {
          ave = (wavelengths_short[iWave] + wavelengths_long[iWave]) / 2.0 * cAtoM;
          wavelengths_energy.push_back(cH * cC / ave);
          // We simply want to initialize these vectors to make them the
          // correct lenght:
          wavelengths_intensity_1au.push_back(0.0);
          wavelengths_intensity_top.push_back(0.0);
        }
      }

      // Read in FISM data - does not need to be "slotted"
      if (input.get_euv_model() == "fism")
        fismData = read_fism(input.get_euv_fismfile());
      // Read in NEUVAC data - also does not need to be "slotted"
      
      // Slot the EUVAC model coefficients:
      if (input.get_euv_model() == "euvac") {
        IsOk = slot_euv("F74113", "", euvac_f74113);
        IsOk = slot_euv("AFAC", "", euvac_afac);
      }

      // Slot the NEUVAC model coefficients:
      if (input.get_euv_model() == "neuvac") {
        IsOk = slot_euv("NEUV_S1", "", neuvac_s1);

        if (IsOk)
          IsOk = slot_euv("NEUV_S2", "", neuvac_s2);

        if (IsOk)
          IsOk = slot_euv("NEUV_S3", "", neuvac_s3);

        if (IsOk)
          IsOk = slot_euv("NEUV_P1", "", neuvac_p1);

        if (IsOk)
          IsOk = slot_euv("NEUV_P2", "", neuvac_p2);

        if (IsOk)
          IsOk = slot_euv("NEUV_I1", "", neuvac_int);
      }

      // Slot the HFG model coefficients:
      if (input.get_euv_model() == "hfg") {
        IsOk = slot_euv("HFGc1", "", solomon_hfg_c1);

        if (IsOk)
          IsOk = slot_euv("HFGc2", "", solomon_hfg_c2);

        if (IsOk)
          IsOk = slot_euv("HFGfref", "", solomon_hfg_fref);
      }
    }
  } else
    doUse = false;

  IsOk = sync_across_all_procs(IsOk);
  return;
}

// ---------------------------------------------------------------------------
// Read in the EUV file that describes all of the wavelengths and
// cross sections
// ---------------------------------------------------------------------------

bool Euv::read_file() {

  waveinfotype tmp;
  std::string line, col;
  precision_t mulfac;
  std::ifstream infile_ptr;
  bool DidWork = true;

  report.print(1, "Reading EUV File : " + input.get_euv_file());

  infile_ptr.open(input.get_euv_file());

  if (!infile_ptr.is_open()) {
    if (iProc == 0)
      std::cout << "Could not open euv file!\n";

    DidWork = false;
  } else {
    nLines = 0;

    if (infile_ptr.good()) {
      int IsFirstTime = 1;

      while (getline(infile_ptr, line)) {
        report.print(5, line);
        line = strip_spaces(line);
        std::stringstream ss(line);

        // This is just to count the number of wavelengths.
        // We assume that all of the lines have the same number of wavelengths.
        if (IsFirstTime) {
          std::stringstream ssdummy(line);
          nWavelengths = 0;

          while (getline(ssdummy, col, ','))
            nWavelengths++;

          // There are 6 extra items in each line:
          nWavelengths -= 6;
        }

        getline(ss, tmp.name, ',');
        report.print(5, tmp.name);
        getline(ss, tmp.to, ',');
        getline(ss, tmp.type, ',');
        getline(ss, col, ',');
        mulfac = stof(col);
        getline(ss, tmp.units, ',');

        for (int iWavelength = 0; iWavelength < nWavelengths; iWavelength++) {
          getline(ss, col, ',');

          if (IsFirstTime)
            tmp.values.push_back(stof(col) * mulfac);
          else
            tmp.values[iWavelength] = stof(col) * mulfac;
        }

        getline(ss, tmp.note, ',');

        waveinfo.push_back(tmp);
        nLines++;
        IsFirstTime = 0;
      }

    } else
      DidWork = false;

    infile_ptr.close();
  }

  return DidWork;
}

// -------------------------------------------------------------------------------
// Read in FISM data. FISM files are created with srcPython/fism.py,
// and the data are read in to an index_file_output_struct.
// Inside the struct, we have time & each of the "variables" correspond to a
// FISM bin. This number of bins should match the number of bins in the EUV file
// -------------------------------------------------------------------------------

index_file_output_struct Euv::read_fism(std::string fism_filename) {

  std::ifstream fismfstream;
  fismfstream.open(fism_filename);
  std::vector<std::vector<std::string>>  fism_file;
  fism_file = read_csv(fismfstream);

  index_file_output_struct fism_contents;

  // one row per time
  fism_contents.nTimes = fism_file.size();
  // first six cols are the YYYY,MM,DD,HH,mm,ss (no ms)
  // the rest are the binned fism data
  fism_contents.nVars = fism_file[0].size() - 6;

  // check that the user provided the correct EUV file
  // The number of bins in euv file should match the number of fism bins ("nVars")
  if (fism_contents.nVars != nWavelengths) {
    report.error("Number of FISM wavelengths does not match the EUV file provided!");
    report.error("Either change EUV file or check your FISM file is correct.");
    IsOk = false;
  }

  std::vector<int> itime(7, 0);
  std::vector<std::vector<float>> values; // holds all values
  std::vector<float> values_tmp(fism_contents.nVars); // holds values in each row

  for (int iLine = 0; iLine < fism_file.size(); iLine ++) {

    itime[0] = stoi(fism_file[iLine][0]);
    itime[1] = stoi(fism_file[iLine][1]);
    itime[2] = stoi(fism_file[iLine][2]);
    itime[3] = stoi(fism_file[iLine][3]);
    itime[4] = stoi(fism_file[iLine][4]);
    itime[5] = stoi(fism_file[iLine][5]);
    itime[6] = 0; // 0 ms
    fism_contents.times.push_back(time_int_to_real(itime));

    for (int iVar = 0; iVar < fism_contents.nVars; iVar++)
      values_tmp[iVar] = stof(fism_file[iLine][iVar + 6]);

    values.push_back(values_tmp);
  }

  fism_contents.values = values;

  return fism_contents;
}

// ---------------------------------------------------------------------------
// Match rows in EUV file to different types of things, such as cross
// sections and spectra
// ---------------------------------------------------------------------------

bool Euv::slot_euv(std::string item,
                   std::string item2,
                   std::vector<float> &values) {

  bool DidWork = true;
  int iLine;
  int IgnoreItem2 = 0;

  report.print(3, "in slot_euv:" + item + ";" + item2);

  if (item2 == "")
    IgnoreItem2 = 1;

  // Find item to move:
  for (iLine = 0; iLine < nLines ; iLine++) {
    if (waveinfo[iLine].name == item) {
      if (IgnoreItem2)
        break;
      else if (waveinfo[iLine].to == item2)
        break;
    }
  }

  if (iLine >= nLines)
    DidWork = false;

  else {
    if (report.test_verbose(2)) {
      std::cout << "Found : " << waveinfo[iLine].name;

      if (!IgnoreItem2)
        std::cout << " with " << waveinfo[iLine].to;

      std::cout << "\n";
    }

    // Move values into the output array (values):
    for (int iWavelength = 0; iWavelength < nWavelengths; iWavelength++)
      values.push_back(waveinfo[iLine].values[iWavelength]);
  }

  return DidWork;
}

//----------------------------------------------------------------------
// This code takes the EUV information that was read in from the EUV
// file and tries to figure out which things are absorbtion/ionization
// cross sections.  It does this by comparing the name of the neutral
// species to the first column in the euv.csv file.  If it finds a
// match, it then checks to see if it is an absorbtion or ionization
// cross section.  If it is an ionization cs, then it tries to figure
// out which ion it is producing (the "to" column).
// ---------------------------------------------------------------------

bool Euv::pair_euv(Neutrals &neutrals,
                   Ions ions) {

  std::string function = "Euv::pair_euv";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool DidWork = true;

  bool includePhotoelectrons = input.get_include_photoelectrons();

  for (int iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++) {

    if (report.test_verbose(5))
      std::cout << neutrals.species[iSpecies].cName << "\n";

    neutrals.species[iSpecies].iEuvAbsId_ = -1;
    neutrals.species[iSpecies].nEuvIonSpecies = 0;
    neutrals.species[iSpecies].nEuvPeiSpecies = 0;

    // Check each row to see if the first column "name" matches:
    int64_t nEuvs = waveinfo.size();

    for (int64_t iEuv = 0; iEuv < nEuvs; iEuv++) {

      if (report.test_verbose(4))
        std::cout << "  " << waveinfo[iEuv].name << "\n";

      // if this matches...
      if (neutrals.species[iSpecies].cName == waveinfo[iEuv].name) {

        // First see if we can find absorbtion:
        if (waveinfo[iEuv].type == "abs") {
          if (report.test_verbose(4))
            std::cout << "  Found absorbtion\n";

          neutrals.species[iSpecies].iEuvAbsId_ = iEuv;
        }

        // Next see if we can find ionizations:
        if (waveinfo[iEuv].type == "ion") {

          // Loop through the ions to see if names match:
          for (int iIon = 0; iIon < ions.nSpecies; iIon++) {
            if (ions.species[iIon].cName == waveinfo[iEuv].to) {
              if (report.test_verbose(4))
                std::cout << "  Found ionization!! --> "
                          << ions.species[iIon].cName << "\n";

              neutrals.species[iSpecies].iEuvIonId_.push_back(iEuv);
              neutrals.species[iSpecies].iEuvIonSpecies_.push_back(iIon);
              neutrals.species[iSpecies].nEuvIonSpecies++;
            }  // if to
          }  // iIon loop
        }  // if ionization

        // Next see if we can find ionizations:
        if (waveinfo[iEuv].type == "pei" &&
            includePhotoelectrons) {

          // Loop through the ions to see if names match:
          for (int iIon = 0; iIon < ions.nSpecies; iIon++) {
            if (ions.species[iIon].cName == waveinfo[iEuv].to) {
              if (report.test_verbose(5))
                std::cout << "  Found photo-electron augmentation!! --> "
                          << ions.species[iIon].cName << "\n";

              neutrals.species[iSpecies].iEuvPeiId_.push_back(iEuv);
              neutrals.species[iSpecies].iEuvPeiSpecies_.push_back(iIon);
              neutrals.species[iSpecies].nEuvPeiSpecies++;
            }  // if to
          }  // iIon loop
        }  // if ionization

      }  // if species is name
    }  // for iEuv
  }  // for iSpecies

  report.exit(function);
  return DidWork;
}

// --------------------------------------------------------------------------
// Scale flux (intensity) at 1 AU to distance from the sun:
// --------------------------------------------------------------------------

void Euv::scale_from_1au(Planets planet,
                         Times time) {
  precision_t d = planet.get_star_to_planet_dist(time);
  precision_t scale = 1.0 / (d * d);

  if (report.test_verbose(7))
    std::cout << "Scale from 1 AU : " << scale << "\n";

  for (int iWave = 0; iWave < nWavelengths; iWave++)
    wavelengths_intensity_top[iWave] = scale * wavelengths_intensity_1au[iWave];

  return;
}

// --------------------------------------------------------------------------
// Calculate EUVAC
// --------------------------------------------------------------------------

bool Euv::euvac(Times time,
                Indices indices) {

  bool didWork = true;
  precision_t slope;

  std::string function = "Euv::euvac";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t f107 = indices.get_f107(time.get_current());
  precision_t f107a = indices.get_f107a(time.get_current());
  precision_t mean_f107 = (f107 + f107a) / 2.0;

  if (report.test_verbose(4))
    std::cout << "F107, f107a, average : "
              << f107 << " " << f107a
              << " -> " << mean_f107 << "\n";

  for (int iWave = 0; iWave < nWavelengths; iWave++) {
    slope = 1.0 + euvac_afac[iWave] * (mean_f107 - 80.0);

    if (slope < 0.8)
      slope = 0.8;

    wavelengths_intensity_1au[iWave] = euvac_f74113[iWave] * slope * pcm2topm2;

    if (report.test_verbose(4))
      std::cout << "     " << iWave << " "
                << wavelengths_short[iWave] << " "
                << wavelengths_long[iWave] << " "
                << wavelengths_intensity_1au[iWave] / 1e12 << " "
                << euvac_afac[iWave] * 100.0 << " "
                << euvac_f74113[iWave] / 1e9 << " "
                << slope << "\n";
  }

  report.exit(function);
  return didWork;
}

// --------------------------------------------------------------------------
// From the FISM file, interpolate the nearest 2 data to the current time
// --------------------------------------------------------------------------

bool Euv::get_fism(Times time) {
  // This is functionally similar to get_indices, however we do not store FISM in
  // the Indices class since it has variable number of bins.
  
  std::string function = "Euv::get_fism";
  static int iFunction = -1;
  report.enter(function, iFunction);

  double time_now = time.get_current();
  bool didWork = true;

  if (fism_prev_index == 0) {
    // This is probably the first time we're "running" fism.
    // Make sure the file covers the entire time range of the run.
    double end_time = time.get_end();

    if (time_now < fismData.times[0] && end_time > fismData.times[-1]) {
      report.error("FISM data does not cover the entire time range!");
      report.error("Please check that your FISM file is correct.");
      didWork = false;
    }
  }

  // Get the index prior to the current time
  while (fismData.times[fism_prev_index + 1] <= time_now)
    fism_prev_index ++;

  // Determine time-interpolation weighting factor
  precision_t dt_fism;
  dt_fism = fismData.times[fism_prev_index + 1] - fismData.times[fism_prev_index];
  precision_t x = (time_now - fismData.times[fism_prev_index]) / dt_fism;

  // store the wavelength:
  for (int iWave = 0; iWave < nWavelengths; iWave ++)
    wavelengths_intensity_1au[iWave] =
      (1.0 - x) * fismData.values[fism_prev_index][iWave]
      + x * fismData.values[fism_prev_index + 1][iWave];

  report.exit(function);
  return didWork;
}

// --------------------------------------------------------------------------
// Calculate EUVAC
// --------------------------------------------------------------------------

bool Euv::neuvac(Times time,
                 Indices indices) {

  bool didWork = true;
  precision_t slope;

  std::string function = "Euv::neuvac";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t f107 = indices.get_f107(time.get_current());
  precision_t f107a = indices.get_f107a(time.get_current());
  precision_t f107_diff = f107a - f107;

  precision_t f107p, f107ap;

  for (int iWave = 0; iWave < nWavelengths; iWave++)
    wavelengths_intensity_1au[iWave] =
      (neuvac_s1[iWave] * pow(f107, neuvac_p1[iWave]) +
       neuvac_s2[iWave] * pow(f107a, neuvac_p2[iWave]) +
       neuvac_s2[iWave] * (f107_diff) +
       neuvac_int[iWave]) / wavelengths_energy[iWave];

  if (report.test_verbose(4)) {
    std::cout << "NEUVAC output : "
              << f107 << " " << f107a
              << " -> " << f107_diff << "\n";

    for (int iWave = 0; iWave < nWavelengths; iWave++) {
      std::cout << "     " << iWave << " "
                << wavelengths_short[iWave] << " "
                << wavelengths_long[iWave] << " "
                << wavelengths_intensity_1au[iWave] << "\n";
    }
  }

  report.exit(function);
  return didWork;
}

// --------------------------------------------------------------------------
// Calculate HFG
// --------------------------------------------------------------------------

bool Euv::solomon_hfg(Times time,
                      Indices indices) {

  std::string function = "Euv::solomon_hfg";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;
  precision_t r1;
  precision_t r2;

  precision_t f107 = indices.get_f107(time.get_current());
  precision_t f107a = indices.get_f107a(time.get_current());

  for (int iWave = 0; iWave < nWavelengths; iWave++) {
    r1 = 0.0138 * (f107 - 71.5) + 0.005 * (f107 - f107a + 3.9);
    r2 = 0.5943 * (f107 - 71.5) + 0.381 * (f107 - f107a + 3.9);
    wavelengths_intensity_1au[iWave] =
      (solomon_hfg_fref[iWave] +
       (r1 * solomon_hfg_c1[iWave]) +
       (r2 * solomon_hfg_c2[iWave])) * pcm2topm2;
  }

  if (report.test_verbose(4)) {
    std::cout << "HFG output : "
              << f107 << " " << f107a << "\n";

    for (int iWave = 0; iWave < nWavelengths; iWave++) {
      std::cout << "     " << iWave << " "
                << wavelengths_short[iWave] << " "
                << wavelengths_long[iWave] << " "
                << wavelengths_intensity_1au[iWave] << "\n";
    }
  }

  report.exit(function);
  return didWork;
}

// --------------------------------------------------------------------------
// check to see if class is ok
// --------------------------------------------------------------------------

bool Euv::is_ok() {
  return IsOk;
}
