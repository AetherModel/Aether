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

Euv::Euv(Inputs args, Report report) {

  precision_t ave;

  IsOk = true;

  // Read in the EUV file:
  IsOk = read_file(args, report);

  if (IsOk) {
    // Slot the short and long wavelengths into their arrays:
    IsOk = slot_euv("Long", "", wavelengths_long, report);

    if (IsOk)
      IsOk = slot_euv("Short", "", wavelengths_short, report);

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

    // Slot the EUVAC model coefficients:
    if (args.get_euv_model() == "euvac") {
      IsOk = slot_euv("F74113", "", euvac_f74113, report);
      IsOk = slot_euv("AFAC", "", euvac_afac, report);
    }

    // Slot the EUVAC model coefficients:
    if (args.get_euv_model() == "neuvac") {
      IsOk = slot_euv("NEUV_S1", "", neuvac_s1, report);
      if (IsOk)
	IsOk = slot_euv("NEUV_S2", "", neuvac_s2, report);
      if (IsOk)
	IsOk = slot_euv("NEUV_S3", "", neuvac_s3, report);
      if (IsOk)
	IsOk = slot_euv("NEUV_P1", "", neuvac_p1, report);
      if (IsOk)
	IsOk = slot_euv("NEUV_P2", "", neuvac_p2, report);
      if (IsOk)
	IsOk = slot_euv("NEUV_I1", "", neuvac_int, report);
    }
    
  }

  IsOk = sync_across_all_procs(IsOk);
}

// ---------------------------------------------------------------------------
// Read in the EUV file that describes all of the wavelengths and
// cross sections
// ---------------------------------------------------------------------------

bool Euv::read_file(Inputs args, Report report) {

  waveinfotype tmp;
  std::string line, col;
  precision_t mulfac;
  std::ifstream infile_ptr;
  bool DidWork = true;

  report.print(1, "Reading EUV File : " + args.get_euv_file());

  infile_ptr.open(args.get_euv_file());

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

// ---------------------------------------------------------------------------
// Match rows in EUV file to different types of things, such as cross
// sections and spectra
// ---------------------------------------------------------------------------

bool Euv::slot_euv(std::string item,
                   std::string item2,
                   std::vector<float> &values,
                   Report report) {

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
		   Ions ions,
		   Inputs input,
		   Report report) {

  bool DidWork = true;

  bool includePhotoelectrons = input.get_include_photoelectrons();
  
  if (report.test_verbose(3))
    std::cout << "Euv::pair_euv \n";

  for (int iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++) {

    if (report.test_verbose(5))
      std::cout << neutrals.species[iSpecies].cName << "\n";

    neutrals.species[iSpecies].iEuvAbsId_ = -1;
    neutrals.species[iSpecies].nEuvIonSpecies = 0;
    neutrals.species[iSpecies].nEuvPeiSpecies = 0;

    // Check each row to see if the first column "name" matches:
    int64_t nEuvs = waveinfo.size();

    for (int64_t iEuv = 0; iEuv < nEuvs; iEuv++) {

      if (report.test_verbose(6))
        std::cout << "  " << waveinfo[iEuv].name << "\n";

      // if this matches...
      if (neutrals.species[iSpecies].cName == waveinfo[iEuv].name) {

        // First see if we can find absorbtion:
        if (waveinfo[iEuv].type == "abs") {
          if (report.test_verbose(5))
            std::cout << "  Found absorbtion\n";

          neutrals.species[iSpecies].iEuvAbsId_ = iEuv;
        }

        // Next see if we can find ionizations:
        if (waveinfo[iEuv].type == "ion") {

          // Loop through the ions to see if names match:
          for (int iIon = 0; iIon < ions.nSpecies; iIon++) {
            if (ions.species[iIon].cName == waveinfo[iEuv].to) {
              if (report.test_verbose(5))
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

  return DidWork;
}


// --------------------------------------------------------------------------
// Scale flux (intensity) at 1 AU to distance from the sun:
// --------------------------------------------------------------------------

int Euv::scale_from_1au(Planets planet,
                        Times time,
                        Report report) {
  int iErr = 0;
  precision_t d = planet.get_star_to_planet_dist(time);
  precision_t scale = 1.0 / (d * d);

  if (report.test_verbose(7))
    std::cout << "Scale from 1 AU : " << scale << "\n";

  for (int iWave = 0; iWave < nWavelengths; iWave++)
    wavelengths_intensity_top[iWave] = scale * wavelengths_intensity_1au[iWave];

  return iErr;
}

// --------------------------------------------------------------------------
// Calculate EUVAC
// --------------------------------------------------------------------------

int Euv::euvac(Times time,
               Indices indices,
               Report &report) {

  int iErr = 0;
  precision_t slope;

  std::string function = "Euv::euvac";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t f107 = indices.get_f107(time.get_current());
  precision_t f107a = indices.get_f107a(time.get_current());
  precision_t mean_f107 = (f107 + f107a) / 2.0;

  for (int iWave = 0; iWave < nWavelengths; iWave++) {
    slope = 1.0 + euvac_afac[iWave] * (mean_f107 - 80.0);

    if (slope < 0.8)
      slope = 0.8;

    wavelengths_intensity_1au[iWave] = euvac_f74113[iWave] * slope * pcm2topm2;
  }

  if (report.test_verbose(4)) {
    std::cout << "EUVAC output : "
              << f107 << " " << f107a
              << " -> " << mean_f107 << "\n";

    for (int iWave = 0; iWave < nWavelengths; iWave++) {
      std::cout << "     " << iWave << " "
                << wavelengths_short[iWave] << " "
                << wavelengths_long[iWave] << " "
                << wavelengths_intensity_1au[iWave] << "\n";
    }
  }

  report.exit(function);
  return iErr;
}

// --------------------------------------------------------------------------
// Calculate EUVAC
// --------------------------------------------------------------------------

int Euv::neuvac(Times time,
		Indices indices,
		Report &report) {

  int iErr = 0;
  precision_t slope;

  std::string function = "Euv::neuvac";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t f107 = indices.get_f107(time.get_current());
  precision_t f107a = indices.get_f107a(time.get_current());
  precision_t f107_diff = f107a - f107;

  precision_t f107p, f107ap;
  
  for (int iWave = 0; iWave < nWavelengths; iWave++)
    wavelengths_intensity_1au[iWave] = (
      neuvac_s1[iWave] * pow(f107, neuvac_p1[iWave]) +
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
  return iErr;
}

// --------------------------------------------------------------------------
// check to see if class is ok
// --------------------------------------------------------------------------

bool Euv::is_ok() {
  return IsOk;
}
