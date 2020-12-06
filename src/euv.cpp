// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <vector>
#include <sstream> 
#include <iostream>

#include "../include/constants.h"
#include "../include/inputs.h"
#include "../include/indices.h"
#include "../include/euv.h"
#include "../include/report.h"

// -----------------------------------------------------------------------------
// Initialize EUV
// -----------------------------------------------------------------------------

Euv::Euv(Inputs args, Report report) {

  int iErr;
  float ave;

  iErr = 0;

  // Read in the EUV file:
  iErr = read_file(args, report);

  if (!iErr) {
    // Slot the short and long wavelengths into their arrays:
    iErr = slot_euv("Long", "", wavelengths_long, report);
    if (!iErr) iErr = slot_euv("Short", "", wavelengths_short, report);

    // This means we found both long and short wavelengths:
    if (!iErr) {
      for (int iWave = 0; iWave < nWavelengths; iWave++) {
	ave = (wavelengths_short[iWave] + wavelengths_long[iWave])/2.0;
	wavelengths_energy.push_back(planck_constant *
				     speed_light /
				     (ave * atom));
	// We simply want to initialize these vectors to make them the
	// correct lenght:
	wavelengths_intensity_1au.push_back(0.0);
	wavelengths_intensity_top.push_back(0.0);
      }
    }

    // Slot the EUVAC model coefficients:
    if (args.get_euv_model() == "euvac") {
      iErr = slot_euv("F74113", "", euvac_f74113, report);
      iErr = slot_euv("AFAC", "", euvac_afac, report);
    }
  }
}

// ---------------------------------------------------------------------------
// Read in the EUV file that describes all of the wavelengths and
// cross sections
// ---------------------------------------------------------------------------

int Euv::read_file(Inputs args, Report report) {

  waveinfotype tmp;
  std::string line, col;
  float mulfac;
  std::ifstream infile_ptr;
  int iErr = 0;

  report.print(1, "Reading EUV File : "+args.get_euv_file());

  infile_ptr.open(args.get_euv_file());

  if (!infile_ptr.is_open()) {
    std::cout << "Could not open euv file!\n";
    iErr = 1;
  } else {

    nLines = 0;

    if (infile_ptr.good()) {

      int IsFirstTime = 1;

      while (getline(infile_ptr,line)) {

	report.print(5, line);
	std::stringstream ss(line);

	// This is just to count the number of wavelengths.
	// We assume that all of the lines have the same number of wavelengths.
	if (IsFirstTime) {
	  std::stringstream ssdummy(line);
	  nWavelengths = 0;
	  while (getline(ssdummy, col, ',')) {
	    nWavelengths++;
	  }
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

	for (int iWavelength=0; iWavelength < nWavelengths; iWavelength++) {
	  getline(ss, col, ',');
	  if (IsFirstTime) tmp.values.push_back(stof(col) * mulfac);
	  else tmp.values[iWavelength] = stof(col) * mulfac;
	}
	getline(ss, tmp.note, ',');

	waveinfo.push_back(tmp);
	nLines++;
	IsFirstTime = 0;

      }

    } else {

      iErr = 1;

    }

    infile_ptr.close();

  }

  return iErr;

}

// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------

int Euv::slot_euv(std::string item,
		  std::string item2,
		  std::vector<float> &values,
		  Report report) {

  int iErr = 0;
  int iLine;
  int IgnoreItem2 = 0;

  report.print(3, "in slot_euv:" + item + ";" + item2);

  if (item2 == "") IgnoreItem2 = 1;

  // Find item to move:
  for (iLine=0; iLine < nLines ; iLine++) {
    if (waveinfo[iLine].name == item) {
      if (IgnoreItem2) break;
      else if (waveinfo[iLine].to == item2) break;
    }
  }

  if (iLine >= nLines) {
    iErr = 1;
  } else {

    if (report.test_verbose(2)) {
      std::cout << "Found : " << waveinfo[iLine].name;
      if (!IgnoreItem2) std::cout << " with " << waveinfo[iLine].to;
      std::cout << "\n";
    }

    // Move values into the output array (values):
    for (int iWavelength=0; iWavelength < nWavelengths; iWavelength++) {
      values.push_back(waveinfo[iLine].values[iWavelength]);
    }

  }

  return iErr;

}



// --------------------------------------------------------------------------
// Scale flux (intensity) at 1 AU to distance from the sun:
// --------------------------------------------------------------------------

int Euv::scale_from_1au(Planets planet,
			Times time) {

  int iErr = 0;
  float d = planet.get_star_to_planet_dist(time);
  float scale = 1.0 / (d*d);

  for (int iWave = 0; iWave < nWavelengths; iWave++)
    wavelengths_intensity_top[iWave] = scale * wavelengths_intensity_1au[iWave];

  return iErr;

}

// --------------------------------------------------------------------------
// EUVAC
// --------------------------------------------------------------------------

int Euv::euvac(Times time,
	       Indices indices,
	       Report report) {

  int iErr = 0;
  float slope;

  float f107 = indices.get_f107(time.get_current());
  float f107a = indices.get_f107a(time.get_current());
  float mean_f107 = (f107 + f107a)/2.0;

  if (report.test_verbose(5))
    std::cout << f107 << " " << f107a << "\n";

  for (int iWave = 0; iWave < nWavelengths; iWave++) {

    slope = 1.0 + euvac_afac[iWave] * (mean_f107 - 80.0);
    if (slope < 0.8) slope = 0.8;
    wavelengths_intensity_1au[iWave] = euvac_f74113[iWave] * slope * pcm2topm2;

  }

  if (report.test_verbose(5)) {

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

  return iErr;

}
