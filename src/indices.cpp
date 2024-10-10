// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "aether.h"

// ----------------------------------------------------------------------
// Initialize the Indices class
// ----------------------------------------------------------------------

Indices::Indices() {

  // Initialize the all_indices_arrays storage structure:

  int iIndex;
  index_time_pair single_index;
  single_index.nValues = 0;
  single_index.name = "";

  std::string lookup_file = input.get_indices_lookup_file();
  indices_lookup = read_json(lookup_file);

  // This is a bit wonky, but we are going to assign names to the indices
  // from the indices lookup file.
  // First, figure out how many indices are included in the file:
  nIndices = 0;

  for (auto it = indices_lookup.begin(); it != indices_lookup.end(); ++it) {
    if (it.value() > nIndices)
      nIndices = it.value();
  }

  nIndices++;

  for (iIndex = 0; iIndex < nIndices; iIndex++) {
    all_indices_arrays.push_back(single_index);

    for (auto it = indices_lookup.begin(); it != indices_lookup.end(); ++it) {
      if (it.value() == iIndex) {
        if (all_indices_arrays[iIndex].name.length() < 1)
          all_indices_arrays[iIndex].name = it.key();
      }
    }
  }
}

// ----------------------------------------------------------------------
// Reads a bunch of different files and then stores all of the data
// into the Indices class.
// Files supported:
//   - NOAA F10.7 files
//   - OMNIWeb files
// ----------------------------------------------------------------------

bool read_and_store_indices(Indices &indices) {

  bool DidWork = true;
  std::string function = "read_and_store_indices";
  static int iFunction = -1;
  report.enter(function, iFunction);
  // ---------------------------------------------------
  // Read F10.7 file (if set):
  // ---------------------------------------------------

  std::string f107_file = input.get_f107_file();

  if (f107_file.length() > 0) {
    report.print(1, "Reading F107 File : " + f107_file);
    index_file_output_struct f107_contents;
    f107_contents = read_f107_file(f107_file, indices);

    if (f107_contents.nTimes > 0)
      indices.set_f107(f107_contents);
    else {
      DidWork = false;
      report.error("ERROR in reading f107 file!!!\n");
      return DidWork;
    }
  }

  // ---------------------------------------------------
  // Read in OMNIWeb files.
  // The user can enter as many as they would like:
  // ---------------------------------------------------

  int nFiles = input.get_number_of_omniweb_files();

  if (nFiles > 0) {
    std::vector<std::string> omniweb_files = input.get_omniweb_files();

    for (int iFile = 0; iFile < nFiles; iFile++) {
      report.print(1, "Reading OMNIWEB File : " + omniweb_files[iFile]);

      index_file_output_struct file_contents;
      file_contents = read_omni_file(omniweb_files[iFile], indices);

      if (report.test_verbose(3))
        print_index_file_output_struct(file_contents);

      int nVars = file_contents.nVars;

      for (int iVar = 0; iVar < nVars ; iVar++) {
        if (file_contents.nTimes > 0) {
          DidWork = indices.set_index(file_contents.var_names[iVar],
                                      file_contents.times,
                                      file_contents.values[iVar],
                                      file_contents.missing_values[iVar]);

          if (!DidWork) {
            report.error("Error setting indices index!!!");
            return DidWork;
          }
        }  // if
      }  // for iVar
    }  // for iFile
  }  // if nFiles

  report.exit(function);
  return DidWork;
}

// ----------------------------------------------------------------------
// Perturb the indices that the user requested
// ----------------------------------------------------------------------

bool Indices::perturb() {
  bool DidWork = true;
  bool DoReport = false;
  int64_t iDebug = 2;

  json perturb_values = input.get_perturb_values();

  if (!perturb_values.empty()) {
    // User has entered some perturb values
    for (auto it = perturb_values.begin(); it != perturb_values.end(); ++it) {
      std::string name = it.key();

      if (name != "Chemistry") {

        if (report.test_verbose(iDebug)) {
          std::cout << "Perturbing Index : " << name << "\n";
          DoReport = true;
        }

        int iIndex = lookup_index_id(name);

        if (iIndex > -1) {
          int seed = input.get_updated_seed();

          if (report.test_verbose(iDebug))
            std::cout << "Index found: " << iIndex
                      << "; seed : " << seed << "\n";

          perturb_index(iIndex, seed, it.value(), DoReport);
        }
      }
    }
  }

  return DidWork;
}

// ----------------------------------------------------------------------
// Perturb a specific index in the way the user requested
// ----------------------------------------------------------------------

void Indices::perturb_index(int iIndex, int seed,
                            json style, bool DoReport) {

  int64_t nValues = all_indices_arrays[iIndex].nValues;
  int64_t nV = nValues;
  precision_t mean = 0.0;
  precision_t std;
  bool add = true;
  bool constant = false;

  if (style.contains("Mean"))
    mean = style["Mean"];

  if (style.contains("Std"))
    std = style["Std"];
  else
    std = standard_deviation(all_indices_arrays[iIndex].values);

  // Add or Multiply the random values
  if (style.contains("Add"))
    add = style["Add"];

  // Only one value for all elements or individual values for elements
  if (style.contains("Constant"))
    constant = style["Constant"];

  if (constant)
    nV = 1;

  std::vector<double> perturbations = get_normal_random_vect(mean,
                                                             std,
                                                             nV,
                                                             seed);
  int64_t iV = 0;

  for (int64_t iValue = 0; iValue < nValues; iValue++) {
    if (!constant)
      iV = iValue;

    if (add) {
      if (DoReport && iValue == 0)
        std::cout << "  ==> Adding " << perturbations[iV] << "\n";

      all_indices_arrays[iIndex].values[iValue] += perturbations[iV];
    } else {
      if (DoReport && iValue == 0)
        std::cout << "  ==> Multiplied by " << perturbations[iV] << "\n";

      all_indices_arrays[iIndex].values[iValue] *= perturbations[iV];
    }
  }
}

// ----------------------------------------------------------------------
// For f10.7 - need to set the 81-day average also.
// ----------------------------------------------------------------------

void Indices::set_f107(index_file_output_struct f107_contents) {

  // The f107 array we can just file away:

  set_index(iF107_,
            f107_contents.times,
            f107_contents.values[0],
            f107_contents.missing_values[0]);

  // We want to set the 81-day average.  This is somewhat complicated,
  // since it seems like the f107 file does have exactly 24 hour
  // spaced data.  It seems like there are 3 points per day.  Let's
  // just ignore this fact.

  // Let's simply start at the start time and then progress forward 24
  // hours at a time

  double currenttime = f107_contents.times[0];
  int64_t nTimes = f107_contents.nTimes, itime = 0;
  double eightone = 81.0 * 86400.0;

  std::vector<double> average_time;
  std::vector<float> average_f107;

  double sumf107 = 0, sumtime = 0;
  int64_t isub, nSubs;

  while (currenttime < f107_contents.times[nTimes - 1] - eightone) {

    sumf107 = 0.0;
    sumtime = 0.0;
    nSubs = 0;
    isub = itime;

    while (f107_contents.times[isub] < currenttime + eightone &&
           f107_contents.values[0][isub] != f107_contents.missing_values[0]) {
      sumf107 += f107_contents.values[0][isub];
      sumtime += f107_contents.times[isub];
      isub++;
      nSubs++;
    }

    average_time.push_back(sumtime / nSubs);
    average_f107.push_back(sumf107 / nSubs);

    itime++;
    currenttime = f107_contents.times[itime];
  }

  // Let's hold the last 81 days constant, which means we just put one
  // single value at the end which is equal to the last average value:

  average_time.push_back(f107_contents.times[nTimes - 1]);
  average_f107.push_back(sumf107 / nSubs);

  set_index(iF107A_,
            average_time,
            average_f107,
            f107_contents.missing_values[0]);
  return;
}

// ----------------------------------------------------------------------
// get functions for getting an index at a specific time
// ----------------------------------------------------------------------

precision_t Indices:: get_f107(double time) {
  return get_index(time, iF107_);
}

precision_t Indices:: get_f107a(double time) {
  return get_index(time, iF107A_);
}

// ----------------------------------------------------------------------
// This is the general function for getting an index
// ----------------------------------------------------------------------

precision_t Indices::get_index(double time, int index) {

  int64_t iLow, iMid, iHigh;

  if (all_indices_arrays[index].nValues <= 0)
    return -1.0e32;

  iLow = 0;
  iHigh = all_indices_arrays[index].nValues - 1;
  iMid = (iHigh + iLow) / 2;

  while (iHigh - iLow > 1) {

    // Break if iMid <= time < iMid+1
    if (all_indices_arrays[index].times[iMid] == time)
      break;

    if (all_indices_arrays[index].times[iMid] <= time &&
        all_indices_arrays[index].times[iMid + 1] > time)
      break;

    // Upper Half:
    if (all_indices_arrays[index].times[iMid] < time) {
      iLow = iMid;
      iMid = (iHigh + iLow) / 2;
    } else {
      iHigh = iMid;
      iMid = (iHigh + iLow) / 2;
    }
  }

  // At this point, time should be between iMid and iMid+1:

  double dt = (all_indices_arrays[index].times[iMid + 1] -
               all_indices_arrays[index].times[iMid]);
  precision_t x = (time - all_indices_arrays[index].times[iMid]) / dt;

  precision_t value = (1.0 - x) * all_indices_arrays[index].values[iMid] +
                      x * all_indices_arrays[index].values[iMid + 1];

  return value;
}

// ----------------------------------------------------------------------
// This function takes a time array and index aray and combines them
// to link them together using the index_name to identify the index
// ----------------------------------------------------------------------

bool Indices::set_index(std::string index_name,
                        std::vector<double> timearray,
                        std::vector<float> indexarray,
                        precision_t missing) {
  bool DidWork = true;
  int id = lookup_index_id(index_name);

  if (id >= 0)
    DidWork = set_index(id, timearray, indexarray, missing);
  else {
    std::cout << "Attempting to set index " << index_name
              << " but can't locate it.  Skipping.\n";
    std::cout << "  -> Modify the file indices_lookup.json file";
    std::cout << " in the UA/inputs directory\n";
    DidWork = false;
  }

  return DidWork;
}

// ----------------------------------------------------------------------
// This function takes a time array and index aray and combines them
// to link them together
// ----------------------------------------------------------------------

bool Indices::set_index(int index,
                        std::vector<double> timearray,
                        std::vector<float> indexarray,
                        precision_t missing) {

  bool DidWork = true;

  if (timearray.size() != indexarray.size()) {
    std::cout << "In set_index. Size of time and index arrays don't match!\n";
    std::cout << "  timearray : " << timearray.size() << "\n";
    std::cout << "  indexarray : " << indexarray.size() << "\n";
    DidWork = false;
  } else {
    int64_t iSize = timearray.size();
    all_indices_arrays[index].nValues = 0;

    for (int64_t i = 0; i < iSize; i++) {
      if (indexarray[i] != missing) {
        all_indices_arrays[index].times.push_back(timearray[i]);
        all_indices_arrays[index].values.push_back(indexarray[i]);
        all_indices_arrays[index].nValues++;
      }
    }
  }

  return DidWork;
}

// ----------------------------------------------------------------------
// Dump the contents of an index_file_output_struct
// ----------------------------------------------------------------------

void print_index_file_output_struct(index_file_output_struct
                                    contents) {

  std::vector<int> iTime(7, 0);
  int iVar;

  std::cout << "Outputing content of index_file_output_struct : \n";
  std::cout << "nVars : " << contents.nVars << "\n";

  for (int iVar = 0; iVar < contents.nVars; iVar++)
    std::cout << "  var_names[" << iVar << "] : "
              << contents.var_names[iVar] << "\n";

  std::cout << "nTimes : " << contents.nTimes << "\n";

  for (int64_t iLine = 0; iLine < 3; iLine++) {
    iTime = time_real_to_int(contents.times[iLine]);
    std::cout << "iLine " << iLine << ": ";
    display_itime(iTime);

    for (iVar = 0; iVar < contents.nVars; iVar++)
      std::cout << " " << contents.values[iLine][iVar];

    std::cout << "\n";
  }

  std::cout << "...\n...\n...\n";

  for (int64_t iLine = contents.nTimes - 3; iLine < contents.nTimes; iLine++) {
    iTime = time_real_to_int(contents.times[iLine]);
    std::cout << "iLine " << iLine << ": ";
    display_itime(iTime);

    for (iVar = 0; iVar < contents.nVars; iVar++)
      std::cout << " " << contents.values[iVar][iLine];

    std::cout << "\n";
  }

}

// ----------------------------------------------------------------------
// Match the name of an index to and index number
// ----------------------------------------------------------------------

int Indices::lookup_index_id(std::string name) {
  std::string name_lower = make_lower(name);
  std::string name_strip = strip_spaces(name_lower);
  int ind = -1;

  if (indices_lookup.contains(name_strip))
    ind = indices_lookup[name_strip];

  else {
    if (iProc == 0) {
      std::cout << "----------------------------------------------------\n";
      std::cout << "  -> Error !!!\n";
      std::cout << "  -> Attempting to set index " << name
                << " but can't locate it.  Skipping.\n";
      std::cout << "  -> Modify the file indices_lookup.json file";
      std::cout << " in the UA/inputs directory\n";
      std::cout << "----------------------------------------------------\n";
    }
  }

  return ind;
}

// ----------------------------------------------------------------------
// These return the index ids for various indices, allowing the user
// to pair up which index goes with which variable
// ----------------------------------------------------------------------

int Indices::get_f107_index_id() {
  return iF107_;
}
int Indices::get_f107a_index_id() {
  return iF107A_;
}
int Indices::get_imfbx_index_id() {
  return iIMFBX_;
}
int Indices::get_imfby_index_id() {
  return iIMFBY_;
}
int Indices::get_imfbz_index_id() {
  return iIMFBZ_;
}
int Indices::get_swvx_index_id() {
  return iSWVX_;
}
int Indices::get_swvy_index_id() {
  return iSWVY_;
}
int Indices::get_swvz_index_id() {
  return iSWVZ_;
}
int Indices::get_swn_index_id() {
  return iSWN_;
}
int Indices::get_swt_index_id() {
  return iSWT_;
}
int Indices::get_ae_index_id() {
  return iAE_;
}
int Indices::get_au_index_id() {
  return iAU_;
}
int Indices::get_al_index_id() {
  return iAL_;
}

std::string Indices::get_name(int iIndex) {
  if (iIndex < 0 || iIndex >= all_indices_arrays.size())
    return std::string();
  else
    return all_indices_arrays[iIndex].name;
}

int Indices::all_indices_array_size() {
  return all_indices_arrays.size();
}
