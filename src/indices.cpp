// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "aether.h"

// ----------------------------------------------------------------------
//
// ----------------------------------------------------------------------

int read_and_store_indices(Indices &indices, Inputs args, Report &report) {

  int iErr = 0;
  std::string function = "read_and_store_indices";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // ---------------------------------------------------
  // Read F10.7 file (if set):
  // ---------------------------------------------------

  std::string f107_file = args.get_f107_file();
  if (f107_file.length() > 0) {
    report.print(1,"Reading F107 File : "+f107_file);
    index_file_output_struct f107_contents;
    f107_contents = read_f107_file(f107_file, indices, report);
    if (f107_contents.nTimes > 0) {
      indices.set_f107(f107_contents);
    } else {
      iErr = 1;
      std::cout << "ERROR in reading f107 file!!!\n";
    }
  }

  // ---------------------------------------------------
  // Read in OMNIWeb files.
  // The user can enter as many as they would like:
  // ---------------------------------------------------
  
  int nFiles = args.get_number_of_omniweb_files();
  if (nFiles > 0) {
    std::vector<std::string> omniweb_files = args.get_omniweb_files();
    for (int iFile = 0; iFile < nFiles; iFile++) {
      report.print(1,"Reading OMNIWEB File : "+omniweb_files[iFile]);

      index_file_output_struct file_contents;
      file_contents = read_omni_file(omniweb_files[iFile], indices, report);
      if (report.test_verbose(3)) print_index_file_output_struct(file_contents);

      int nVars = file_contents.nVars;

      for (int iVar = 0; iVar < nVars ; iVar++) {
	if (file_contents.index_id[iVar] > -1 &&
	    file_contents.nTimes > 0) {
	  indices.set_index(file_contents.index_id[iVar],
			    file_contents.times,
			    file_contents.values[iVar],
			    file_contents.missing_values[iVar]);
	}  // if
      }  // for iVar
    }  // for iFile
  }  // if nFiles

  report.exit(function);
  return iErr;
}

// ----------------------------------------------------------------------
// Initialize
// ----------------------------------------------------------------------

Indices::Indices(Inputs args) {

  // Initialize the all_indices_arrays storage structure:
  
  int iIndex;
  index_time_pair single_index;
  single_index.nValues = 0;
  single_index.name = "Not set";
  for (iIndex = 0; iIndex < nIndices; iIndex++) {
    all_indices_arrays.push_back(single_index);
  }

}

// ----------------------------------------------------------------------
// For f10.7 - need to set the 81-day average also.
// ----------------------------------------------------------------------

void Indices::set_f107(index_file_output_struct f107_contents) {

  // The f107 array we can just file away:

  set_index(f107_,
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
  
  while (currenttime < f107_contents.times[nTimes-1]-eightone) {

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

    average_time.push_back(sumtime/nSubs);
    average_f107.push_back(sumf107/nSubs);

    itime++;
    currenttime = f107_contents.times[itime];
  }

  // Let's hold the last 81 days constant, which means we just put one
  // single value at the end which is equal to the last average value:

  average_time.push_back(f107_contents.times[nTimes-1]);
  average_f107.push_back(sumf107/nSubs);

  set_index(f107a_,
	    average_time,
	    average_f107,
	    f107_contents.missing_values[0]);
  
  return;
}

// ----------------------------------------------------------------------
// get functions for getting an index at a specific time
// ----------------------------------------------------------------------

float Indices:: get_f107(double time) {
  return get_index(time, f107_);
}

float Indices:: get_f107a(double time) {
  return get_index(time, f107a_);
}

// ----------------------------------------------------------------------
// This is the general function for getting an index
// ----------------------------------------------------------------------

float Indices::get_index(double time, int index) {

  int64_t iLow, iMid, iHigh;

  iLow = 0;
  iHigh = all_indices_arrays[index].nValues-1;
  iMid = (iHigh+iLow)/2;

  while (iHigh-iLow > 1) {

    // Break if iMid <= time < iMid+1
    if (all_indices_arrays[index].times[iMid] == time) break;
    if (all_indices_arrays[index].times[iMid] <= time &&
	all_indices_arrays[index].times[iMid+1] > time) break;
    // Upper Half:
    if (all_indices_arrays[index].times[iMid] < time) {
      iLow = iMid;
      iMid = (iHigh+iLow)/2;
    } else {
      iHigh = iMid;
      iMid = (iHigh+iLow)/2;
    }
  }

  // At this point, time should be between iMid and iMid+1:

  double dt = (all_indices_arrays[index].times[iMid+1] -
	       all_indices_arrays[index].times[iMid]);
  float x = (time - all_indices_arrays[index].times[iMid]) / dt;
  
  float value = (1.0 - x) * all_indices_arrays[index].values[iMid] +
    x * all_indices_arrays[index].values[iMid+1];

  return value;
}

// ----------------------------------------------------------------------
// This function takes a time array and index aray and combines them
// to link them together
// ----------------------------------------------------------------------

void Indices::set_index(int index,
			std::vector<double> timearray,
			std::vector<float> indexarray,
			float missing) {

  if (timearray.size() != indexarray.size()) {
    std::cout << "In set_index. Size of time and index arrays don't match!\n";
    std::cout << "  timearray : " << timearray.size() << "\n";
    std::cout << "  indexarray : " << indexarray.size() << "\n";
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

  return;
}

// ----------------------------------------------------------------------
// Dump the contents of an index_file_output_struct
// ----------------------------------------------------------------------

void print_index_file_output_struct(index_file_output_struct
				    contents) {

  std::vector<int> iTime(7,0);
  int iVar;
  
  std::cout << "Outputing content of index_file_output_struct : \n";
  std::cout << "nVars : " << contents.nVars << "\n";
  for (int iVar = 0; iVar < contents.nVars; iVar++)
    std::cout << "  var_names[" << iVar << "] : "
	      << contents.var_names[iVar] << "\n";
  std::cout << "nTimes : " << contents.nTimes << "\n";

  for (int64_t iLine = 0; iLine < 3; iLine++) {
    time_real_to_int(contents.times[iLine], iTime);
    std::cout << "iLine " << iLine << ": ";
    display_itime(iTime);
    for (iVar = 0; iVar < contents.nVars; iVar++) 
      std::cout << " " << contents.values[iLine][iVar];
    std::cout << "\n";
  }
  std::cout << "...\n...\n...\n";
  for (int64_t iLine = contents.nTimes-3; iLine < contents.nTimes; iLine++) {
    time_real_to_int(contents.times[iLine], iTime);
    std::cout << "iLine " << iLine << ": ";
    display_itime(iTime);
    for (iVar = 0; iVar < contents.nVars; iVar++) 
      std::cout << " " << contents.values[iVar][iLine];
    std::cout << "\n";
  }
  
}

int Indices::get_f107_index_id() { return f107_; }
int Indices::get_f107a_index_id() { return f107a_; }
int Indices::get_imf_bx_index_id() { return imf_bx_; }
int Indices::get_imf_by_index_id() { return imf_by_; }
int Indices::get_imf_bz_index_id() { return imf_bz_; }
int Indices::get_sw_vx_index_id() { return sw_vx_; }
int Indices::get_sw_vy_index_id() { return sw_vy_; }
int Indices::get_sw_vz_index_id() { return sw_vz_; }
int Indices::get_sw_n_index_id() { return sw_n_; }
int Indices::get_sw_t_index_id() { return sw_t_; }
int Indices::get_ae_index_id() { return ae_; }
int Indices::get_au_index_id() { return au_; }
int Indices::get_al_index_id() { return al_; }

