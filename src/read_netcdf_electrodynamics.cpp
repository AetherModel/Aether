// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// -----------------------------------------------------------------------------
// The first few functions are general readers for one "record" of data in
// a binary file.
// (These files are a bit weird, because of history.  They are binary files
//  that used to be produced in fortran, so they have a record-length at
//  the beginning and end of a record. So, the general format is:
//  record-length data data data record-length
//  record-length data data data record-length
//  record-length data data data record-length
//  The specific format on what the data is can be found in the python code.)
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Reads in a 2d matrix of floats from a binary file.
// -----------------------------------------------------------------------------

arma_mat read_in_fmat_array(FILE *infile, int nRows, int nCols) {

  arma_mat values;
  float dummy;
  int reclen;
  values.set_size(nCols, nRows);
  int64_t iDummy;

  iDummy = fread(&reclen, sizeof(reclen), 1, infile);

  for (int i = 0; i < nCols; i++) {
    for (int j = 0; j < nRows; j++) {
      iDummy = fread(&dummy, sizeof(dummy), 1, infile);
      values(i, j) = dummy;
    }
  }

  iDummy = fread(&reclen, sizeof(reclen), 1, infile);

  return values;
}

// -----------------------------------------------------------------------------
// Reads in a 1d vector of floats from a binary file.
// -----------------------------------------------------------------------------

std::vector<float> read_in_float_array(FILE *infile, int nPoints) {

  std::vector<float> values;
  float dummy;
  int reclen;
  int64_t iDummy;
  iDummy = fread(&reclen, sizeof(reclen), 1, infile);

  for (int i = 0; i < nPoints; i++) {
    iDummy = fread(&dummy, sizeof(dummy), 1, infile);
    values.push_back(dummy);
  }

  iDummy = fread(&reclen, sizeof(reclen), 1, infile);

  return values;
}

// -----------------------------------------------------------------------------
// Read in a 1d vector of integers from a binary file.
// -----------------------------------------------------------------------------

std::vector<int> read_in_int_array(FILE *infile, int nPoints) {

  std::vector<int> values;
  int dummy;
  int reclen;
  int64_t iDummy;

  iDummy = fread(&reclen, sizeof(reclen), 1, infile);

  for (int i = 0; i < nPoints; i++) {
    iDummy = fread(&dummy, sizeof(dummy), 1, infile);
    values.push_back(dummy);
  }

  iDummy = fread(&reclen, sizeof(reclen), 1, infile);

  return values;
}

// -----------------------------------------------------------------------------
// Read in a string from a binary file
// -----------------------------------------------------------------------------

std::string read_in_string(FILE *infile) {

  std::string value;
  char dummy;
  int reclen;
  int64_t iDummy;

  iDummy = fread(&reclen, sizeof(reclen), 1, infile);

  for (int i = 0; i < reclen; i++) {
    iDummy = fread(&dummy, sizeof(dummy), 1, infile);
    value.push_back(dummy);
  }

  iDummy = fread(&reclen, sizeof(reclen), 1, infile);

  return value;
}

// -----------------------------------------------------------------------------
// Read in an amie-type of electrodynamics file. There is code in the
// aetherpy distribution that will produce these types of files, with
// a write function.
// -----------------------------------------------------------------------------

void Electrodynamics::read_netcdf_electrodynamics_file(std::string filename) {

  std::string function = "Electrodynamics::read_netcdf_electrodynamics_file";
  static int iFunction = -1;
  report.enter(function, iFunction);

  if (filename == "") {
    report.print(1, "No Electrodynamics File Specified");
    report.exit(function);
    return;
  }

  int reclen;
  int nLats, nMlts, nTimes, nVars;
  std::vector<double> real_times;
  arma_vec mlats_struct, mlts_struct;
  std::vector<arma_mat> potential_struct;
  std::vector<arma_mat> energy_flux_struct, average_energy_struct;
  std::vector<arma_mat> ion_energy_flux_struct, ion_average_energy_struct;

  report.print(0, "Reading Electrodynamics file : " + filename);
  FILE *infile;

  char* char_arr;
  char_arr = &filename[0];
  infile = fopen(char_arr, "rb");

  if (infile == NULL)
    std::cout << "can't find the electrodynamice file!\n";

  else {

    int64_t iDummy;
    // read in number of lats, mlts, times:
    iDummy = fread(&reclen, sizeof(reclen), 1, infile);
    iDummy = fread(&nLats, sizeof(nLats), 1, infile);
    iDummy = fread(&nMlts, sizeof(nMlts), 1, infile);
    iDummy = fread(&nTimes, sizeof(nTimes), 1, infile);
    iDummy = fread(&reclen, sizeof(reclen), 1, infile);

    std::vector<float> lats, mlts;
    lats = read_in_float_array(infile, nLats);

    for (int iLat = 0; iLat < nLats; iLat++)
      lats[iLat] = 90.0 - lats[iLat];

    mlts = read_in_float_array(infile, nMlts);

    mlats_struct = conv_to<arma_mat>::from(lats);
    mlts_struct = conv_to<arma_mat>::from(mlts);

    // read in number of variables
    iDummy = fread(&reclen, sizeof(reclen), 1, infile);
    iDummy = fread(&nVars, sizeof(nVars), 1, infile);
    iDummy = fread(&reclen, sizeof(reclen), 1, infile);

    std::vector<std::string> Vars;

    for (int i = 0; i < nVars; i++) {
      Vars.push_back(read_in_string(infile));
      report.print(2, "Reading Var : " + Vars[i]);
    }

    std::vector<int> itime;
    std::vector<float> indices;

    arma_mat values;
    std::vector<arma_mat> values_one_time;
    std::vector<std::vector<arma_mat>> all_values;

    for (int it = 0; it < nTimes; it++) {

      itime = read_in_int_array(infile, 6);
      itime.erase(itime.begin());
      itime.push_back(0);
      itime.push_back(0);
      double real_time = time_int_to_real(itime);
      real_times.push_back(real_time);

      indices = read_in_float_array(infile, 13);

      for (int iv = 0; iv < nVars; iv++) {
        values = read_in_fmat_array(infile, nMlts, nLats);

        if (it == 0)
          values_one_time.push_back(values);

        else
          values_one_time[iv] = values;
      }

      potential_struct.push_back(values_one_time[0]);
      energy_flux_struct.push_back(values_one_time[1]);
      average_energy_struct.push_back(values_one_time[2]);

      if (nVars == 5) {
        ion_energy_flux_struct.push_back(values_one_time[3]);
        ion_average_energy_struct.push_back(values_one_time[4]);
      }

      all_values.push_back(values_one_time);
    }
  }

  fclose(infile);

  HaveElectrodynamicsFile = true;

  input_electrodynamics_struct obj = input_electrodynamics_struct();
  obj.nLats = nLats;
  obj.nMlts = nMlts;
  obj.times = real_times;
  obj.mlats = mlats_struct;
  obj.mlts = mlts_struct;

  obj.potential = potential_struct;
  obj.energy_flux = energy_flux_struct;
  obj.average_energy = average_energy_struct;
  obj.ion_energy_flux = ion_energy_flux_struct;
  obj.ion_average_energy = ion_average_energy_struct;

  if (nVars == 5)
    obj.DoesIncludeIonPrecip = 1;

  input_electrodynamics.push_back(obj);

  report.exit(function);
  return;
}
