

#include <iostream>
#include <string>
#include <stdlib.h>
#include <netcdf>
#include "../include/electrodynamics.h"
#include <vector>

//g++ read_amie_bin.cpp -o -std=c++11
fmat read_in_fmat_array(FILE *infile, int nRows, int nCols) {

  // Read in latitudes
  fmat values;
  float dummy;
  int reclen;
  values.set_size(nCols, nRows);
  
  fread(&reclen, sizeof(reclen), 1, infile);
  for (int i=0; i < nCols; i++) {
    for (int j=0; j < nRows; j++) {
      fread(&dummy, sizeof(dummy), 1, infile);
      values(i,j) = dummy;
    }
  }
  fread(&reclen, sizeof(reclen), 1, infile);

  return values;
  
}

std::vector<float> read_in_float_array(FILE *infile, int nPoints) {

  // Read in latitudes
  std::vector<float> values;
  float dummy;
  int reclen;
  fread(&reclen, sizeof(reclen), 1, infile);
  for (int i=0; i < nPoints; i++) {
    fread(&dummy, sizeof(dummy), 1, infile);
    values.push_back(dummy);
    std::cout << dummy << " ";
  }
  fread(&reclen, sizeof(reclen), 1, infile);

  return values;
  
}

std::vector<int> read_in_int_array(FILE *infile, int nPoints) {

  // Read in latitudes
  std::vector<int> values;
  int dummy;
  int reclen;

  fread(&reclen, sizeof(reclen), 1, infile);
  //std::cout << "\niTime vec:\n";
  for (int i=0; i < nPoints; i++) {
    fread(&dummy, sizeof(dummy), 1, infile);
    values.push_back(dummy);
    //std::cout << dummy << " ";
  }
  fread(&reclen, sizeof(reclen), 1, infile);
  //std::cout << "\nEND\n";

  return values;
  
}

std::string read_in_string(FILE *infile) {

  // Read in latitudes
  std::string value;
  char dummy;
  int reclen;

  fread(&reclen, sizeof(reclen), 1, infile);
  for (int i=0; i < reclen; i++) {
    fread(&dummy, sizeof(dummy), 1, infile);
    value.push_back(dummy);
  }
  fread(&reclen, sizeof(reclen), 1, infile);
  std::cout << value << "\n";

  return value;
  
}

