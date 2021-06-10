

#include "../include/aether.h"
#include <stdlib.h>

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
  std::cout << "\n";

  return values;
  
}

std::vector<int> read_in_int_array(FILE *infile, int nPoints) {

  // Read in latitudes
  std::vector<int> values;
  int dummy;
  int reclen;

  fread(&reclen, sizeof(reclen), 1, infile);
  for (int i=0; i < nPoints; i++) {
    fread(&dummy, sizeof(dummy), 1, infile);
    values.push_back(dummy);
    std::cout << dummy << " ";
  }
  fread(&reclen, sizeof(reclen), 1, infile);
  std::cout << "\n";

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

int main() {

  int reclen;
  int nLats, nMlts, nTimes, nVars;

  std::string filename = "test_ed.bin";

  std::cout << filename << "\n";

  FILE *infile;

  char* char_arr;
  char_arr = &filename[0];
  infile = fopen(char_arr, "rb");
  if (infile == NULL) {
    std::cout << "can't find the file!\n";
  } else {

    // read in number of lats, mlts, times:
    fread(&reclen, sizeof(reclen), 1, infile);
    fread(&nLats, sizeof(nLats), 1, infile);
    fread(&nMlts, sizeof(nMlts), 1, infile);
    fread(&nTimes, sizeof(nTimes), 1, infile);
    std::cout << nLats << " "
	      << nMlts << " "
	      << nTimes << "\n";
    fread(&reclen, sizeof(reclen), 1, infile);

    std::vector<float> lats, mlts;
    lats = read_in_float_array(infile, nLats);
    mlts = read_in_float_array(infile, nMlts);

    // read in number of variables
    fread(&reclen, sizeof(reclen), 1, infile);
    fread(&nVars, sizeof(nVars), 1, infile);
    std::cout << nVars << "\n";
    fread(&reclen, sizeof(reclen), 1, infile);

    std::vector<std::string> Vars;
    for (int i = 0; i < nVars; i++) {
      Vars.push_back(read_in_string(infile));
    }

    std::vector<int> itime;
    std::vector<float> indices;

    fmat values;
    std::vector<fmat> values_one_time;
    std::vector<std::vector<fmat>> all_values;

    for (int it = 0; it < nTimes; it++) {

      itime = read_in_int_array(infile, 6);
      indices = read_in_float_array(infile, 13);
    
      for (int iv = 0; iv < nVars; iv++) { 
	values = read_in_fmat_array(infile, nMlts, nLats);
	if (it == 0) {
	  values_one_time.push_back(values);
	} else {
	  values_one_time[iv] = values;
	}
      }
      all_values.push_back(values_one_time);
    }      
  }
  fclose(infile);
}
