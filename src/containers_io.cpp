// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

/* ---------------------------------------------------------------------

   General Methods for the output containers

 -------------------------------------------------------------------- */

// -----------------------------------------------------------------------------
// gets the number of elements in the output container
// -----------------------------------------------------------------------------

int64_t OutputContainer::get_nElements() {
  return elements.size();
}

// -----------------------------------------------------------------------------
// gets the ith element arma_cube in the output container
// -----------------------------------------------------------------------------

arma_cube OutputContainer::get_element_value(int64_t iElement) {
  arma_cube val;

  if (iElement < elements.size())
    val = elements[iElement].value;

  return val;
}

// -----------------------------------------------------------------------------
// gets the element arma_cube in the output container given name
// -----------------------------------------------------------------------------

arma_cube OutputContainer::get_element_value(std::string var_to_get) {
  arma_cube val;
  int64_t iElement = find_element(var_to_get);

  if (iElement >= 0)
    val = elements[iElement].value;

  return val;
}

// -----------------------------------------------------------------------------
// gets the ith element name in the output container
// -----------------------------------------------------------------------------

std::string OutputContainer::get_element_name(int64_t iElement) {
  std::string val = "";

  if (iElement < elements.size())
    val = elements[iElement].cName;

  return val;
}

// -----------------------------------------------------------------------------
// find the element number for the given variable name
// -----------------------------------------------------------------------------

int64_t OutputContainer::find_element(std::string var_to_find) {
  int64_t nVars = elements.size();
  int64_t iVarSave = -1;

  for (int64_t iVar = 0; iVar < nVars; iVar++)
    if (elements[iVar].cName.compare(var_to_find) == 0)
      iVarSave = iVar;

  return iVarSave;
}

// -----------------------------------------------------------------------------
// sets the output type to be binary
// -----------------------------------------------------------------------------

void OutputContainer::set_binary() {
  output_type = binary_type;
}

// -----------------------------------------------------------------------------
// sets the output type to be netcdf
// -----------------------------------------------------------------------------

void OutputContainer::set_netcdf() {
  output_type = netcdf_type;
}

// -----------------------------------------------------------------------------
// sets the output type to be hdf5
// -----------------------------------------------------------------------------

void OutputContainer::set_hdf5() {
  output_type = hdf5_type;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

void OutputContainer::set_directory(std::string in_dir) {
  directory = in_dir;
}

// -----------------------------------------------------------------------------
// Set the filename for file output.  We can increase the
// sophistication of this system by adding on the timestamp.
// -----------------------------------------------------------------------------

void OutputContainer::set_filename(std::string in_filename) {
  filename = in_filename;
}

// -----------------------------------------------------------------------------
// Store variable information into the output container
// -----------------------------------------------------------------------------

void OutputContainer::store_variable(std::string name,
                                     std::string unit,
                                     arma_cube value) {
  var_struct single;
  single.cName = name;
  single.cLongName = "";
  single.cUnit = unit;
  single.value = value;
  elements.push_back(single);
}

// -----------------------------------------------------------------------------
// Store variable information into the output container
//   - add long_name (needed for netcdf standardization)
// -----------------------------------------------------------------------------

void OutputContainer::store_variable(std::string name,
                                     std::string long_name,
                                     std::string unit,
                                     arma_cube value) {
  var_struct single;
  single.cName = name;
  single.cLongName = long_name;
  single.cUnit = unit;
  single.value = value;
  elements.push_back(single);
}

// -----------------------------------------------------------------------------
// Set the time of the output file
// -----------------------------------------------------------------------------

void OutputContainer::set_time(double time) {
  itime = time_real_to_int(time);
}

// -----------------------------------------------------------------------------
// Set the version number of the output. This, in theory, should be
// the version of Aether.
// -----------------------------------------------------------------------------

void OutputContainer::set_version(float in_version) {
  version = in_version;
}

// -----------------------------------------------------------------------------
// Set the number of ghostcells.
// -----------------------------------------------------------------------------

void OutputContainer::set_nGhostCells(int in_nGCs) {
  nGCs = in_nGCs;
}

// -----------------------------------------------------------------------------
// Clears the elements vector within the output container
// -----------------------------------------------------------------------------

void OutputContainer::clear_variables() {
  elements.clear();
}

// -----------------------------------------------------------------------------
// Initialize the output container
// -----------------------------------------------------------------------------

OutputContainer::OutputContainer() {
  // Set default output type to netCDF
  nGCs = 0;
#ifdef NETCDF
  output_type = netcdf_type;
#else
  output_type = binary_type;
#endif
}

// -----------------------------------------------------------------------------
// This is the write method. Here we look at which type of file output
// the user wants, and then output that particular type.
// -----------------------------------------------------------------------------

bool OutputContainer::write() {

  bool didWork = true;

  if (output_type == binary_type) {
    didWork = write_container_header();

    if (didWork)
      didWork = write_container_binary();
  }

  if (output_type == netcdf_type)
    didWork = write_container_netcdf();

  return didWork;
}

// -----------------------------------------------------------------------------
// This is the read method. Here we look at which type of file output
// the user wants, and then read that particular type.
// -----------------------------------------------------------------------------

bool OutputContainer::read() {

  bool didWork = true;

  if (output_type == binary_type)
    didWork = read_container_binary();

  if (output_type == netcdf_type)
    didWork = read_container_netcdf();

  return didWork;
}



// -----------------------------------------------------------------------------
// display the contents of an output container to the screen
// -----------------------------------------------------------------------------

void OutputContainer::display() {
  std::cout << "Displaying Container Information:\n";
  std::cout << "  time : ";
  display_itime(itime);
  std::cout << "  nX : " << elements[0].value.n_rows << "\n";
  std::cout << "  nY : " << elements[0].value.n_cols << "\n";
  std::cout << "  nZ : " << elements[0].value.n_slices << "\n";
  std::cout << "  nGCs : " << nGCs << "\n";
  int64_t nVars = elements.size();
  std::cout << "  Number of Variables : " << nVars << "\n";

  for (int64_t iVar = 0; iVar < nVars; iVar++) {
    std::cout << "  Variable " << iVar << ": " << elements[iVar].cName << "\n";
    std::cout << "      Unit  : " << elements[iVar].cUnit << "\n";

    if (elements[iVar].cLongName.length() > 0)
      std::cout << "      Long Name  : " << elements[iVar].cLongName;
  }
}

/* ---------------------------------------------------------------------

   Output methods for binary files

 -------------------------------------------------------------------- */

// -----------------------------------------------------------------------------
// Write a header file for the container.  This outputs into a json
// formatted file.
// -----------------------------------------------------------------------------

bool OutputContainer::write_container_header() {

  bool didWork = true;
  int64_t nVars = elements.size();
  int64_t nX = elements[0].value.n_rows;
  int64_t nY = elements[0].value.n_cols;
  int64_t nZ = elements[0].value.n_slices;

  std::vector<std::string> variables;
  std::vector<std::string> units;

  for (int64_t iVar = 0; iVar < nVars; iVar++) {
    variables.push_back(elements[iVar].cName);
    units.push_back(elements[iVar].cUnit);
  }

  json header = json::object({ {"version", version},
    {"time", itime},
    {"nVars", nVars},
    {"nX", nX},
    {"nY", nY},
    {"nZ", nZ},
    {"nGCs", nGCs},
    {"nLons", nX},
    {"nLats", nY},
    {"nAlts", nZ},
    {"variables", variables},
    {"units", units} });
  std::string whole_filename = directory + "/" + filename + ".json";

  try {
    std::ofstream file(whole_filename);
    file << header.dump(4) << "\n";
    file.close();
  } catch (...) {
    report.error("Error writing header file : " + whole_filename);
    didWork = false;
  }

  return didWork;
}

// -----------------------------------------------------------------------------
// Write a header file for the container.  This outputs into a json
// formatted file.
// -----------------------------------------------------------------------------

bool OutputContainer::read_container_binary() {

  std::string function = "OutputContainer::read_container_binary";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;
  std::string bin_filename = directory + "/" + filename + ".bin";
  std::string json_filename = directory + "/" + filename + ".json";

  json header = read_json(json_filename);

  if (report.test_verbose(2)) {
    std::cout << "reading binary restart\n --> json header is here:\n";
    std::cout << std::setw(2) << header;
  }

  int64_t iVar, nVars = header["nVars"];
  int64_t iX, nX = header["nX"];
  int64_t iY, nY = header["nY"];
  int64_t iZ, nZ = header["nZ"];
  int64_t iTotalSize = nX * nY * nZ;
  nGCs = header["nGCs"];

  float *variable_array = new float[iTotalSize];
  arma_cube value_scgc;
  value_scgc.set_size(nX, nY, nZ);
  int64_t index;

  std::ifstream binary;
  binary.open(bin_filename, ios::binary | ios::in);

  // Now, read and store variable-by-variable

  for (iVar = 0; iVar < nVars; iVar++) {

    // Read from the binary file
    binary.read((char *) variable_array, iTotalSize * sizeof(float));

    for (iZ = 0; iZ < nZ; iZ++) {
      for (iY = 0; iY < nY; iY++) {
        for (iX = 0; iX < nX; iX++) {
          // Python ordering!
          index = iX + iY * nX + iZ * nY * nX;
          value_scgc(iX, iY, iZ) = variable_array[index];
        }
      }
    }

    // Store in the container:
    if (report.test_verbose(2)) {
      std::cout << "Storing Variable : ";
      std::cout << header["variables"][iVar] << " : " << value_scgc(int(nX / 2),
                int(nY / 2), int(nZ / 2)) << "\n";
    }

    store_variable(header["variables"][iVar], header["units"][iVar], value_scgc);
  }

  report.exit(function);
  return didWork;
}

//----------------------------------------------------------------------
// Output a given arma_cube to the binary file.
// ----------------------------------------------------------------------

void output_binary_3d(std::ofstream &binary,
                      arma_cube value) {

  // Get the size of the cube:

  int64_t nX = value.n_rows;
  int64_t nY = value.n_cols;
  int64_t nZ = value.n_slices;
  int64_t iX, iY, iZ, index;

  int64_t nPts = nX * nY * nZ;
  int64_t iTotalSize = nPts * sizeof(float);

  // Create a temporary c-array to use to output the variable
  float *tmp_s3gc = static_cast<float*>(malloc(iTotalSize));

  // Move the data from the cube to the c-array

  for (iZ = 0; iZ < nZ; iZ++) {
    for (iY = 0; iY < nY; iY++) {
      for (iX = 0; iX < nX; iX++) {
        // Python ordering!
        index = iX + iY * nX + iZ * nY * nX;
        tmp_s3gc[index] = value(iX, iY, iZ);
      }
    }
  }

  // Output the data to the binary file
  binary.write((char *) tmp_s3gc, iTotalSize);

  // delete the c-array
  free(tmp_s3gc);
}

// -----------------------------------------------------------------------------
// dump the contents of the container out into a binary file
// -----------------------------------------------------------------------------

bool OutputContainer::write_container_binary() {

  bool didWork = true;
  std::ofstream binary;
  std::string whole_filename = directory + "/" + filename + ".bin";

  try {
    binary.open(whole_filename, ios::binary | ios::out);

    int64_t nVars = elements.size();

    for (int64_t iVar = 0; iVar < nVars; iVar++)
      output_binary_3d(binary, elements[iVar].value);

    return didWork;
  } catch (...) {
    report.error("Error writing binary file : " + whole_filename);
    didWork = false;
  }

  return didWork;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
