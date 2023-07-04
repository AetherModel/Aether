// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// #include <netcdf>

#include "aether.h"
//#include "output.cpp"

/* ---------------------------------------------------------------------

   Fill output containers for certain output types.  Supported types:
   states - neutral states, ion den., bulk ion vel., temp, elec. temp
   neutrals - neutral states
   ions - Ion densites, temperatures, par & perp velocities, elec temp
   bfield - magnetic coordinates, b-field vector

 -------------------------------------------------------------------- */

// -----------------------------------------------------------------------------
//  Fills output containers and outputs them for common output types
// -----------------------------------------------------------------------------

int output(const Neutrals &neutrals,
           const Ions &ions,
           const Grid &grid,
           Times time,
           const Planets &planet,
           Inputs args,
           Report &report) {

  std::string function = "output";
  static int iFunction = -1;
  report.enter(function, iFunction);

  static bool IsFirstTime = true;

  int iErr = 0;

  int nOutputs = args.get_n_outputs();
  static std::vector<OutputContainer> AllOutputContainers;

  if (IsFirstTime) {
    // Initialize all of the output containers for all of the output
    // types requested
    OutputContainer DummyOutputContainer;
    std::string output_dir = "UA/output/";
    DummyOutputContainer.set_directory(output_dir);
    DummyOutputContainer.set_version(0.1);

    for (int iOutput = 0; iOutput < nOutputs; iOutput++)
      AllOutputContainers.push_back(DummyOutputContainer);

    IsFirstTime = false;
  }

  report.student_checker_function_name(args.get_is_student(),
				       args.get_student_name(),
				       3, "");
  
  for (int iOutput = 0; iOutput < nOutputs; iOutput++) {

    if (time.check_time_gate(args.get_dt_output(iOutput))) {

      // ------------------------------------------------------------
      // Store time in all of the files:

      AllOutputContainers[iOutput].set_time(time.get_current());

      std::string type_output = args.get_type_output(iOutput);

      // ------------------------------------------------------------
      // Put Lon, Lat, Alt into all output containers:

      if (type_output == "corners") {
        // Cell Corners:
        AllOutputContainers[iOutput].
        store_variable("lon",
                       "longitude",
                       "degrees_east",
                       grid.geoLon_Corner * cRtoD);
        AllOutputContainers[iOutput].
        store_variable("lat",
                       "latitude",
                       "degrees_north",
                       grid.geoLat_Corner * cRtoD);
        AllOutputContainers[iOutput].
        store_variable("z",
                       "height above mean sea level",
                       "m",
                       grid.geoAlt_Corner);
      } else {
        // Cell Centers:
        AllOutputContainers[iOutput].
        store_variable("lon",
                       "longitude",
                       "degrees_east",
                       grid.geoLon_scgc * cRtoD);
        AllOutputContainers[iOutput].
        store_variable("lat",
                       "latitude",
                       "degrees_north",
                       grid.geoLat_scgc * cRtoD);
        AllOutputContainers[iOutput].
        store_variable("z",
                       "height above mean sea level",
                       "m",
                       grid.geoAlt_scgc);
      }

      // ------------------------------------------------------------
      // Put certain variables into each file type

      // Neutral Densities:
      if (type_output == "neutrals" ||
          type_output == "states")
        for (int iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
          AllOutputContainers[iOutput].
          store_variable(neutrals.species[iSpecies].cName,
                         neutrals.density_unit,
                         neutrals.species[iSpecies].density_scgc);

      // Neutral Temperature:
      if (type_output == "neutrals" ||
          type_output == "states")
        AllOutputContainers[iOutput].
        store_variable(neutrals.temperature_name,
                       neutrals.temperature_unit,
                       neutrals.temperature_scgc);

      // Neutral Winds:
      if (type_output == "neutrals" ||
          type_output == "states")
        for (int iDir = 0; iDir < 3; iDir++)
          AllOutputContainers[iOutput].
          store_variable(neutrals.velocity_name[iDir],
                         neutrals.velocity_unit,
                         neutrals.velocity_vcgc[iDir]);

      // Ion Densities:
      if (type_output == "ions" ||
          type_output == "states")
        for (int iSpecies = 0; iSpecies <= ions.nSpecies; iSpecies++)
          AllOutputContainers[iOutput].
          store_variable(ions.species[iSpecies].cName,
                         ions.density_unit,
                         ions.species[iSpecies].density_scgc);

      // Ion Temperatures:
      if (type_output == "ions" ||
          type_output == "states")
        for (int iSpecies = 0; iSpecies <= ions.nSpecies; iSpecies++)
          AllOutputContainers[iOutput].
          store_variable(ions.species[iSpecies].cName + " " + ions.temperature_name,
                         ions.temperature_unit,
                         ions.species[iSpecies].temperature_scgc);

      // Bulk Ion Temperature:
      if (type_output == "ions" ||
          type_output == "states")
        AllOutputContainers[iOutput].store_variable("Bulk Ion " +
                                                    ions.temperature_name,
                                                    ions.temperature_unit,
                                                    ions.temperature_scgc);

      // Bulk Ion Drifts:
      if (type_output == "states")
        for (int iDir = 0; iDir < 3; iDir++)
          AllOutputContainers[iOutput].store_variable("Bulk " +
                                                      ions.velocity_name[iDir],
                                                      ions.velocity_unit,
                                                      ions.velocity_vcgc[iDir]);

      // Electric Potential:
      if (type_output == "ions" ||
          type_output == "states")
        AllOutputContainers[iOutput].store_variable(ions.potential_name,
                                                    ions.potential_unit,
                                                    ions.potential_scgc);

      if (type_output == "bfield") {
        AllOutputContainers[iOutput].store_variable("mlat",
                                                    "Magnetic Latitude",
                                                    "degrees",
                                                    grid.magLat_scgc * cRtoD);
        AllOutputContainers[iOutput].store_variable("mlon",
                                                    "Magnetic Longitude",
                                                    "degrees",
                                                    grid.magLon_scgc * cRtoD);
        AllOutputContainers[iOutput].store_variable("mlt",
                                                    "Magnetic Local Time",
                                                    "hours",
                                                    grid.magLocalTime_scgc);
        AllOutputContainers[iOutput].store_variable("Beast",
                                                    "nT",
                                                    grid.bfield_vcgc[0]);
        AllOutputContainers[iOutput].store_variable("Bnorth",
                                                    "nT",
                                                    grid.bfield_vcgc[1]);
        AllOutputContainers[iOutput].store_variable("Bvertical",
                                                    "nT",
                                                    grid.bfield_vcgc[2]);
      }

      if (type_output == "cent_acc") {
        AllOutputContainers[iOutput].store_variable("radius",
                                                    "Magnetic Latitude",
                                                    "degrees",
                                                    grid.radius_scgc);
        AllOutputContainers[iOutput].store_variable("geoLat",
                                                    "Magnetic Longitude",
                                                    "degrees",
                                                    grid.geoLat_scgc);
      }

      // ------------------------------------------------------------
      // Set output file names

      std::string filename;

      if (type_output == "neutrals")
        filename = "3DNEU_";

      if (type_output == "states")
        filename = "3DALL_";

      if (type_output == "ions")
        filename = "3DION_";

      if (type_output == "bfield")
        filename = "3DBFI_";
        
      if (type_output == "cent_acc")
        filename = "3DFBI_";

      if (type_output == "corners")
        filename = "3DCOR_";

      filename = filename + time.get_YMD_HMS();

      if (nMembers > 1)
        filename = filename + "_" + cMember;

     
      filename = filename + "_" + cGrid;
      

      report.print(0, "Writing file : " + filename);
      AllOutputContainers[iOutput].set_filename(filename);

      // ------------------------------------------------------------
      // write output container

      AllOutputContainers[iOutput].write();

      // ------------------------------------------------------------
      // Clear variables for next time

      AllOutputContainers[iOutput].clear_variables();
    }
  }

  report.exit(function);
  return iErr;
}


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

void OutputContainer::write() {

  int iErr = 0;

  if (output_type == binary_type) {
    iErr = write_container_header();

    if (iErr == 0)
      iErr = write_container_binary();
  }

  if (output_type == netcdf_type)
    iErr = write_container_netcdf();
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

int OutputContainer::write_container_header() {

  int iErr = 0;
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
    std::cout << "Error writing header file : "
              << whole_filename << "\n";
    iErr = 1;
  }

  return iErr;
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

int OutputContainer::write_container_binary() {

  int iErr = 0;
  std::ofstream binary;
  std::string whole_filename = directory + "/" + filename + ".bin";

  try {
    binary.open(whole_filename, ios::binary | ios::out);

    int64_t nVars = elements.size();

    for (int64_t iVar = 0; iVar < nVars; iVar++)
      output_binary_3d(binary, elements[iVar].value);

    return iErr;
  } catch (...) {
    std::cout << "Error writing header file : "
              << whole_filename << "\n";
    iErr = 1;
  }

  return iErr;
}




// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------


