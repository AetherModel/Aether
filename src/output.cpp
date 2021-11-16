// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <netcdf>

#include "aether.h"

using namespace netCDF;
using namespace netCDF::exceptions;

//----------------------------------------------------------------------
// This is at the top so it doesn't have to go in the include file!
// Output a given variable to the netCDF file.  The netCDF system
// doesn't work with Armadillo cubes, so we have to transform the cube
// to a C-array, and then output
// ----------------------------------------------------------------------

void output_netcdf_3d(std::vector<size_t> count_start,
		      std::vector<size_t> count_end,
		      arma_cube value,
		      NcVar variable) {

  // Get the size of the cube:

  int64_t nX = value.n_rows;
  int64_t nY = value.n_cols;
  int64_t nZ = value.n_slices;
  int64_t iX, iY, iZ, iTotal, index;

  iTotal = nX * nY * nZ;

  // Create a temporary c-array to use to output the variable

  float *tmp_s3gc = static_cast<float*>(malloc(iTotal * sizeof(float)));

  // Move the data from the cube to the c-array

  for (iX = 0; iX < nX; iX++) {
    for (iY = 0; iY < nY; iY++) {
      for (iZ = 0; iZ < nZ; iZ++) {
        index = iX * nY * nZ + iY * nZ + iZ;
        tmp_s3gc[index] = value(iX, iY, iZ);
      }
    }
  }

  // Output the data to the netCDF file
  variable.putVar(count_start, count_end, tmp_s3gc);

  // delete the c-array
  free(tmp_s3gc);
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
//
// -----------------------------------------------------------------------------

void OutputContainer::store_variable(std::string name,
				     std::string unit,
				     arma_cube value) {
  var_struct single;
  single.cName = name;
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
  output_type = netcdf_type;
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
  std::cout << "Writing file : " << whole_filename << "\n";
  
  try{
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
  
// -----------------------------------------------------------------------------
// dump the contents of the container out into a binary file
// -----------------------------------------------------------------------------

int OutputContainer::write_container_binary() {

  int iErr = 0;
  std::ofstream binary;
  std::string whole_filename = directory + "/" + filename + ".bin";

  try{
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
// dump the contents of the container out into a binary file
// -----------------------------------------------------------------------------

int OutputContainer::write_container_netcdf() {

  int iErr = 0;
  std::string whole_filename = directory + "/" + filename + ".nc";
  std::string UNITS = "units";

  try{
    std::cout << "Writing File : " << whole_filename << "\n";
    NcFile ncdf_file(whole_filename, NcFile::replace);
    // Add dimensions:
    NcDim xDim = ncdf_file.addDim("x", elements[0].value.n_rows);
    NcDim yDim = ncdf_file.addDim("y", elements[0].value.n_cols);
    NcDim zDim = ncdf_file.addDim("z", elements[0].value.n_slices);
    NcDim tDim = ncdf_file.addDim("Time", 1);

    // Define the netCDF variables for the 3D data.
    // First create a vector of dimensions:

    std::vector<NcDim> dimVector{xDim, yDim, zDim};
    std::vector<size_t> startp{ 0, 0, 0};
    std::vector<size_t> countp{elements[0].value.n_rows,
			       elements[0].value.n_cols,
			       elements[0].value.n_slices};

    // Output time:
    NcVar timeVar = ncdf_file.addVar("Time", ncDouble, tDim);
    double time_array[1];
    time_array[0] = time_int_to_real(itime);
    timeVar.putVar(time_array);
    
    // Output all objects in the container:
    std::vector<NcVar> Var;
    int64_t nVars = elements.size();
    
    for (int64_t iVar = 0; iVar < nVars; iVar++) {
      Var.push_back(ncdf_file.addVar(elements[iVar].cName, ncFloat, dimVector));
      Var[iVar].putAtt(UNITS, elements[iVar].cUnit);
      output_netcdf_3d(startp, countp, elements[iVar].value, Var[iVar]);
    }

    ncdf_file.close();
  } catch (...) {
    std::cout << "Error writing header file : "
	      << whole_filename << "\n";
    iErr = 1;
  }
  return iErr;
}

// -----------------------------------------------------------------------------
//  Fills output containers and outputs them for common output types
// -----------------------------------------------------------------------------

int output(Neutrals neutrals,
	   Ions ions,
	   Grid grid,
	   Times time,
	   Planets planet,
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
    OutputContainer DummyOutputContainer;
    std::string output_dir = "UA/output/";
    DummyOutputContainer.set_netcdf();
    DummyOutputContainer.set_directory(output_dir);
    DummyOutputContainer.set_version(0.1);
    for (int iOutput = 0; iOutput < nOutputs; iOutput++)
      AllOutputContainers.push_back(DummyOutputContainer);
    IsFirstTime = false;
  }

  for (int iOutput = 0; iOutput < nOutputs; iOutput++) {

    if (time.check_time_gate(args.get_dt_output(iOutput))) {

      // Store time in all of the files:
      AllOutputContainers[iOutput].set_time(time.get_current());

      // Put Lon, Lat, Alt into all files:
      AllOutputContainers[iOutput].
	store_variable("Longitude", "(radians)", grid.geoLon_scgc);
      AllOutputContainers[iOutput].
	store_variable("Latitude", "(radians)", grid.geoLat_scgc);
      AllOutputContainers[iOutput].
	store_variable("Altitude", "(m)", grid.geoAlt_scgc);
      
      std::string type_output = args.get_type_output(iOutput);

      // Put certain variables into each file type

      // Neutral Densities:
      if (type_output == "neutrals" ||
	  type_output == "states") 
	for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) 
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
	for (int iSpecies = 0; iSpecies < nIons + 1; iSpecies++) 
	  AllOutputContainers[iOutput].
	    store_variable(ions.species[iSpecies].cName,
			   ions.density_unit,
			   ions.species[iSpecies].density_scgc);

      // Bulk Ion Drifts:
      if (type_output == "states") 
	for (int iDir = 0; iDir < 3; iDir++) 
	  AllOutputContainers[iOutput].
	    store_variable("Bulk" + ions.velocity_name[iDir],
			   ions.velocity_unit,
			   ions.velocity_vcgc[iDir]);

      // Electric Potential:
      if (type_output == "ions" ||
	  type_output == "states") 
	  AllOutputContainers[iOutput].
	    store_variable(ions.potential_name,
			   ions.potential_unit,
			   ions.potential_scgc);

      if (type_output == "bfield") {
	AllOutputContainers[iOutput].
	  store_variable("Magnetic Latitude", "radians", grid.magLat_scgc);
	AllOutputContainers[iOutput].
	  store_variable("Magnetic Longitude", "radians", grid.magLon_scgc);
	AllOutputContainers[iOutput].
	  store_variable("Magnetic Local Time",
			 "hours",
			 grid.magLocalTime_scgc);
	AllOutputContainers[iOutput].
	  store_variable("Beast", "nT", grid.bfield_vcgc[0]);
	AllOutputContainers[iOutput].
	  store_variable("Bnorth", "nT", grid.bfield_vcgc[1]);
	AllOutputContainers[iOutput].
	  store_variable("Bvertical", "nT", grid.bfield_vcgc[2]);
      }
      
      if (type_output == "neutrals")
	AllOutputContainers[iOutput].set_filename("3DNEU_"+time.get_YMD_HMS());
      if (type_output == "states")
	AllOutputContainers[iOutput].set_filename("3DALL_"+time.get_YMD_HMS());
      if (type_output == "ions")
	AllOutputContainers[iOutput].set_filename("3DION_"+time.get_YMD_HMS());
      if (type_output == "bfield")
	AllOutputContainers[iOutput].set_filename("3DBFI_"+time.get_YMD_HMS());
      
      AllOutputContainers[iOutput].write();
      AllOutputContainers[iOutput].clear_variables();
    }
  }
  return iErr;
}  


//----------------------------------------------------------------------
// Output a given variable to the binary file.
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
//
// -----------------------------------------------------------------------------




// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------


