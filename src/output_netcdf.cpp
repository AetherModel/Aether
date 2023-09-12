// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// Added by A. Ridley - Apr. 25, 2023

#include "aether.h"

#ifdef NETCDF

/* ---------------------------------------------------------------------

   Output and Input methods for netcdf files

 -------------------------------------------------------------------- */

#include <netcdf>

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
// read contents of a netcdf file into an output container
// -----------------------------------------------------------------------------

int OutputContainer::read_container_netcdf() {

  int iErr = 0;
  std::string whole_filename = directory + "/" + filename + ".nc";
  std::string UNITS = "units";

  try {
    std::cout << "Reading NetCDF file into container : "
              << whole_filename << "\n";
    NcFile ncdf_file_in(whole_filename, NcFile::read);

    std::multimap<std::string, NcVar> variables_in_file;
    std::multimap<std::string, NcVar>::iterator iter;

    // Declare a string to store the variable name.
    std::string variable_name;
    std::string variable_unit;
    // Declare a netCDF variable attribute object.
    NcVarAtt attribute;

    // Declare a vector of netCDF dimension objects.
    std::vector <NcDim> dimensions;
    std::string dimension_name;
    std::vector<int> nPts(3);

    // Assign the variables in the netCDF file to the multimap.
    variables_in_file = ncdf_file_in.getVars();

    // Use the iterator to loop through the multimap.
    for (iter = variables_in_file.begin();
         iter != variables_in_file.end(); iter++) {

      variable_name = iter->first;

      if (variable_name.compare("time") != 0) {

        attribute = iter->second.getAtt("units");
        attribute.getValues(variable_unit);
        dimensions = iter->second.getDims();
        int nDims =  dimensions.size();
        int iTotal = 1;

        // For this specific app, we only want the 3d arrays.
        if (nDims == 3) {

          for (int iDim = 0; iDim < nDims; iDim++) {
            dimension_name = dimensions[iDim].getName();
            nPts[iDim] = dimensions[iDim].getSize();
            iTotal = iTotal * nPts[iDim];
          }

          float *variable_array = new float[iTotal];
          iter->second.getVar(variable_array);

          arma_cube value_scgc;
          value_scgc.set_size(nPts[0], nPts[1], nPts[2]);
          int64_t index;

          // NetCDF ordering.
          for (int64_t iX = 0; iX < nPts[0]; iX++) {
            for (int64_t iY = 0; iY < nPts[1]; iY++) {
              for (int64_t iZ = 0; iZ < nPts[2]; iZ++) {
                index = iX * nPts[1] * nPts[2] + iY * nPts[2] + iZ;
                value_scgc(iX, iY, iZ) = variable_array[index];
              }
            }
          }

          // Store in the container:
          store_variable(variable_name, variable_unit, value_scgc);
        }
      } else {
        double *time_array = new double[1], t;
        iter->second.getVar(time_array);
        t = time_array[0];
        set_time(t);
      }
    }

    ncdf_file_in.close();
  } catch (...) {
    std::cout << "Error reading netcdf file : "
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
  std::string LONG_NAME = "long_name";

  try {
    NcFile ncdf_file(whole_filename, NcFile::replace);
    // Add dimensions:
    NcDim xDim = ncdf_file.addDim("x", elements[0].value.n_rows);
    NcDim yDim = ncdf_file.addDim("y", elements[0].value.n_cols);
    NcDim zDim = ncdf_file.addDim("z", elements[0].value.n_slices);
    NcDim tDim = ncdf_file.addDim("time", 1);

    // Define the netCDF variables for the 3D data.
    // First create a vector of dimensions:

    std::vector<NcDim> dimVector{xDim, yDim, zDim};
    std::vector<size_t> startp{ 0, 0, 0};
    std::vector<size_t> countp{elements[0].value.n_rows,
                               elements[0].value.n_cols,
                               elements[0].value.n_slices};

    // Output time:
    NcVar timeVar = ncdf_file.addVar("time", ncDouble, tDim);
    double time_array[1];
    time_array[0] = time_int_to_real(itime);
    timeVar.putVar(time_array);

    // Output all objects in the container:
    std::vector<NcVar> Var;
    int64_t nVars = elements.size();

    for (int64_t iVar = 0; iVar < nVars; iVar++) {
      Var.push_back(ncdf_file.addVar(elements[iVar].cName, ncFloat, dimVector));
      Var[iVar].putAtt(UNITS, elements[iVar].cUnit);

      if (elements[iVar].cLongName.length() > 0)
        Var[iVar].putAtt(LONG_NAME, elements[iVar].cLongName);

      output_netcdf_3d(startp, countp, elements[iVar].value, Var[iVar]);
    }

    ncdf_file.close();
  } catch (...) {
    std::cout << "Error writing netcdf container file : "
              << whole_filename << "\n";
    iErr = 1;
  }

  return iErr;
}

#else

/* ---------------------------------------------------------------------

   These are dummy functions for compiling without netcdf libraries.

 -------------------------------------------------------------------- */

int OutputContainer::read_container_netcdf() {
  int iErr = 1;
  std::cout << "read_container_netcdf is not working!\n";
  return iErr;
}

int OutputContainer::write_container_netcdf() {
  int iErr = 1;
  std::cout << "write_container_netcdf is not working!\n";
  return iErr;
}

#endif
