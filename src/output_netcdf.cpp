// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <netcdf>

#include "../include/aether.h"

#include <algorithm>

using namespace netCDF;
using namespace netCDF::exceptions;

//----------------------------------------------------------------------
// Output a given variable to the netCDF file.  The netCDF system
// doesn't work with Armadillo cubes, so we have to transform the cube
// to a C-array, and then output
// ----------------------------------------------------------------------

void output_variable_3d(std::vector<size_t> count_start,
                        std::vector<size_t> count_end,
                        fcube value,
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
        index = iX*nY*nZ + iY*nZ + iZ;
        tmp_s3gc[index] = value(iX, iY, iZ);
      }
    }
  }

  // Output the data to the netCDF file
  variable.putVar(count_start, count_end, tmp_s3gc);

  // delete the c-array
  free(tmp_s3gc);
}

//----------------------------------------------------------------------
// Output the different file types to netCDF files. 
//----------------------------------------------------------------------

int output(Neutrals neutrals,
           Ions ions,
           Grid grid,
           Times time,
           Planets planet,
           Inputs args,
           Report &report) {

  int iErr = 0;

  int nOutputs = args.get_n_outputs();

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  double time_array[1];

  std::string function = "output";
  static int iFunction = -1;
  report.enter(function, iFunction);

  for (int iOutput = 0; iOutput < nOutputs; iOutput++) {

    if (time.check_time_gate(args.get_dt_output(iOutput))) {

      grid.calc_sza(planet, time, report);
      grid.calc_gse(planet, time, report);
      grid.calc_mlt(report);

      std::string time_string;
      std::string file_name;
      std::string file_ext = ".nc";
      std::string UNITS = "units";
      std::string DESC = "description";
      std::string MIN_VAL = "minimum value";
      std::string MAX_VAL = "maximum value";
      std::string file_pre;

      std::string type_output = args.get_type_output(iOutput);

      if (type_output == "neutrals") file_pre = "3DNEU";
      if (type_output == "states") file_pre = "3DALL";
      if (type_output == "bfield") file_pre = "3DBFI";

      time_string = time.get_YMD_HMS();
      file_name = file_pre + "_" + time_string + file_ext;

      // Create the file:
      report.print(0, "Writing file : " + file_name);
      NcFile ncdf_file(file_name, NcFile::replace);

      // Add dimensions:
      NcDim lonDim = ncdf_file.addDim("longitude", nLons);
      NcDim latDim = ncdf_file.addDim("latitude", nLats);
      NcDim altDim = ncdf_file.addDim("altitude", nAlts);

      NcDim timeDim = ncdf_file.addDim("Time", 1);

      // Define the Coordinate Variables

      // Define the netCDF variables for the 3D data.
      // First create a vector of dimensions:

      std::vector<NcDim> dimVector;
      dimVector.push_back(lonDim);
      dimVector.push_back(latDim);
      dimVector.push_back(altDim);

      NcVar timeVar = ncdf_file.addVar("Time", ncDouble, timeDim);
      NcVar lonVar = ncdf_file.addVar("Geographic_Longitude", ncFloat, dimVector);
      NcVar latVar = ncdf_file.addVar("Geographic_Latitude", ncFloat, dimVector);
      NcVar altVar = ncdf_file.addVar("Altitude", ncFloat, dimVector);

      timeVar.putAtt(UNITS, "seconds");
      lonVar.putAtt(UNITS, "radians");
      //lonVar.putAtt(UNITS, "degrees");
      latVar.putAtt(UNITS, "radians");
      //latVar.putAtt(UNITS, "degrees");
      altVar.putAtt(UNITS, "meters");

      timeVar.putAtt(DESC, "Time in s from XXX");
      lonVar.putAtt(DESC, "Geographic longitude in radians");
      latVar.putAtt(DESC, "Geographic latitude in radians");
      altVar.putAtt(DESC, "Altitude in m");

      std::vector<size_t> startp, countp;
      startp.push_back(0);
      startp.push_back(0);
      startp.push_back(0);

      countp.push_back(nLons);
      countp.push_back(nLats);
      countp.push_back(nAlts);

      // Output time:

      time_array[0] = time.get_current();
      timeVar.putVar(time_array);


      // Output longitude, latitude, altitude 3D arrays:

      output_variable_3d(startp, countp, grid.geoLon_scgc, lonVar);
      output_variable_3d(startp, countp, grid.geoLat_scgc, latVar);
      output_variable_3d(startp, countp, grid.geoAlt_scgc, altVar);


      // ----------------------------------------------
      // Neutral Densities and Temperature
      // ----------------------------------------------

      if (type_output == "neutrals" ||
          type_output == "states") {

        // Output all species densities:
        std::vector<NcVar> denVar;
        for (int iSpecies=0; iSpecies < nSpecies; iSpecies++) {
          if (report.test_verbose(3))
            std::cout << "Outputting Var : "
                      << neutrals.species[iSpecies].cName << "\n";
          denVar.push_back(ncdf_file.addVar("Neutral_Density_" + neutrals.species[iSpecies].cName,
                                            ncFloat, dimVector));
          denVar[iSpecies].putAtt(UNITS, neutrals.density_unit);
          denVar[iSpecies].putAtt(DESC, "Specie density of " + neutrals.species[iSpecies].cName + " in " + neutrals.density_unit);

          // Work in Progress: calculating min and max values for MIN_VAL and MAX_VAL
          //int min = std::min_element(neutrals.species[iSpecies].density_scgc.begin(), neutrals.species[iSpecies].density_scgc.end());
          //denVar[iSpecies].putAtt(MIN_VAL, std::min_element(neutrals.species[iSpecies].density_scgc.begin(), neutrals.species[iSpecies].density_scgc.end()));

          output_variable_3d(startp, countp,
                             neutrals.species[iSpecies].density_scgc,
                             denVar[iSpecies]);
        }

        // Output bulk temperature:
        NcVar tempVar = ncdf_file.addVar("Neutral_" + neutrals.temperature_name,
                                         ncFloat, dimVector);
        tempVar.putAtt(UNITS, neutrals.temperature_unit);
        tempVar.putAtt(DESC, "Neutral temperature in K");
        output_variable_3d(startp, countp, neutrals.temperature_scgc, tempVar);

        // Output SZA
        NcVar szaVar = ncdf_file.addVar("Solar_Zenith_Angle",
          ncFloat, dimVector);
        szaVar.putAtt(UNITS, "degrees");
        szaVar.putAtt(DESC, "Solar zenith angle in degrees");
        output_variable_3d(startp, countp, grid.sza_scgc * cRtoD, szaVar);

      }

      // ----------------------------------------------
      // Ion Densities and Ion Temperature and Electron Temperature
      // ----------------------------------------------

      if (type_output == "ions" ||
          type_output == "states") {

        // Output all species densities:
        std::vector<NcVar> ionVar;
        for (int iSpecies=0; iSpecies < nIons; iSpecies++) {
          if (report.test_verbose(3))
            std::cout << "Outputting Var : "
                      << ions.species[iSpecies].cName << "\n";
          ionVar.push_back(ncdf_file.addVar("Ion_Density_" + ions.species[iSpecies].cName,
                                            ncFloat, dimVector));
          ionVar[iSpecies].putAtt(UNITS, neutrals.density_unit);
          ionVar[iSpecies].putAtt(DESC, "Specie density of " + ions.species[iSpecies].cName + " in " + neutrals.density_unit);
          output_variable_3d(startp, countp,
                             ions.species[iSpecies].density_scgc,
                             ionVar[iSpecies]);
        }

        ionVar.push_back(ncdf_file.addVar("Electron_Density_e-", ncFloat, dimVector));
        ionVar[nIons].putAtt(UNITS, neutrals.density_unit);
        ionVar[nIons].putAtt(DESC, "Specie density of e- in " + neutrals.density_unit);
        output_variable_3d(startp, countp, ions.density_scgc, ionVar[nIons]);
      }

      // ----------------------------------------------
      // Magnetic field
      // ----------------------------------------------

      if (type_output == "bfield") {
        NcVar mLatVar = ncdf_file.addVar("Magnetic_Latitude",
                                         ncFloat, dimVector);
        mLatVar.putAtt(UNITS, "radians");
        output_variable_3d(startp, countp, grid.magLat_scgc, mLatVar);

        NcVar mLonVar = ncdf_file.addVar("Magnetic_Longitude",
                                         ncFloat, dimVector);
        mLonVar.putAtt(UNITS, "radians");
        output_variable_3d(startp, countp, grid.magLat_scgc, mLonVar);

        NcVar mLTVar = ncdf_file.addVar("Magnetic_Local_Time",
                                         ncFloat, dimVector);
        mLTVar.putAtt(UNITS, "hours");
        output_variable_3d(startp, countp, grid.magLocalTime_scgc, mLTVar);

        // Output magnetic field components:

        NcVar bxVar = ncdf_file.addVar("Bx", ncFloat, dimVector);
        bxVar.putAtt(UNITS, "nT");
        output_variable_3d(startp, countp, grid.bfield_vcgc[0], bxVar);

        NcVar byVar = ncdf_file.addVar("By", ncFloat, dimVector);
        byVar.putAtt(UNITS, "nT");
        output_variable_3d(startp, countp, grid.bfield_vcgc[1], bxVar);

        NcVar bzVar = ncdf_file.addVar("Bz", ncFloat, dimVector);
        bzVar.putAtt(UNITS, "nT");
        output_variable_3d(startp, countp, grid.bfield_vcgc[2], bxVar);
      }  // if befield

      ncdf_file.close();
    }  // if time check
  }  // for iOutput

  report.exit(function);
  return iErr;
}
