// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <netcdf>

#include "../include/neutrals.h"
#include "../include/grid.h"
#include "../include/times.h"
#include "../include/planets.h"
#include "../include/inputs.h"
#include "../include/earth.h"
#include "../include/report.h"
#include "../include/transform.h"

using namespace netCDF;
using namespace netCDF::exceptions;

void output_variable_3d(std::vector<size_t> count_start,
			std::vector<size_t> count_end,
			fcube value,
			NcVar variable) {

  long nX = value.n_rows;
  long nY = value.n_cols;
  long nZ = value.n_slices;
  long iX, iY, iZ, iTotal, index;

  iTotal = nX * nY * nZ; 
  
  float *tmp_s3gc = (float*) malloc( iTotal * sizeof(float) );

  for (iX = 0; iX < nX; iX++) {
    for (iY = 0; iY < nY; iY++) {
      for (iZ = 0; iZ < nZ; iZ++) {
	index = iX*nY*nZ + iY*nZ + iZ;
	tmp_s3gc[index] = value(iX,iY,iZ);
      }
    }
  }

  variable.putVar(count_start, count_end, tmp_s3gc);
    
}
			

int output(Neutrals neutrals,
	   Ions ions,
	   Grid grid,
	   Times time,
	   Planets planet,
	   Inputs args,
	   Report &report) {

  int iErr = 0;

  int nOutputs = args.get_n_outputs();
  int IsGeoGrid = grid.get_IsGeoGrid();

  long nLons = grid.get_nLons();
  long nLats = grid.get_nLats();
  long nAlts = grid.get_nAlts();
  long iTotal = long(nLons) * long(nLats) * long(nAlts);

  float *tmp_s3gc = (float*) malloc( iTotal * sizeof(float) );
  
  std::string function="output";
  static int iFunction = -1;
  report.enter(function, iFunction);  

  for (int iOutput = 0; iOutput < nOutputs; iOutput++) {

    if (time.check_time_gate(args.get_dt_output(iOutput))) {
 
      std::string time_string;
      std::string file_name;
      std::string file_ext = ".nc";
      std::string UNITS = "units";
      std::string file_pre;

      std::string type_output = args.get_type_output(iOutput);

      if (type_output == "neutrals") file_pre = "3DNEU";
      if (type_output == "states") file_pre = "3DALL";
      if (type_output == "bfield") file_pre = "3DBFI";

      time_string = time.get_YMD_HMS();
      file_name = file_pre + "_" + time_string + file_ext;
  
      // Create the file:
      NcFile ncdf_file(file_name, NcFile::replace);

      // Add dimensions:
      NcDim lonDim = ncdf_file.addDim("Longitude", nLons); 
      NcDim latDim = ncdf_file.addDim("Latitude", nLats); 
      NcDim altDim = ncdf_file.addDim("Altitude", nAlts); 

      // If we wanted 1D variables, we would do something like this, but
      // since all of out variables will be 3d, skip this:
      // Define the Coordinate Variables
      //NcVar altVar = ncdf_file.addVar("Altitude", ncFloat, altDim);
  
      // Define the netCDF variables for the 3D data.
      // First create a vector of dimensions:
      std::vector<NcDim> dimVector;
      dimVector.push_back(lonDim);
      dimVector.push_back(latDim);
      dimVector.push_back(altDim);

      NcVar lonVar = ncdf_file.addVar("Longitude", ncFloat, dimVector);
      NcVar latVar = ncdf_file.addVar("Latitude", ncFloat, dimVector);
      NcVar altVar = ncdf_file.addVar("Altitude", ncFloat, dimVector);

      lonVar.putAtt(UNITS,"radians");
      latVar.putAtt(UNITS,"radians");
      altVar.putAtt(UNITS,"meters");

      std::vector<size_t> startp,countp;
      startp.push_back(0);
      startp.push_back(0);
      startp.push_back(0);

      countp.push_back(nLons);
      countp.push_back(nLats);
      countp.push_back(nAlts);

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
		      << neutrals.neutrals[iSpecies].cName << "\n";
	  denVar.push_back(ncdf_file.addVar(neutrals.neutrals[iSpecies].cName, ncFloat, dimVector));
	  denVar[iSpecies].putAtt(UNITS,neutrals.density_unit);
	  output_variable_3d(startp, countp, neutrals.neutrals[iSpecies].density_scgc, denVar[iSpecies]);
	}
  
	// Output bulk temperature:
	NcVar tempVar = ncdf_file.addVar(neutrals.temperature_name, ncFloat, dimVector);
	tempVar.putAtt(UNITS,neutrals.temperature_unit);
	output_variable_3d(startp, countp, neutrals.temperature_scgc, tempVar);

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
	  ionVar.push_back(ncdf_file.addVar(ions.species[iSpecies].cName, ncFloat, dimVector));
	  ionVar[iSpecies].putAtt(UNITS,neutrals.density_unit);
	  output_variable_3d(startp, countp, ions.species[iSpecies].density_scgc, ionVar[iSpecies]);
	}
  
	ionVar.push_back(ncdf_file.addVar("e-", ncFloat, dimVector));
	ionVar[nIons].putAtt(UNITS,neutrals.density_unit);
	output_variable_3d(startp, countp, ions.density_scgc, ionVar[nIons]);

	// // Output bulk temperature:
	// NcVar tempVar = ncdf_file.addVar(neutrals.temperature_name, ncFloat, dimVector);
	// tempVar.putAtt(UNITS,neutrals.temperature_unit);
	// tempVar.putVar(startp, countp, neutrals.temperature_s3gc);

      }

      // ----------------------------------------------
      // Magnetic field
      // ----------------------------------------------

      if (type_output == "bfield") {
	NcVar mLatVar = ncdf_file.addVar("Magnetic Latitude", ncFloat, dimVector);
	mLatVar.putAtt(UNITS,"radians");
	output_variable_3d(startp, countp, grid.magLat_scgc, mLatVar);

	NcVar mLonVar = ncdf_file.addVar("Magnetic Longitude", ncFloat, dimVector);
	mLonVar.putAtt(UNITS,"radians");
	output_variable_3d(startp, countp, grid.magLat_scgc, mLonVar);

	// Output magnetic field components:
	
	NcVar bxVar = ncdf_file.addVar("Bx", ncFloat, dimVector);
	bxVar.putAtt(UNITS,"nT");
	output_variable_3d(startp, countp, grid.bfield_vcgc[0], bxVar);
	
	NcVar byVar = ncdf_file.addVar("By", ncFloat, dimVector);
	byVar.putAtt(UNITS,"nT");
	output_variable_3d(startp, countp, grid.bfield_vcgc[1], bxVar);
	
	NcVar bzVar = ncdf_file.addVar("Bz", ncFloat, dimVector);
	bzVar.putAtt(UNITS,"nT");
	output_variable_3d(startp, countp, grid.bfield_vcgc[2], bxVar);
	
      }
      
      ncdf_file.close();

    }
  }
  
  report.exit(function);
  return iErr;
  
}
