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
      copy_cube_to_array(grid.geoLon_scgc, tmp_s3gc);
      lonVar.putVar(startp, countp, tmp_s3gc);
      copy_cube_to_array(grid.geoLat_scgc, tmp_s3gc);
      latVar.putVar(startp, countp, tmp_s3gc);
      copy_cube_to_array(grid.geoAlt_scgc, tmp_s3gc);
      altVar.putVar(startp, countp, tmp_s3gc);

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
	  copy_cube_to_array(neutrals.neutrals[iSpecies].density_scgc, tmp_s3gc);
	  denVar[iSpecies].putVar(startp, countp, tmp_s3gc);
	}
  
	// Output bulk temperature:
	NcVar tempVar = ncdf_file.addVar(neutrals.temperature_name, ncFloat, dimVector);
	tempVar.putAtt(UNITS,neutrals.temperature_unit);
	copy_cube_to_array(neutrals.temperature_scgc, tmp_s3gc);
	tempVar.putVar(startp, countp, tmp_s3gc);

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
	  copy_cube_to_array(ions.species[iSpecies].density_scgc, tmp_s3gc);
	  ionVar[iSpecies].putVar(startp, countp, tmp_s3gc);
	}
  
	ionVar.push_back(ncdf_file.addVar("e-", ncFloat, dimVector));
	ionVar[nIons].putAtt(UNITS,neutrals.density_unit);
	copy_cube_to_array(ions.density_scgc, tmp_s3gc);
	ionVar[nIons].putVar(startp, countp, tmp_s3gc);

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
	mLatVar.putVar(startp, countp, grid.magLat_s3gc);
	NcVar mLonVar = ncdf_file.addVar("Magnetic Longitude", ncFloat, dimVector);
	mLonVar.putAtt(UNITS,"radians");
	mLonVar.putVar(startp, countp, grid.magLon_s3gc);

	// Output magnetic field components:
	float *bfield_component_s3gc;
	long nPointsTotal = grid.get_nPointsInGrid();
	bfield_component_s3gc = (float*) malloc( nPointsTotal * sizeof(float) );
	
	NcVar bxVar = ncdf_file.addVar("Bx", ncFloat, dimVector);
	bxVar.putAtt(UNITS,"nT");
	get_vector_component(grid.bfield_v3gc, 0, IsGeoGrid, bfield_component_s3gc);
	bxVar.putVar(startp, countp, bfield_component_s3gc);
	
	NcVar byVar = ncdf_file.addVar("By", ncFloat, dimVector);
	byVar.putAtt(UNITS,"nT");
	get_vector_component(grid.bfield_v3gc, 1, IsGeoGrid, bfield_component_s3gc);
	byVar.putVar(startp, countp, bfield_component_s3gc);
	
	NcVar bzVar = ncdf_file.addVar("Bz", ncFloat, dimVector);
	bzVar.putAtt(UNITS,"nT");
	get_vector_component(grid.bfield_v3gc, 2, IsGeoGrid, bfield_component_s3gc);
	bzVar.putVar(startp, countp, bfield_component_s3gc);
	
      }
      
      ncdf_file.close();

    }
  }
  
  report.exit(function);
  return iErr;
  
}
