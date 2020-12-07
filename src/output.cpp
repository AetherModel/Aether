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

using namespace netCDF;
using namespace netCDF::exceptions;

int output(Neutrals neutrals,
	   Grid grid,
	   Times time,
	   Planets planet,
	   Inputs args,
	   Report &report) {

  int iErr = 0;
    
  if (time.check_time_gate(args.get_dt_output(0))) {
 
    std::string function="output";
    report.enter(function);

    std::string time_string;
    std::string file_name;
    std::string file_pre = "3DALL";
    std::string file_ext = ".nc";
    std::string UNITS = "units";

    time_string = time.get_YMD_HMS();
    file_name = file_pre + "_" + time_string + file_ext;
  
    // Create the file:
    NcFile ncdf_file(file_name, NcFile::replace);

    // Add dimensions:
    NcDim lonDim = ncdf_file.addDim("Longitude", nGeoLonsG); 
    NcDim latDim = ncdf_file.addDim("Latitude", nGeoLatsG); 
    NcDim altDim = ncdf_file.addDim("Altitude", nGeoAltsG); 

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

    countp.push_back(nGeoLonsG);
    countp.push_back(nGeoLatsG);
    countp.push_back(nGeoAltsG);

    // Output longitude, latitude, altitude 3D arrays:
    lonVar.putVar(startp, countp, grid.geoLon_s3gc);
    latVar.putVar(startp, countp, grid.geoLat_s3gc);
    altVar.putVar(startp, countp, grid.geoAlt_s3gc);

    // Output all species densities:
    std::vector<NcVar> denVar;
    for (int iSpecies=0; iSpecies < nSpecies; iSpecies++) {
      if (report.test_verbose(3))
	std::cout << "Outputting Var : "
		  << neutrals.neutrals[iSpecies].cName << "\n";
      denVar.push_back(ncdf_file.addVar(neutrals.neutrals[iSpecies].cName, ncFloat, dimVector));
      denVar[iSpecies].putAtt(UNITS,neutrals.density_unit);
      denVar[iSpecies].putVar(startp, countp, neutrals.neutrals[iSpecies].density_s3gc);
  }
  
    // Output bulk temperature:
    NcVar tempVar = ncdf_file.addVar(neutrals.temperature_name, ncFloat, dimVector);
    tempVar.putAtt(UNITS,neutrals.temperature_unit);
    tempVar.putVar(startp, countp, neutrals.temperature_s3gc);

    ncdf_file.close();

    report.exit(function);

  }
  
  return iErr;
  
}
