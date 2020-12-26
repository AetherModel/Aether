// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <string>

#include "../include/neutrals.h"
#include "../include/times.h"
#include "../include/inputs.h"
#include "../include/report.h"

void Neutrals::add_sources( Times time, Report &report) {

  std::string function="add_sources";
  static int iFunction = -1;
  report.enter(function, iFunction);  

  long iLon, iLat, iAlt, index;

  float dt = time.get_dt();
  
  for (iLon = 0; iLon < nGeoLonsG; iLon++) {
    for (iLat = 0; iLat < nGeoLatsG; iLat++) {
      for (iAlt = 0; iAlt < nGeoAltsG; iAlt++) {
	
	index = ijk_geo_s3gc(iLon,iLat,iAlt);

	temperature_s3gc[index] =
	  temperature_s3gc[index] +
	  dt * ( heating_euv_s3gc[index] +
		 conduction_s3gc[index]);

      }
    }
  }
  
  report.exit(function);  
  return;
  
}
