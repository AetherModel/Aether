// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/times.h"
#include "../include/inputs.h"
#include "../include/neutrals.h"
#include "../include/euv.h"
#include "../include/grid.h"
#include "../include/planets.h"
#include "../include/ions.h"
#include "../include/calc_euv.h"
#include "../include/report.h"
#include "../include/output.h"


int advance( Planets &planet,
	     Grid &gGrid,
	     Times &time,
	     Euv &euv,
	     Neutrals &neutrals,
	     Ions &ions,
	     Indices &indices,
	     Inputs &input,
	     Report &report) {

  int iErr=0;

  std::string function="advance";
  report.enter(function);

  time.display();

  gGrid.calc_sza(planet, time, report);
  neutrals.calc_mass_density();
  neutrals.calc_specific_heat();
  time.calc_dt();
  
  report.enter("euv+cond");
  iErr = calc_euv(planet,
		  gGrid,
		  time,
		  euv,
		  neutrals,
		  ions,
		  indices,
		  input,
		  report);
  
  neutrals.calc_conduction(gGrid, time, report);
  report.exit("euv+cond");

  neutrals.add_sources(time, report);

  time.increment_time();

  iErr = output(neutrals, gGrid, time, planet, input, report);

  report.exit(function);
  return iErr;

}
  
