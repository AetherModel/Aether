// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>


#include "../include/times.h"
#include "../include/inputs.h"
#include "../include/report.h"


#include "../include/neutrals.h"
#include "../include/euv.h"
#include "../include/grid.h"
#include "../include/planets.h"
#include "../include/sizes.h"
#include "../include/ions.h"
#include "../include/output.h"
#include "../include/advance.h"

int main() {

  int iErr = 0;

  Times time;
  Report report;
  Inputs input(time, report);
  Euv euv(input, report);
  Planets planet(input, report);
  Indices indices(input);

  // Geo grid stuff:
  Grid gGrid(nGeoLonsG, nGeoLatsG, nGeoAltsG);
  gGrid.init_geo_grid(planet, input, report);
  gGrid.fill_grid(planet, report);

  // Magnetic grid stuff:
  Grid mGrid(nMagLonsG, nMagLatsG, nMagAltsG);
  
  Neutrals neutrals(gGrid, input, report);
  Ions ions(input, report);
  neutrals.pair_euv(euv, ions, report);  

  // This is for the initial output.  If it is not a restart, this will go:
  if (time.check_time_gate(input.get_dt_output(0))) {
    iErr = output(neutrals, gGrid, time, planet, input, report);
  }

  // This is advancing now...

  double dt_couple = 1800.0;
  
  // The way most codes are set up in the SWMF is that there are two
  // times, an end time which ends the simulation, and an intermediate
  // time, which allows coupling or something to happen.  So, typically
  // the advance functions should only go to this intermediate time,
  // then a loop around that goes to the end time.  Then, the code can
  // be made into a library and run externally.
  
  while (time.get_current() < time.get_end()) {

    time.increment_intermediate(dt_couple);

    while (time.get_current() < time.get_intermediate())
      iErr = advance(planet,
		     gGrid,
		     time,
		     euv,
		     neutrals,
		     ions,
		     indices,
		     input,
		     report);

    // Do some coupling here. But we have no coupling to do. Sad.
    
  }

  report.times();
    
  return iErr;

}
