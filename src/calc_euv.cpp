// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <vector>

#include "../include/times.h"
#include "../include/inputs.h"
#include "../include/indices.h"
#include "../include/euv.h"
#include "../include/report.h"
#include "../include/neutrals.h"
#include "../include/ions.h"

int calc_euv( Planets planet,
	      Grid grid,
	      Times time,
	      Euv &euv,
	      Neutrals &neutrals,
	      Ions &ions,
	      Indices indices,
	      Inputs args,
	      Report &report) {
  
  int iErr=0;
  
  if (time.check_time_gate(args.get_dt_euv())) {

    std::string function="Euv::calc_euv";
    report.enter(function);
    
    // Chapman integrals for EUV energy deposition:
    neutrals.calc_chapman(grid, report);
  
    iErr = euv.euvac(time, indices, report);
    iErr = euv.scale_from_1au(planet, time);
  
    neutrals.calc_ionization_heating(euv, ions, report);

    report.exit(function);
  }

  return iErr;

}

