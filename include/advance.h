// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_ADVANCE_H_
#define INCLUDE_ADVANCE_H_

/**************************************************************
 *
 * advance:
 *
 * - Function that advances the states in Aether by one time step
 *
 *   Pretty much all classes have to be passed into this function
 *   this function calls a bunch of functions that alters all
 *   of the states in the system.
 **************************************************************/

#include "../include/aether.h"

/**********************************************************************
  * brief advance the states in the model
  *
  * param planet all of the planet infomation about the simulation
  * param gGrid the geographic grid class to use
  * param time the time class to use
  * param euv the euv class to use
  * param neutrals the class that contains all info about the neutrals
  * param ions the class that contains all info about the ions
  * param chemistry the class that contains all of the info about chemistry
  * param electrodynamics the class that contains all of the electrodynamics
  * param indices the class that contains all of the indices
 **/


int advance(Planets &planet,
            Grid &gGrid,
            Times &time,
            Euv &euv,
            Neutrals &neutrals,
            Ions &ions,
            Chemistry &chemistry,
            Electrodynamics &electrodynamics,
            Indices &indices,
            Logfile &logfile);

#endif // INCLUDE_ADVANCE_H_
