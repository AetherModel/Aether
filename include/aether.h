#ifndef INCLUDE_AETHER_H_
#define INCLUDE_AETHER_H_

/*! \mainpage Aether: Thermosphere-Ionosphere model
 *
 * \section intro_sec Introduction
 *
 * Aether is a Thermosphere and Ionosphere model where the two regions are
 * coupled through different grids.
 *
 * \section install_sec Installation
 *
 * Run make.
 *
 */

/*! \file aether.h
    \brief Main include file that includes other parts.
    
    Top level include file.
*/
// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// The armadillo library is to allow the use of 3d cubes and other
// array types, with array math built in. This eliminates loops!

#include <armadillo>
using namespace arma;

#include "times.h"
#include "report.h"
#include "inputs.h"
#include "constants.h"
#include "sizes.h"
#include "time_conversion.h"
#include "file_input.h"
#include "indices.h"
#include "read_f107_file.h"
#include "planets.h"
#include "grid.h"

#include "neutrals.h"
#include "ions.h"

#include "bfield.h"

#include "euv.h"
#include "calc_euv.h"
#include "chemistry.h"
#include "output.h"
#include "advance.h"

#include "solvers.h"
#include "transform.h"

#endif  // INCLUDE_AETHER_H_
