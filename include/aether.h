// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_AETHER_H_
#define INCLUDE_AETHER_H_

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

#include "electrodynamics.h"

#include "aurora.h"
#include "euv.h"
#include "calc_euv.h"
#include "chemistry.h"
#include "read_collision_file.h"
#include "output.h"
#include "advance.h"

#include "solvers.h"
#include "tools.h"
#include "transform.h"

#include "calc_grid_derived.h"

#endif  // INCLUDE_AETHER_H_
