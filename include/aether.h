// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_AETHER_H_
#define INCLUDE_AETHER_H_

/// The armadillo library is to allow the use of 3d cubes and other
/// array types, with array math built in. This eliminates loops!
#include <armadillo>

/// This is used for timing and the random seed generator:
#include <chrono>

/// This is for generating random numbers:
#include <random>

// Types
// Precision compile-time aliasing
#ifdef AETHER_USE_PRECISION_DOUBLE
/// Precision type chosen to be `double` through `AETHER_USE_PRECISION_DOUBLE`
using precision_t = double;
#else
/// Precision type compile-time default to float.
using precision_t = float;
#endif

/// Armadillo type vector (single column) with compile-time precision.
using arma_vec = arma::Col<precision_t>;
/// Armadillo type matrix (two dimension) with compile-time precision.
using arma_mat = arma::Mat<precision_t>;
/// Armadillo type cube (three dimension) with compile-time precision.
using arma_cube = arma::Cube<precision_t>;


#include <nlohmann/json.hpp>
using json = nlohmann::json;

// This is for manipulating IO
#include <iomanip>

// This is for manipulating strings
#include <sstream>

// Aether includes
#include "earth.h"

// Contains all information about time in the code and wall time:
#include "times.h"

// Contains a reporting system for the model include verbose level and timing:
#include "report.h"

// not done
#include "inputs.h"

// Defines physical and conversion constants
#include "constants.h"

// Needs to be eliminated (but not right now!)
#include "sizes.h"

// These are functions that convert time between different systems.
#include "time_conversion.h"

// Functions that assist in the reading of files.
#include "file_input.h"

// A class for keeping track of indices (1d vectors w/time)
#include "indices.h"

// Read indices (f107, omni) file types
#include "read_indices_files.h"

// A class for keeping track of all of the planetary characteristics
#include "planets.h"

// not done
#include "grid.h"

// Contains the neutral states and derived quantities
#include "neutrals.h"

// not done
#include "ions.h"

// not done
#include "bfield.h"

// Contains the electrodynamic states (potential and aurora)
#include "electrodynamics.h"

// not done
#include "aurora.h"

// Defines the Extreme Ultraviolet radiation above the atmosphere
#include "euv.h"

// not done
#include "calc_euv.h"
// not done
#include "chemistry.h"
// not done
#include "read_collision_file.h"
// not done
#include "output.h"
// not done
#include "advance.h"

// not done
#include "solvers.h"
// not done
#include "tools.h"
// not done
#include "transform.h"

// not done
#include "calc_grid_derived.h"

// not done
#include "parallel.h"

#endif  // INCLUDE_AETHER_H_
