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

#ifndef AETHER_INCLUDE_H_
#define AETHER_INCLUDE_H_

// std libraries
#include <vector>

// Aether libraries
include "constants.h"
include "file_input.h"
include "inputs.h"
include "sizes.h"
include "time_conversion.h"
include "times.h"

#endif // AETHER_INCLUDE_H_
