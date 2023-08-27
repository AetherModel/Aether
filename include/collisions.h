// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_COLLISIONS_H_
#define INCLUDE_COLLISIONS_H_

#include <vector>
#include <string>

#include "neutrals.h"
#include "ions.h"

void calc_ion_neutral_coll_freq(Neutrals &neutrals,
                 Ions &ions);

#endif  // INCLUDE_COLLISIONS_H_
