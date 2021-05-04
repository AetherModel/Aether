// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// -----------------------------------------------------------------------
// This is a super-simple chemistry solver that takes an implicit time-step.
// We should create more sophisticated ones, but this is ok for now.
// -----------------------------------------------------------------------

fcube solver_chemistry(fcube density, fcube source, fcube loss, float dt) {
  fcube normalized_loss = loss / (density + 1e-6);
  fcube new_density = (density + dt * source) / (1.0 + dt * normalized_loss);
  return new_density;
}
