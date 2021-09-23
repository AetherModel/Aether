// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// -----------------------------------------------------------------------
// This is a super-simple chemistry solver that takes an implicit time-step.
// We should create more sophisticated ones, but this is ok for now.
// -----------------------------------------------------------------------

arma_cube solver_chemistry(arma_cube density, arma_cube source, arma_cube loss, precision_t dt) {
  arma_cube normalized_loss = loss / (density + 1e-6);
  arma_cube new_density = (density + dt * source) / (1.0 + dt * normalized_loss);
  return new_density;
}
