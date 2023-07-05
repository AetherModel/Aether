// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// Initial version: A. Ridley, June 2023

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Need to find the minimum of three numbers
// --------------------------------------------------------------------------

precision_t min_of_three(precision_t one,
			 precision_t two,
			 precision_t three) {
  if (one <= two && one <= three)
    return one;
  // one is not the smallest!
  if (two <= three)
    return two;
  else
    return three;
}

// --------------------------------------------------------------------------
// Need to find the maximum of three numbers
// --------------------------------------------------------------------------

precision_t max_of_three(precision_t one,
			 precision_t two,
			 precision_t three) {
  if (one >= two && one >= three)
    return one;
  // one is not the biggest!
  if (two >= three)
    return two;
  else
    return three;
}

// --------------------------------------------------------------------------
// MC Limiter, beta as input
// --------------------------------------------------------------------------

precision_t limiter_mc(precision_t dUp,
		       precision_t dDown,
		       precision_t beta) {
  precision_t limited = 0.0;
  precision_t ave;
  if (dUp > 0) {
    if (dDown > 0) {
      ave = 0.5 * (dUp + dDown);
      limited = min_of_three(beta * dUp, beta * dDown, ave);
    } // else limited = 0.0
  } else {
    if (dDown < 0) {
      ave = 0.5 * (dUp + dDown);
      limited = max_of_three(beta * dUp, beta * dDown, ave);
    } // else limited = 0.0
  }
  return limited;
}

