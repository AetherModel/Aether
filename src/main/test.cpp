// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/aether.h"
// -----------------------------------------------------------------------------
// Main file for testing the Aether model.
// This needs more work!
// -----------------------------------------------------------------------------

int main() {

  int iErr = 0;

  // ------------------------------------------------------------
  // Test time routines:
  // ------------------------------------------------------------

  iErr = test_time_routines();

  if (iErr == 0)
    std::cout << "Passed test_time_routines!\n";
  else
    std::cout << "Failed test_time_routines!\n";

  return iErr;
}
