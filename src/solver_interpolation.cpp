// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <vector>

#include "../include/aether.h"

// ----------------------------------------------------------------------
// This is the general function for getting an index within an ordered
// vector:
// - Returns the index and the linear interpolation coefficient as a float
// - If time is below the first time, returns 0.0
// - If time is beyond the last value, returns (n-1).0
// - Otherwise returns N.x, where:
//   value = (1.0 - x) * array[N] + x * array[N+1]
// ----------------------------------------------------------------------

double interpolate_1d_get_index_doubles(double intime,
					std::vector<double> times) {

  int64_t iLow, iMid, iHigh, N;
  double interpolation_index, x, dt;

  // Check to see if the time is below the bottom time in the vector:
  iLow = 0;
  if (intime <= times[iLow]) {
    interpolation_index = 0.0;
    return interpolation_index;
  }

  // Check to see if the time is above the top time in the vector:
  iHigh = sizeof(times)-1;
  if (intime >= times[iHigh]) {
    interpolation_index = iHigh;
    return interpolation_index;
  }

  // At this point, we know that it is somewhere between the highest
  // and lowest values:

  iMid = (iHigh+iLow)/2;

  while (iHigh-iLow > 1) {
    // Break if iMid <= time < iMid+1
    if (times[iMid] == intime) break;
    if (times[iMid] <= intime &&
	times[iMid+1] > intime) break;
    // Upper Half:
    if (times[iMid] < intime) {
      iLow = iMid;
      iMid = (iHigh+iLow)/2;
    } else {
      iHigh = iMid;
      iMid = (iHigh+iLow)/2;
    }
  }

  // At this point, time should be between iMid and iMid+1:

  dt = (times[iMid+1] - times[iMid]);
  x = (intime - times[iMid]) / dt;

  interpolation_index = iMid + x;
  return interpolation_index;
}

// ----------------------------------------------------------------------
//
// ----------------------------------------------------------------------

double interpolate_1d_w_index(std::vector<double> values,
			      double interpolation_index,
			      int interpolation_type) {
  
  int64_t n = interpolation_index;
  double x = (interpolation_index - n);
  
  if (interpolation_type == iPrevious_) {
    // Ignore x completely:
    return values[n];
  }
  if (interpolation_type == iNext_) {
    // If x is 0, stick with the current number, else go to the next one:
    if (x > 0) n++;
    return values[n];
  }
  if (interpolation_type == iClosest_) {
    // Rounding is done by adding 0.5 and taking the int:
    n = n + (x+0.5);
    return values[n];
  }
  // Interpolate properly:
  return (1.0-x) * values[n] + x * values[n+1];
}

// ----------------------------------------------------------------------
//
// ----------------------------------------------------------------------

double interpolate_1d_w_index(std::vector<float> values,
			      double interpolation_index,
			      int interpolation_type) {
  
  int64_t n = interpolation_index;
  double x = (interpolation_index - n);
  
  if (interpolation_type == iPrevious_) {
    // Ignore x completely:
    return values[n];
  }
  if (interpolation_type == iNext_) {
    // If x is 0, stick with the current number, else go to the next one:
    if (x > 0) n++;
    return values[n];
  }
  if (interpolation_type == iClosest_) {
    // Rounding is done by adding 0.5 and taking the int:
    n = n + (x+0.5);
    return values[n];
  }
  // Interpolate properly:
  return (1.0-x) * values[n] + x * values[n+1];
}

// ----------------------------------------------------------------------
//
// ----------------------------------------------------------------------

double interpolate_1d_w_index(std::vector<float> values,
			      float interpolation_index,
			      int interpolation_type) {
  
  int64_t n = interpolation_index;
  double x = (interpolation_index - n);
  
  if (interpolation_type == iPrevious_) {
    // Ignore x completely:
    return values[n];
  }
  if (interpolation_type == iNext_) {
    // If x is 0, stick with the current number, else go to the next one:
    if (x > 0) n++;
    return values[n];
  }
  if (interpolation_type == iClosest_) {
    // Rounding is done by adding 0.5 and taking the int:
    n = n + (x+0.5);
    return values[n];
  }
  // Interpolate properly:
  return (1.0-x) * values[n] + x * values[n+1];
}

// ----------------------------------------------------------------------
//
// ----------------------------------------------------------------------

double interpolate_1d_w_index(fvec values,
			      double interpolation_index,
			      int interpolation_type) {
  
  int64_t n = interpolation_index;
  double x = (interpolation_index - n);
  
  if (interpolation_type == iPrevious_) {
    // Ignore x completely:
    return values(n);
  }
  if (interpolation_type == iNext_) {
    // If x is 0, stick with the current number, else go to the next one:
    if (x > 0) n++;
    return values(n);
  }
  if (interpolation_type == iClosest_) {
    // Rounding is done by adding 0.5 and taking the int:
    n = n + (x+0.5);
    return values(n);
  }
  // Interpolate properly:
  return (1.0-x) * values(n) + x * values(n+1);
}

// ----------------------------------------------------------------------
// This is for a series of 2d slices (matrix) as an array.
// ----------------------------------------------------------------------

fmat interpolate_3d_w_index(std::vector<fmat> values,
			    double interpolation_index,
			    int interpolation_type) {
  
  int64_t n = interpolation_index;
  double x = (interpolation_index - n);
  
  if (interpolation_type == iPrevious_) {
    // Ignore x completely:
    return values[n];
  }
  if (interpolation_type == iNext_) {
    // If x is 0, stick with the current number, else go to the next one:
    if (x > 0) n++;
    return values[n];
  }
  if (interpolation_type == iClosest_) {
    // Rounding is done by adding 0.5 and taking the int:
    n = n + (x+0.5);
    return values[n];
  }
  // Interpolate properly:
  return (1.0-x) * values[n] + x * values[n+1];
}
