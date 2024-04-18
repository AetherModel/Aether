// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <vector>

#include "../include/aether.h"

double interpolate_1d(double outX,
                      std::vector<double> inXs,
                      std::vector<double> inValues) {
  double ind = interpolate_1d_get_index_doubles(outX, inXs);
  double outValue = interpolate_1d_w_index(inValues, ind, iInterp_);
  return outValue;
}

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
  iHigh = times.size() - 1;

  if (intime >= times[iHigh]) {
    interpolation_index = iHigh;
    return interpolation_index;
  }

  // At this point, we know that it is somewhere between the highest
  // and lowest values:

  iMid = (iHigh + iLow) / 2;

  while (iHigh - iLow > 1) {
    // Break if iMid <= time < iMid+1
    if (times[iMid] == intime)
      break;

    if (times[iMid] <= intime &&
        times[iMid + 1] > intime)
      break;

    // Upper Half:
    if (times[iMid] < intime) {
      iLow = iMid;
      iMid = (iHigh + iLow) / 2;
    } else {
      iHigh = iMid;
      iMid = (iHigh + iLow) / 2;
    }
  }

  // At this point, time should be between iMid and iMid+1:

  dt = (times[iMid + 1] - times[iMid]);
  x = (intime - times[iMid]) / dt;

  interpolation_index = iMid + x;
  return interpolation_index;
}

// ----------------------------------------------------------------------
// These functions are overloaded, in that they are called the same
// thing but take different inputs - mostly a mixture of float/double
// arrays, vectors, and armadillo types.
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
// Conducts interpolation using one of 4 methods
// (vector is double, index is double)
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
    if (x > 0)
      n++;

    return values[n];
  }

  if (interpolation_type == iClosest_) {
    // Rounding is done by adding 0.5 and taking the int:
    n = n + (x + 0.5);
    return values[n];
  }

  // Interpolate properly:
  return (1.0 - x) * values[n] + x * values[n + 1];
}

// ----------------------------------------------------------------------
// Conducts interpolation using one of 4 methods
// (vector is float, index is double)
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
    if (x > 0)
      n++;

    return values[n];
  }

  if (interpolation_type == iClosest_) {
    // Rounding is done by adding 0.5 and taking the int:
    n = n + (x + 0.5);
    return values[n];
  }

  // Interpolate properly:
  return (1.0 - x) * values[n] + x * values[n + 1];
}

// ----------------------------------------------------------------------
// Conducts interpolation using one of 4 methods
// (vector is float, index is float)
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
    if (x > 0)
      n++;

    return values[n];
  }

  if (interpolation_type == iClosest_) {
    // Rounding is done by adding 0.5 and taking the int:
    n = n + (x + 0.5);
    return values[n];
  }

  // Interpolate properly:
  return (1.0 - x) * values[n] + x * values[n + 1];
}

// ----------------------------------------------------------------------
// Conducts interpolation using one of 4 methods
// (vector is fvec, index is double)
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
    if (x > 0)
      n++;

    return values(n);
  }

  if (interpolation_type == iClosest_) {
    // Rounding is done by adding 0.5 and taking the int:
    n = n + (x + 0.5);
    return values(n);
  }

  // Interpolate properly:
  return (1.0 - x) * values(n) + x * values(n + 1);
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
    if (x > 0)
      n++;

    return values[n];
  }

  if (interpolation_type == iClosest_) {
    // Rounding is done by adding 0.5 and taking the int:
    n = n + (x + 0.5);
    return values[n];
  }

  // Interpolate properly:
  return (1.0 - x) * values[n] + x * values[n + 1];
}

// ----------------------------------------------------------------------
// Fundamental calculation formula of linear interpolation
// If we have (x0, y0) and (x1, y1), the value at x within [x0, x1) is estimated to be
// (x1 - x) / (x1 - x0) * y0 + (x - x0) / (x1 - x0) * y1
// = (1 - ratio) * y0 + ratio * y1
// where ratio = (x - x0) / (x1 - x0)
// ----------------------------------------------------------------------

precision_t linear_interpolation(const precision_t y0,
                                 const precision_t y1,
                                 const precision_t ratio) {
  return (1.0 - ratio) * y0 + ratio * y1;
}

// ----------------------------------------------------------------------
// Estimate the value of a point inside the cube using linear interpolation
// Requirements: xRatio, yRatio, zRatio within [0, 1]
// In the documentation of arma, 3d cube is ordered slice by slice, then column by column
// i.e data[1] is the value of (0, 1, 0), data[3] is the value of (1, 1, 0)
// ----------------------------------------------------------------------

precision_t interpolate_unit_cube(const arma_cube &data,
                                  const precision_t xRatio,
                                  const precision_t yRatio,
                                  const precision_t zRatio) {

  // check the number of elements
  if (data.n_rows != 2 || data.n_cols != 2 || data.n_slices != 2)
    return std::numeric_limits<precision_t>::quiet_NaN();

  // interpolate along the x axis, calculate the value at
  // (xRatio, 0, 0), (xRatio, 1, 0), (xRatio, 0, 1), (xRatio, 1, 1)
  precision_t yzPlane[4];

  for (int64_t i = 0; i < 4; ++i)
    yzPlane[i] = linear_interpolation(data[2 * i], data[2 * i + 1], xRatio);

  // interpolate along the y axis, calculate the value at
  // (xRatio, yRatio, 0), (xRatio, yRatio, 1)
  precision_t zLine[2];

  for (int64_t i = 0; i < 2; ++i)
    zLine[i] = linear_interpolation(yzPlane[2 * i], yzPlane[2 * i + 1], yRatio);

  // return the value at (xRatio, yRatio, zRatio)
  return linear_interpolation(zLine[0], zLine[1], zRatio);
}
