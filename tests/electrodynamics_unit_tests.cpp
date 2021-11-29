// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <string>
#include <netcdf>
#include "../include/aether.h"
#include "../include/electrodynamics.h"
#include "../include/report.h"
#include <vector>
#include <cassert>
using namespace std;

// -----------------------------------------------------------------------------
// To run these unit tests perform:
// g++ electrodynamics_unit_tests.cpp -o e_test.exe -std=c++11
//
//
//  Binomial search to find interpolation indices.
// -----------------------------------------------------------------------------

fmat get_interpolation_indices(fmat vals, fvec search) {
  fmat res(vals.n_rows, vals.n_cols, fill::zeros);

  for (int i = 0; i < vals.n_rows; ++i) {
    for (int j = 0; j < vals.n_cols; ++j) {

      float in = vals(i, j);
      int64_t iLow, iMid, iHigh, N;
      double interpolation_index, x, dx;


      // Check to see if the time is below the bottom time in the vector:
      iLow = 0;

      if (in <= search(iLow))
        interpolation_index = 0.0;

      // Check to see if the time is above the top time in the vector:
      iHigh = search.n_rows - 1;

      if (in >= search(iHigh))
        interpolation_index = iHigh;

      // At this point, we know that it is somewhere between the highest
      // and lowest values:

      iMid = (iHigh + iLow) / 2;

      while (iHigh - iLow > 1) {
        // Break if iMid <= time < iMid+1
        if (search[iMid] == in)
          break;

        if (search[iMid] <= in &&
            search[iMid + 1] > in)
          break;

        // Upper Half:
        if (search[iMid] < in) {
          iLow = iMid;
          iMid = (iHigh + iLow) / 2;
        } else {
          iHigh = iMid;
          iMid = (iHigh + iLow) / 2;
        }
      }

      // At this point, time should be between iMid and iMid+1:

      dx = (search[iMid + 1] - search[iMid]);
      x = (in - search[iMid]) / dx;

      interpolation_index = iMid + x;
      res(i, j) = interpolation_index;
    }
  }

  return res;
}

// -----------------------------------------------------------------------------
// Test 2d interpolation
// -----------------------------------------------------------------------------

float two_d_interpolate(fmat vals, float row, float col) {
  fmat e_potentials = vals;
  float h_pos = col;
  float v_pos = row;
  int r_start, c_start;

  if (h_pos < 0)
    c_start = 0;
  else if (h_pos >= vals.n_cols - 1)
    c_start = vals.n_cols - 2;
  else
    c_start = h_pos;

  if (v_pos < 0)
    r_start = 0;
  else if (v_pos >= vals.n_rows - 1)
    r_start = vals.n_rows - 2;
  else
    r_start = v_pos;

  float first_row_slope = e_potentials(r_start,
                                       c_start + 1) - e_potentials(r_start, c_start);
  float second_row_slope = e_potentials(r_start + 1,
                                        c_start + 1) - e_potentials(r_start + 1, c_start);
  float first_row_val = e_potentials(r_start,
                                     c_start) + first_row_slope * (h_pos - c_start);
  float second_row_val = e_potentials(r_start + 1,
                                      c_start) + second_row_slope * (h_pos - c_start);
  //time to vertically linear interpolate these two values
  float vertical_slope = second_row_val - first_row_val;
  float final_value = first_row_val + vertical_slope * (v_pos - r_start);
  return final_value;
}

//----------------------------------------------------------------------------
// Tests start below

//integer indices
bool test_interpolation_indices_1() {
  fmat vals = { {1, 3},
    {2, 4}
  };
  fvec search = {1, 2, 3, 4, 5};
  fmat res = get_interpolation_indices(vals, search);

  //results verification
  for (int i = 0; i < vals.n_rows; ++i) {
    for (int j = 0; j < vals.n_cols; ++j) {
      int ind = res(i, j);

      if (vals(i, j) != search(ind))
        return false;
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
// Test for interpolation indices
// -----------------------------------------------------------------------------

bool test_interpolation_indices_2() {
  fmat vals = {};
  fvec search = {1, 2, 3, 4, 5};
  fmat res = get_interpolation_indices(vals, search);

  //results verification
  for (int i = 0; i < vals.n_rows; ++i) {
    for (int j = 0; j < vals.n_cols; ++j) {
      int ind = res(i, j);

      if (vals(i, j) != search(ind))
        return false;
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
// Test for interpolation indices
// -----------------------------------------------------------------------------

bool test_interpolation_indices_3() {
  fmat vals = {{1.5, 3}, {4, 5}, {7, 9}};
  fvec search = {1, 3, 4, 8};
  fmat res = get_interpolation_indices(vals, search);

  fmat answer = {{0.25, 1}, {2, 2.25}, {2.75, 3.25}};

  //results verification
  for (int i = 0; i < vals.n_rows; ++i) {
    for (int j = 0; j < vals.n_cols; ++j) {
      if (res(i, j) != answer(i, j))
        return false;
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
// Test for interpolation indices
// -----------------------------------------------------------------------------

bool test_interpolation_indices_4() {
  fmat vals = {{0, 3}, {4, 5}};
  fvec search = {2, 4, 8};
  fmat res = get_interpolation_indices(vals, search);

  fmat answer = {{-1, 0.5}, {1, 1.25}};

  //results verification
  for (int i = 0; i < vals.n_rows; ++i) {
    for (int j = 0; j < vals.n_cols; ++j) {
      if (res(i, j) != answer(i, j))
        return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
// Testing a bunch of different functions below
// -----------------------------------------------------------------------------

bool test_two_d_interpolate_1() {
  fmat vals = { {0, 1, 2},
    {0, 1, 2}
  };
  float res = two_d_interpolate(vals, 0.5, 0.5);

  return res == 0.5;
}

bool test_two_d_interpolate_2() {
  fmat vals = { {0, 1, 2},
    {0, 2, 3}
  };
  float res = two_d_interpolate(vals, 0.5, 0.5);
  return res == 0.75;
}

bool test_two_d_interpolate_3() {
  fmat vals = { {2, 4, 5},
    {6, 7, 8}
  };
  float res = two_d_interpolate(vals, 0, 1.5);
  return res == 4.5;
}

bool test_two_d_interpolate_4() {
  fmat vals = { {2, 4, 5},
    {6, 7, 8}
  };
  float res = two_d_interpolate(vals, -1, 0);
  return res == -2;
}

bool test_two_d_interpolate_5() {
  fmat vals = { {2, 4, 5},
    {6, 7, 8}
  };
  float res = two_d_interpolate(vals, 0.5, -1);
  return res == 2.5;
}

bool test_two_d_interpolate_6() {
  fmat vals = { {2, 4, 5},
    {6, 7, 8}
  };
  float res = two_d_interpolate(vals, 2, 2);
  return res == 11;
}

bool test_two_d_interpolate_7() {
  fmat vals = { {2, 4, 5},
    {6, 7, 8}
  };
  float res = two_d_interpolate(vals, 2, 3);
  return res == 12;
}

bool test_two_d_interpolate_8() {
  fmat vals = { {2, 4, 5},
    {6, 7, 8}
  };
  float res = two_d_interpolate(vals, -1, 4);
  return res == 4;
}

int main() {
  assert(test_interpolation_indices_1());
  assert(test_interpolation_indices_2());
  assert(test_interpolation_indices_3());
  assert(test_interpolation_indices_4());

  assert(test_two_d_interpolate_1());
  assert(test_two_d_interpolate_2());
  assert(test_two_d_interpolate_3());
  assert(test_two_d_interpolate_4());
  assert(test_two_d_interpolate_5());
  assert(test_two_d_interpolate_6());
  assert(test_two_d_interpolate_7());
  assert(test_two_d_interpolate_8());
}

