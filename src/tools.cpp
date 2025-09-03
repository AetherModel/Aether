// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"


// -----------------------------------------------------------------------------
// Set all of the ghost cells to a constant value that is fed in.
//   This is primarily for testing of message passing.
// -----------------------------------------------------------------------------

void set_gcs_to_value(arma_cube &var_scgc,
                      precision_t value,
                      int64_t nGCs) {

  std::string function = "set_gcs_to_value";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iX, nX = var_scgc.n_rows;
  int64_t iY, nY = var_scgc.n_cols;
  int64_t iZ, nZ = var_scgc.n_slices;

  for (iZ = 0; iZ < nGCs; iZ++) {
    var_scgc.slice(iZ).fill(value);
    var_scgc.slice(nZ - 1 - iZ).fill(value);
  }

  // bottom:
  var_scgc.tube(0, 0, nX - 1, nGCs - 1).fill(value);
  // top:
  var_scgc.tube(0, nY - nGCs, nX - 1, nY - 1).fill(value);
  // left:
  var_scgc.tube(0, 0, nGCs - 1, nY - 1).fill(value);
  // right:
  var_scgc.tube(nX - nGCs, 0, nX - 1, nY - 1).fill(value);

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// find interpolation coefficients for a 1D interpolator
//   inX is the grid you are interpolating FROM
//   outX is the position you want to interpolate TO
//   outIndex and outRatio are the interpolation coefficents
// -----------------------------------------------------------------------------

bool find_interpolation_coefficients(arma_vec inX,
                                     arma_vec outX,
                                     arma_vec &outIndex,
                                     arma_vec &outRatio) {

  bool didWork = true;

  // Assume inX and outX are defined the same...
  int64_t iXo, iXi, nX = outX.n_rows;

  outIndex.set_size(nX);
  outRatio.set_size(nX);

  bool isFound;

  for (iXo = 0; iXo < nX; iXo++) {
    iXi = 0;
    isFound = false;

    while (!isFound && iXi < nX - 1) {
      if (inX[iXi] <= outX[iXo] &&
          inX[iXi + 1] > outX[iXo])
        isFound = true;
      else
        iXi++;
    }

    if (isFound) {
      outIndex[iXo] = iXi;
      outRatio[iXo] =
        (outX[iXo] - inX[iXi]) /
        (inX[iXi + 1] - inX[iXi]);
    } else {
      didWork = false;
      outIndex[iXo] = -1;
      outRatio[iXo] = 0.0;
    }
  }

  return didWork;
}

// -----------------------------------------------------------------------------
// This takes the index and ratio determined in the above function and
// uses them to interpolate.
// -----------------------------------------------------------------------------

arma_vec interpolate1d(arma_vec inY,
                       arma_vec &index,
                       arma_vec &ratio) {
  int64_t iY, iy_, nY = inY.n_rows;
  precision_t r_;
  arma_vec outY(nY);

  for (iY = 0; iY < nY; iY++) {
    iy_ = index(iY);
    r_ = ratio(iY);

    if (iy_ > -1)
      outY(iY) = (1.0 - r_) * inY(iy_) + r_ * inY(iy_);
  }

  return outY;
}

// ----------------------------------------------------------------------------
// Fix corners in an arma cube
//   - basically fill in the corners with values near them
// ----------------------------------------------------------------------------

void fill_horizontal_ghostcels(arma_cube &values, int64_t nGCs) {

  int64_t nXs = values.n_rows, iX;
  int64_t nYs = values.n_cols, iY;
  int64_t nZs = values.n_slices, iZ;
  int64_t iGCx, iGCy, iGCz;

  for (iGCx = 0; iGCx < nGCs; iGCx++) {
    for (iY = 0; iY < nYs; iY++) {
      // Bottom:
      values.tube(iGCx, iY) = values.tube(nGCs, iY);
      values.tube(nXs - iGCx - 1, iY) = values.tube(nXs - nGCs - 1, iY);
    }
  }

  for (iX = 0; iX < nXs; iX++) {
    for (iGCy = 0; iGCy < nGCs; iGCy++) {
      // Bottom:
      values.tube(iX, iGCy) = values.tube(iX, nGCs);
      values.tube(iX, nYs - iGCy - 1) = values.tube(iX, nYs - nGCs - 1);
    }
  }

  //fill_corners(values, nGCs);

  return;

}

// ----------------------------------------------------------------------------
// Fix corners in an arma cube
//   - basically fill in the corners with values near them
// ----------------------------------------------------------------------------

void fill_corners(arma_cube &values, int64_t nGCs) {

  int64_t nXs = values.n_rows, iX;
  int64_t nYs = values.n_cols, iY;
  int64_t nZs = values.n_slices, iZ;
  int64_t iGCx, iGCy, iGCz;

  // Bottom:
  for (iGCz == 0; iGCz < nGCs; iGCz++) {
    for (iGCx = 0; iGCx < nGCs; iGCx++) {
      for (iY = 0; iY < nYs; iY++) {
        // Bottom:
        values(iGCx, iY, iGCz) =
          values(nGCs, iY, nGCs);
        values(nXs - iGCx - 1, iY, iGCz) =
          values(nXs - nGCs - 1, iY, nGCs);
        // top:
        values(iGCx, iY, nZs - iGCz - 1) =
          values(nGCs, iY, nZs - nGCs - 1);
        values(nXs - iGCx - 1, iY, nZs - iGCz - 1) =
          values(nXs - nGCs - 1, iY, nZs - nGCs - 1);
      }
    }
  }

  for (iGCz = 0; iGCz < nGCs; iGCz++) {
    for (iGCy = 0; iGCy < nGCs; iGCy++) {
      for (iX = 0; iX < nXs; iX++) {
        // Bottoms:
        values(iX, iGCy, iGCz) =
          values(iX, nGCs, nGCs);
        values(iX, nYs - iGCy - 1, iGCz) =
          values(iX, nYs - nGCs - 1, nGCs);
        // tops:
        values(iX, iGCy, nZs - iGCz - 1) =
          values(iX, nGCs, nZs - nGCs - 1);
        values(iX, nYs - iGCy - 1, nZs - iGCz - 1) =
          values(iX, nYs - nGCs - 1, nZs - nGCs - 1);

      }
    }
  }

  for (iGCx = 0; iGCx < nGCs; iGCx++) {
    for (iGCy = 0; iGCy < nGCs; iGCy++) {
      // lower left:
      values.tube(iGCx, iGCy) = 0.5 * (
                                  values.tube(iGCx, nGCs) +
                                  values.tube(nGCs, iGCy));
      // lower right:
      values.tube(nXs - iGCx - 1, iGCy) = 0.5 * (
                                            values.tube(nXs - iGCx - 1, nGCs) +
                                            values.tube(nXs - nGCs - 1, iGCy));
      // upper left:
      values.tube(iGCx, nYs - iGCy - 1) = 0.5 * (
                                            values.tube(iGCx, nYs - nGCs - 1) +
                                            values.tube(nGCs, nYs - iGCy - 1));
      // upper right:
      values.tube(nXs - iGCx - 1, nYs - iGCy - 1) = 0.5 * (
                                                      values.tube(nXs - iGCx - 1, nYs - nGCs - 1) +
                                                      values.tube(nXs - nGCs - 1, nYs - iGCy - 1));
    }
  }

  return;
}


// ----------------------------------------------------------------------------
// Neatly display an armadillo vector
// ----------------------------------------------------------------------------

void display_vector(arma_vec vec) {
  for (int64_t i = 0; i < vec.n_rows; i++)
    std::cout << vec(i) << " ";

  std::cout << "\n";
}

// ----------------------------------------------------------------------------
// Neatly display an armadillo matrix with a name
// ----------------------------------------------------------------------------

void display_cube(std::string name, arma_cube values) {
  std::cout << name << " ";

  for (int64_t i = 0; i < values.n_slices; i++) {
    std::cout << "Slice : " << i << ":\n";
    display_matrix(" ", values.slice(i));
  }

}

// ----------------------------------------------------------------------------
// Neatly display an armadillo matrix with a name
// ----------------------------------------------------------------------------

void display_matrix(std::string name, arma_mat mat) {
  std::cout << name << "\n";

  for (int64_t i = 0; i < mat.n_cols; i++)
    display_vector(" ", mat.col(i));

  std::cout << "\n";
}



// ----------------------------------------------------------------------------
// Neatly display an armadillo vector with a name
// ----------------------------------------------------------------------------

void display_vector(std::string name, arma_vec vec) {
  std::cout << name << " ";

  for (int64_t i = 0; i < vec.n_rows; i++)
    std::cout << vec(i) << " ";

  std::cout << "\n";
}

// ----------------------------------------------------------------------------
// Neatly display a c++ vector with a name
// ----------------------------------------------------------------------------

void display_vector(std::string name, std::vector<precision_t> vec) {
  std::cout << name << " ";

  for (int64_t i = 0; i < vec.size(); i++)
    std::cout << vec[i] << " ";

  std::cout << "\n";
}


// ----------------------------------------------------------------------------
// synchronize a (boolean) variable across all processors
// ----------------------------------------------------------------------------

bool sync_across_all_procs(bool value) {
  bool global_value;
  MPI_Allreduce(&value, &global_value, 1, MPI_C_BOOL, MPI_LAND, aether_comm);
  return global_value;
}

// ----------------------------------------------------------------------------
// Find min across all processors and return value to everyone
// ----------------------------------------------------------------------------

precision_t sync_min_across_all_procs(precision_t value) {
  precision_t global_value;
  double vSend, vReceive;
  vSend = value;
  MPI_Allreduce(&vSend, &vReceive, 1, MPI_DOUBLE, MPI_MIN, aether_comm);
  global_value = vReceive;
  return global_value;
}

// ----------------------------------------------------------------------------
// Find max across all processors and return value to everyone
// ----------------------------------------------------------------------------

precision_t sync_max_across_all_procs(precision_t value) {
  precision_t global_value;
  double vSend, vReceive;
  vSend = value;
  MPI_Allreduce(&vSend, &vReceive, 1, MPI_DOUBLE, MPI_MAX, aether_comm);
  global_value = vReceive;
  return global_value;
}

// ----------------------------------------------------------------------------
// Calculate the average value across all processors
// ----------------------------------------------------------------------------

precision_t sync_mean_across_all_procs(precision_t value) {
  precision_t global_value;
  double vSend, vReceive;
  double nSend, nReceive;
  vSend = value;
  nSend = 1.0;
  MPI_Allreduce(&vSend, &vReceive, 1, MPI_DOUBLE, MPI_SUM, aether_comm);
  MPI_Allreduce(&nSend, &nReceive, 1, MPI_DOUBLE, MPI_SUM, aether_comm);
  global_value = vReceive / nReceive;
  return global_value;
}

// ----------------------------------------------------------------------------
// Calculate the average value across all processors
//   - this is the same as sync_mean_across_all_procs, but is limited to
//     processors in a given member
// ----------------------------------------------------------------------------

precision_t sync_mean_across_member(precision_t value) {
  precision_t global_value;
  double vSend, vReceive;
  double nSend, nReceive;
  vSend = value;
  nSend = 1.0;
  MPI_Allreduce(&vSend, &vReceive, 1, MPI_DOUBLE, MPI_SUM, aether_member_comm);
  MPI_Allreduce(&nSend, &nReceive, 1, MPI_DOUBLE, MPI_SUM, aether_member_comm);
  global_value = vReceive / nReceive;
  return global_value;
}
// ----------------------------------------------------------------------------
// Generate a vector of normally distributed random doubles
// ----------------------------------------------------------------------------

std::vector<double> get_normal_random_vect(double mean,
                                           double std,
                                           int64_t nValues,
                                           int seed) {
  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution(mean, std);
  std::vector<double> values(nValues);

  for (int64_t iVal = 0; iVal < nValues; iVal++)
    values[iVal] = distribution(generator);

  return values;
}

// ----------------------------------------------------------------------------
// Generate a vector of uniformly distributed random unsigned ints
// ----------------------------------------------------------------------------

std::vector<unsigned int> get_random_unsigned_vect(int64_t nValues,
                                                   int seed) {
  std::default_random_engine get_random(seed);
  std::vector<unsigned int> values(nValues);

  for (int64_t iVal = 0; iVal < nValues; iVal++)
    values[iVal] = get_random();

  return values;
}

// -----------------------------------------------------------------------------
// Compare two numbers and fail if the difference is too large
// -----------------------------------------------------------------------------

bool compare(precision_t value1, precision_t value2) {
  precision_t diff = fabs(value1 - value2);

  if (diff <= cSmall * (fabs(value1) + fabs(value2) + cSmall))
    return true;
  else
    return false;
}

// -----------------------------------------------------------------------------
// add cMember into a string just before last period
// -----------------------------------------------------------------------------

std::string add_cmember(std::string inString) {
  std::string outString = inString;
  std::size_t found = outString.rfind(".");

  if (found != std::string::npos)
    outString.replace(found, 1, "_" + cMember + ".");

  return outString;
}

// -----------------------------------------------------------------------------
// Convert an integer to a zero-padded string
// -----------------------------------------------------------------------------

std::string tostr(int64_t num_to_convert, int64_t zero_padding_len) {
  std::ostringstream ss;
  ss << std::setw( zero_padding_len ) << std::setfill( '0' ) << num_to_convert;
  return ss.str();
}

// -----------------------------------------------------------------------------
// Convert a number to a float/double
//    - Can convert scientific notation
// -----------------------------------------------------------------------------

precision_t str_to_num(std::string input) {
  std::stringstream ss(input);
  precision_t output = 0;
  ss >> output;
  return output;
}

// -----------------------------------------------------------------------
// Read json file
// -----------------------------------------------------------------------

json read_json(std::string json_file) {

  int iErr = 0;

  json json_input;
  std::ifstream infile_ptr;
  infile_ptr.open(json_file);

  if (!infile_ptr.is_open())
    std::cout << "Could not open input file: " << json_file << "!!!\n";

  else
    infile_ptr >> json_input;

  return json_input;
}

// -----------------------------------------------------------------------
// Write json file
// -----------------------------------------------------------------------

bool write_json(std::string json_file, json json_output) {

  bool DidWork = true;

  std::ofstream outfile_ptr;
  outfile_ptr.open(json_file);

  if (!outfile_ptr.is_open()) {
    std::cout << "Could not open output json file: " << json_file << "!!!\n";
    DidWork = false;
  } else
    outfile_ptr << std::setw(2) << json_output << "\n";

  return DidWork;
}

// -----------------------------------------------------------------------------
// Translate an arma_vec into a vector
// -----------------------------------------------------------------------------

std::vector<precision_t> make_vector_from_fvec(arma_vec in_fvec) {

  int64_t nPts = in_fvec.n_elem;
  std::vector<precision_t> out_vector(nPts);

  for (int64_t iPt = 0; iPt < nPts; iPt++)
    out_vector[iPt] = in_fvec(iPt);

  return out_vector;
}

// -----------------------------------------------------------------------------
// Translate a vector into an fvec
// -----------------------------------------------------------------------------

arma_vec make_fvec_from_vector(std::vector<precision_t> in_vector) {

  int64_t nPts = in_vector.size();
  arma_vec out_fvec(nPts);

  for (int64_t iPt = 0; iPt < nPts; iPt++)
    out_fvec(iPt) = in_vector[iPt];

  return out_fvec;
}

// -----------------------------------------------------------------------------
// Make a vector of arma_cubes
// -----------------------------------------------------------------------------

std::vector<arma_cube> make_cube_vector(int64_t nLons,
                                        int64_t nLats,
                                        int64_t nAlts,
                                        int64_t nComps) {
  std::vector<arma_cube> vec;
  arma_cube one_component(nLons, nLats, nAlts);
  one_component.zeros();

  for (int64_t iComp = 0; iComp < nComps; iComp++)
    vec.push_back(one_component);

  return vec;
}

// -----------------------------------------------------------------------------
// Dot product
// This assumes a 3-component vector of arma_cubes:
// -----------------------------------------------------------------------------

arma_cube dot_product(std::vector<arma_cube> vec1,
                      std::vector<arma_cube> vec2) {
  // create the dot product:
  arma_cube dot = vec1[0];
  dot.zeros();

  for (int64_t iComp = 0; iComp < 3; iComp++)
    dot = dot + vec1[iComp] % vec2[iComp];

  return dot;
}

// -----------------------------------------------------------------------------
// Cross product
// This assumes a 3-component vector of arma_cubes:
// -----------------------------------------------------------------------------

std::vector<arma_cube> cross_product(std::vector<arma_cube> vec1,
                                     std::vector<arma_cube> vec2) {
  std::vector<arma_cube> cross;
  // East:
  cross.push_back(vec1[1] % vec2[2] - vec1[2] % vec2[1]);
  // North:
  cross.push_back(-(vec1[0] % vec2[2] - vec1[2] % vec2[0]));
  // Vertical:
  cross.push_back(vec1[0] % vec2[1] - vec1[1] % vec2[0]);
  return cross;
}

// -----------------------------------------------------------------------------
// calculate mean of vector
// -----------------------------------------------------------------------------

precision_t mean(std::vector<precision_t> values) {
  int64_t nValues = values.size();
  precision_t m = 0.0;

  for (int64_t iValue = 0; iValue < nValues; iValue++)
    m = m + values[iValue];

  m = m / nValues;
  return m;
}

// --------------------------------------------------------------------
// calculate standard deviation of vector
// --------------------------------------------------------------------

precision_t standard_deviation(std::vector<precision_t> values) {
  int64_t nValues = values.size();
  precision_t m = mean(values);
  precision_t s = 0;

  for (int64_t iValue = 0; iValue < nValues; iValue++)
    s = s + (m - values[iValue]) * (m - values[iValue]);

  s = sqrt(s / nValues);
  return s;
}

//----------------------------------------------------------------------
// Get min, mean, and max of an arma_cube
//----------------------------------------------------------------------

std::vector<precision_t> get_min_mean_max(const arma_cube &value) {
  std::vector<precision_t> mmm(3);
  mmm[0] = value.min();
  mmm[1] = arma::accu(value) / value.n_elem;
  mmm[2] = value.max();
  return mmm;
}

//----------------------------------------------------------------------
// Find the name of given species in neutrals and ions.
// Throw exception if not found
//----------------------------------------------------------------------

const arma_cube& find_species_density(const std::string &name,
                                      Neutrals &neutrals,
                                      Ions &ions) {
  // Try to find the name in neutrals
  int id = neutrals.get_species_id(name);

  if (id > -1)
    return neutrals.species[id].density_scgc;

  id = ions.get_species_id(name);

  if (id > -1)
    return ions.species[id].density_scgc;

  // Throw an exception if the species is not found
  throw std::string("Can not find species named " + name);
}

//----------------------------------------------------------------------
// Get min, mean, and max of either a neutral or ion species
//----------------------------------------------------------------------

std::vector<precision_t> get_min_mean_max_density(const std::string &name,
                                                  Neutrals &neutrals,
                                                  Ions &ions) {
  return get_min_mean_max(find_species_density(name, neutrals, ions));
}

//-------------------------------------------------------------
// Checks whether two arma vectors are approximately equal
//-------------------------------------------------------------
bool is_approx_equal(arma_vec &vec1, arma_vec &vec2, precision_t tol) {
  // Check for absolute largest relative difference
  // if max diff is beyond tol, return false
  precision_t max_diff = 0.;

  // Find maximum value
  precision_t vec1_max = abs(vec1).max();
  precision_t vec2_max = abs(vec2).max();
  precision_t vec_max = std::max(vec1_max, vec2_max);

  // Check whether vectors are the same size
  // if not, return false
  if (vec1.size() != vec2.size())
    return false;

  // Loop through every member of vector
  for (int64_t i = 0; i < vec1.size(); i++) {
    precision_t curr_diff = abs(vec1(i) - vec2(i)) / vec_max;

    if (curr_diff > max_diff)
      max_diff = curr_diff;
  }

  if (max_diff > tol)
    return false;

  return true;
}

//-------------------------------------------------------------
// Overload col vector function with row vec
//-------------------------------------------------------------
bool is_approx_equal(Row<precision_t> &vec1, Row<precision_t> &vec2,
                     precision_t tol) {
  // Check for absolute largest relative difference
  // if max diff is beyond tol, return false
  precision_t max_diff = 0.;

  // Find maximum value
  precision_t vec1_max = abs(vec1).max();
  precision_t vec2_max = abs(vec2).max();
  precision_t vec_max = std::max(vec1_max, vec2_max);

  // Check whether vectors are the same size
  // if not, return false
  if (vec1.size() != vec2.size())
    return false;

  // Loop through every member of vector
  for (int64_t i = 0; i < vec1.size(); i++) {
    precision_t curr_diff = abs(vec1(i) - vec2(i)) / vec_max;

    if (curr_diff > max_diff)
      max_diff = curr_diff;
  }

  if (max_diff > tol)
    return false;

  return true;
}

//-------------------------------------------------------------
// Checks whether a vector is constant (all values the same)
// Method uses variance as evaluating factor
//-------------------------------------------------------------
bool is_approx_constant(arma_vec &vec, precision_t tol) {
  // Find variance (normalize with vector 2-norm)
  precision_t vec_norm = arma::norm(vec, 2);

  precision_t vec_var = arma::var(vec) / vec_norm;

  if (vec_var > tol)
    return false;

  return true;
}

// --------------------------------------------------------------------------
// Convert spherical vector (velocities) to reference (contravariant) vector
// Units of the velocities and transformation laws must be the same
// u and v are spherical velocities
// u1 and u2 are contravariant velocities
// --------------------------------------------------------------------------
void sphvect2ref(arma_mat& u, arma_mat& v, arma_mat& u1, arma_mat& u2,
                 mat_2x2 &A_inv_mat) {
  u1 = u % A_inv_mat.A11 + v % A_inv_mat.A12;
  u2 = u % A_inv_mat.A21 + v % A_inv_mat.A22;
}

// --------------------------------------------------------------------------
// Convert spherical vector (velocities) to reference (contravariant) vector
// Units of the velocities and transformation laws must be the same
// u and v are spherical velocities
// u1 and u2 are contravariant velocities
// --------------------------------------------------------------------------
void refvect2sph(arma_mat &u1, arma_mat &u2, arma_mat &u, arma_mat &v,
                 mat_2x2 &A_mat) {
  u = u1 % A_mat.A11 + u2 % A_mat.A12;
  v = u1 % A_mat.A21 + u2 % A_mat.A22;
}

//----------------------------------------------------------------------
// Takes a single index and finds the i, j, k position in an arma_cube
//----------------------------------------------------------------------

std::vector<int> index_to_ijk(arma_cube cube, int index) {
  uvec u = ind2sub(size(cube), index);
  int iLon = u(0);
  int iLat = u(1);
  int iAlt = u(2);
  return std::vector<int> {iLon, iLat, iAlt};
}

//----------------------------------------------------------------------
// This will find NaNs or Inf in an arma_cube and will return
// and error message if found. To be used for scalar values
//----------------------------------------------------------------------

bool all_finite(arma_cube cube, std::string name) {
  // if cube has not inf or nans, then do nothing
  if (is_finite(cube))
    return true;
  else {
    // Report where NaNs and Infs were found:
    std::vector<int> locations = indef_vector(cube);
    std::vector<int> loc = index_to_ijk(cube, locations[0]);
    std::string position =
      "(" + std::to_string(loc[0]) +
      "," + std::to_string(loc[1]) +
      "," + std::to_string(loc[2]) + ")";
    int size = locations.size();
    std::cout << "all_finite (" << name << "): " << cube(loc[0], loc[1],
                                                         loc[2]) << "\n";
    std::string error_message =
      std::to_string(size) +
      " Nonfinite values exist in " + name +
      " on iProc " + cProc +
      " and iMember " + cMember +
      " starting at: " + position;
    report.error(error_message);
    return false;
  }
}

//----------------------------------------------------------------------
// This will find NaNs or Inf in a VECTOR of arma_cubes and will return
// and error message if found. To be used for things like velocities.
//----------------------------------------------------------------------

bool all_finite(std::vector<arma_cube> cubes, std::string name) {
  bool no_nans = true;

  for (int i = 0; i < cubes.size(); ++i) {
    std::string new_name = name + "[" + std::to_string(i) + "] ";

    if (!all_finite(cubes.at(i), new_name))
      no_nans = false;
  }

  return no_nans;
}

//----------------------------------------------------------------------
// Insert a bunch of nans and inf in random places for testing
//----------------------------------------------------------------------

std::vector<int> insert_indefinites(arma_cube &cube) {
  int size = cube.n_elem;
  std::vector<int> locations;

  while (locations.size() < 6) {
    int random = rand() % size;

    if (std::find(locations.begin(),
                  locations.end(),
                  random) == locations.end())
      locations.push_back(random);
  }

  std::vector<int> nan_locations(locations.begin(), locations.begin() + 3);
  std::vector<int> indef_locations(locations.begin() + 3, locations.end());

  for (int i = 0; i < nan_locations.size(); i++) {
    cube.at(nan_locations.at(i)) = datum::nan;
    cube.at(indef_locations.at(i)) = datum::inf;
  }

  return locations;
}

//----------------------------------------------------------------------
// Loop through arma_cube and check individual cells to see if they
// are valid.  Needed this, since the included function doesn't work
// all of the time.
//----------------------------------------------------------------------

bool is_finite(arma_cube &cube) {
  for (int i = 0; i < cube.n_elem; i++) {
    if (is_nan_inf(cube.at(i)))
      return false;
  }

  return true;
}

//----------------------------------------------------------------------
// Check to see if a quantity is a NaN.  This is an explicit check
// for the bits!
//----------------------------------------------------------------------

bool is_nan(double value) {
  uint64_t bits = *reinterpret_cast<uint64_t*>(&value);
  uint64_t expMask = 0x7FF0000000000000ULL;
  uint64_t fracMask = 0x000FFFFFFFFFFFFFULL;

  return ((bits & expMask) == expMask) && ((bits & fracMask) != 0);
}

//----------------------------------------------------------------------
// Check to see if a quantity is a Inf.  This is an explicit check
// for the bits!
//----------------------------------------------------------------------

bool is_inf(double value) {
  uint64_t bits = *reinterpret_cast<uint64_t*>(&value);
  uint64_t expMask = 0x7FF0000000000000ULL;
  return (bits & expMask) == expMask;
}

//----------------------------------------------------------------------
// Check whether a value is NaN or Inf
//----------------------------------------------------------------------

bool is_nan_inf(double value) {
  return (is_nan(value) || is_inf(value));
}

//----------------------------------------------------------------------
// If the value is NaN or Inf, report its position
//----------------------------------------------------------------------

std::string print_nan_vector(std::vector<int> input, arma_cube cube) {
  std::string output("nans exist at ");
  std::vector<int> loc;

  for (int i = 0; i < 3; i++) {
    loc = index_to_ijk(cube, input.at(i));
    output += ("(" + std::to_string(loc.at(0)) +
               "," + std::to_string(loc.at(1)) +
               ","  + std::to_string(loc.at(2)) + ") ");
  }

  output += "infs exist at ";

  for (int i = 3; i < 6; i++) {
    loc = index_to_ijk(cube, input.at(i));
    output += ("(" + std::to_string(loc.at(0)) +
               "," + std::to_string(loc.at(1)) +
               ","  + std::to_string(loc.at(2)) + ") ");
  }

  output += "\n";
  return output;
}

//----------------------------------------------------------------------
//
//----------------------------------------------------------------------

std::vector<int> indef_vector(arma_cube cube) {
  std::vector<int> locations;

  for (int i = 0; i < cube.n_elem; i++) {
    if (is_nan_inf(cube.at(i)))
      locations.push_back(i);
  }

  if (locations.size() > 0)
    return locations;
  else {
    locations.push_back(-1);
    return locations;
  }
}

// --------------------------------------------------------------------------
// Project a point described by lon and lat to a point on a surface of the 2-2-2 cube
// --------------------------------------------------------------------------

arma_vec sphere_to_cube(precision_t lon_in, precision_t lat_in) {
  // See init_geo_grid.cpp:126. The offset for lon is subtracted
  lon_in = lon_in - 3 * cPI / 4;

  // Transfer polar coordinate to cartesian coordinate
  precision_t xy_temp;
  arma_vec ans(3);
  ans[2] = sin(lat_in);
  xy_temp = cos(lat_in);
  ans[1] = xy_temp * sin(lon_in);
  ans[0] = xy_temp * cos(lon_in);

  // Project this point onto the surface of cube
  precision_t coef = 1.0 / std::max({std::abs(ans[0]), std::abs(ans[1]), std::abs(ans[2])});
  ans *= coef;

  // Round the number if it is close to 1 or -1, otherwise the == and != operator
  // won't behave as expected because of the accuracy problem of floating point numbers
  for (int64_t i = 0; i < 3; ++i) {
    if (std::abs(ans[i] + 1) < cSmall)
      ans[i] = -1;

    else if (std::abs(ans[i] - 1) < cSmall)
      ans[i] = 1;
  }

  return ans;
}

////////////////////////////////////////////
// convert cell coordinates to geographic //
////////////////////////////////////////////
std::vector <arma_cube> mag_to_geo(arma_cube magLon, arma_cube magLat,
                                   arma_cube magAlt,
                                   Planets planet) {
  std::string function = "Grid::mag_to_geo";
  static int iFunction = -1;
  report.enter(function, iFunction);

  std::vector<arma_cube> llr, xyz_mag, xyz_geo, xyzRot1, xyzRot2;
  llr.push_back(magLon);
  llr.push_back(magLat);
  llr.push_back(magAlt);
  xyz_mag = transform_llr_to_xyz_3d(llr);

  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();
  std::vector<precision_t> dipole_center = planet.get_dipole_center();

  // Reverse our dipole rotations:
  xyzRot1 = rotate_around_y_3d(xyz_mag, magnetic_pole_tilt);
  xyzRot2 = rotate_around_z_3d(xyzRot1, magnetic_pole_rotation);

  // offset dipole (not fully suported yet, so will be zero)
  if ((dipole_center[0] != 0.0) || (dipole_center[1] != 0.0) ||
      (dipole_center[2] != 0.0)) {

    dipole_center = {0.0, 0.0, 0.0};
  }

  xyz_geo.push_back(xyzRot2[0] + dipole_center[0]);
  xyz_geo.push_back(xyzRot2[1] + dipole_center[1]);
  xyz_geo.push_back(xyzRot2[2] + dipole_center[2]);

  // transform back to lon, lat, radius:
  llr = transform_xyz_to_llr_3d(xyzRot2);

  report.exit(function);
  return llr;
}

////////////////////////////////////////////
//  convert cell coordinates to magnetic  //
////////////////////////////////////////////

std::vector<arma_cube> geo_to_mag(arma_cube glon,
                                  arma_cube glat,
                                  arma_cube radius,
                                  Planets &planet) {

  std::string function = "Grid::geo_to_gmag";
  static int iFunction = -1;
  report.enter(function, iFunction);

  std::vector<arma_cube> llr, xyz_mag, xyz_geo, xyzRot1, xyzRot2;
  llr.push_back(glon);
  llr.push_back(glat);
  llr.push_back(radius);
  xyz_mag = transform_llr_to_xyz_3d(llr);

  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();
  std::vector<precision_t> dipole_center = planet.get_dipole_center();

  // Reverse our dipole rotations:
  xyzRot1 = rotate_around_z_3d(xyz_mag, -magnetic_pole_tilt);
  xyzRot2 = rotate_around_y_3d(xyzRot1, -magnetic_pole_rotation);

  // offset dipole (not fully suported yet, so will be zero)
  if ((dipole_center[0] != 0.0) || (dipole_center[1] != 0.0) ||
      (dipole_center[2] != 0.0)) {

    dipole_center = {0.0, 0.0, 0.0};
  }

  xyz_geo.push_back(xyzRot2[0] - dipole_center[0]);
  xyz_geo.push_back(xyzRot2[1] - dipole_center[1]);
  xyz_geo.push_back(xyzRot2[2] - dipole_center[2]);

  // transform back to lon, lat, radius:
  llr = transform_xyz_to_llr_3d(xyzRot2);

  report.exit(function);
  return llr;
}


std::vector<precision_t> mag_to_ijk(precision_t mlon,
                                    precision_t mLat,
                                    precision_t radius,
                                    precision_t planet_radius) {

  precision_t i_lon, j_p, k_q;

  // precision_t planet_radius = planet.get_radius();

  i_lon = mlon;
  j_p = radius / planet_radius / pow(cos(mLat), 2);
  k_q = sin(mLat) / pow(radius / planet_radius, 2.);

  return {i_lon, j_p, k_q};
}

// -----------------------------------------------------------------------
// Transform a flat vector into a 1D cube (avoids having to overload everything)
// - this is overloaded for one vec/cube
// -----------------------------------------------------------------------

arma_cube vec2cube(std::vector<precision_t> ivec) {
  arma_cube outvec;
  int sizei = ivec.size();
  arma_cube I;

  I.set_size(sizei, 1, 1);

  for (int i = 0; i < sizei; i++)
    I[i] = ivec[i];

  return I;
}