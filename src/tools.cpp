// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"


// ----------------------------------------------------------------------------
// Neatly display an armadillo vector
// ----------------------------------------------------------------------------

void display_vector(arma_vec vec) {
  for (int64_t i = 0; i < vec.n_rows; i++)
    std::cout << vec(i) << " ";
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
  global_value = vReceive/nReceive;
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
  if (id > -1) {
    return neutrals.species[id].density_scgc;
  }

  // Try to find the name in ions
  id = ions.get_species_id(name);
  if (id > -1) {
    return ions.species[id].density_scgc;
  }

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

//----------------------------------------------------------------------
// Takes a single index and finds the i, j, k position in an arma_cube
//----------------------------------------------------------------------

std::vector<int> index_to_ijk(arma_cube cube, int index){
  arma::uword x = cube.n_rows;
  arma::uword y = cube.n_cols;
  int altitude = index/(x*y);
  int remainder = index%(x*y);
  int lattitude = remainder/y;
  int longitude = remainder%y;
  return std::vector<int> {lattitude, longitude, altitude};
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
  for(int i=0; i<cubes.size(); ++i){
    std::string new_name = name + "[" + std::to_string(i) + "] ";
    if (!all_finite(cubes.at(i), new_name)){
      no_nans = false;
    }
  }
  return no_nans;
}

//----------------------------------------------------------------------
// Insert a bunch of nans and inf in random places for testing
//----------------------------------------------------------------------

std::vector<int> insert_indefinites(arma_cube &cube) {
  int size = cube.n_elem;
  std::vector<int> locations;
  while (locations.size() < 6){
    int random = rand() % size;
    if (std::find(locations.begin(),
		  locations.end(),
		  random) == locations.end())
      locations.push_back(random);
  }
  std::vector<int> nan_locations(locations.begin(), locations.begin() + 3);
  std::vector<int> indef_locations(locations.begin() + 3, locations.end());
  for (int i = 0; i < nan_locations.size(); i++){
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
  for (int i=0; i<cube.n_elem; i++){
    if (is_nan_inf(cube.at(i))){
      return false;
    }
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

std::string print_nan_vector(std::vector<int> input, arma_cube cube){
  std::string output("nans exist at ");
  std::vector<int> loc;
  for (int i = 0; i < 3; i++){
    loc = index_to_ijk(cube, input.at(i));
    output += ("(" + std::to_string(loc.at(0)) +
	       "," + std::to_string(loc.at(1)) +
	       ","  + std::to_string(loc.at(2)) + ") ");
  }
  output += "infs exist at ";
  for (int i = 3; i < 6; i++){
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

std::vector<int> indef_vector(arma_cube cube){
  std::vector<int> locations;
  for (int i=0; i<cube.n_elem; i++){
    if (is_nan_inf(cube.at(i))){
      locations.push_back(i);
    }
  }
  if (locations.size() > 0){
    return locations;
  }
  else{
    locations.push_back(-1);
    return locations;
  }
}
