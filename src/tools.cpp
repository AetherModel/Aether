// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"


// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------

void display_vector(arma_vec vec) {

  for (int64_t i = 0; i < vec.n_rows; i++)
    std::cout << vec(i) << " ";

  std::cout << "\n";

}


// ----------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------

bool sync_across_all_procs(bool value) {
  bool global_value;
  MPI_Allreduce(&value, &global_value, 1, MPI_C_BOOL, MPI_LAND, aether_comm);
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

// -----------------------------------------------------------------------------
// calculate standard deviation of vector
// -----------------------------------------------------------------------------

precision_t standard_deviation(std::vector<precision_t> values) {
  int64_t nValues = values.size();
  precision_t m = mean(values);
  precision_t s = 0;

  for (int64_t iValue = 0; iValue < nValues; iValue++)
    s = s + (m - values[iValue]) * (m - values[iValue]);

  s = sqrt(s / nValues);
  return s;
}

//-------------------------------------------------------------
// Get min, mean, and max of an arma_cube
//-------------------------------------------------------------

std::vector<precision_t> get_min_mean_max(const arma_cube &value) {
  std::vector<precision_t> mmm(3);
  mmm[0] = value.min();
  mmm[1] = arma::accu(value) / value.n_elem;
  mmm[2] = value.max();
  return mmm;
}

//-------------------------------------------------------------
// Find the name of given species in neutrals and ions. Throw exception if not found
//-------------------------------------------------------------

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

//-------------------------------------------------------------
// Get min, mean, and max of either a neutral or ion species
//-------------------------------------------------------------

std::vector<precision_t> get_min_mean_max_density(const std::string &name,
                                                  Neutrals &neutrals,
                                                  Ions &ions) {
  return get_min_mean_max(find_species_density(name, neutrals, ions));
}
