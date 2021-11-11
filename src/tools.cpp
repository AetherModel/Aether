// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Convert an integer to a zero-padded string
// -----------------------------------------------------------------------------

std::string tostr(int64_t num_to_convert, int64_t zero_padding_len) {
  std::ostringstream ss;
  ss << std::setw( zero_padding_len ) << std::setfill( '0' ) << num_to_convert;
  return ss.str();
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
