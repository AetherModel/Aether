// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_TOOLS_H_
#define INCLUDE_TOOLS_H_

// ----------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------

void display_vector(arma_vec vec);

// -----------------------------------------------------------------------------
// synchronize a variable across all processors
// -----------------------------------------------------------------------------

bool sync_across_all_procs(bool value);

// -----------------------------------------------------------------------------
// Generate a vector of normally distributed random doubles
// -----------------------------------------------------------------------------

std::vector<double> get_normal_random_vect(double mean,
					   double std,
					   int64_t nValues,
					   int seed);

// -----------------------------------------------------------------------------
// Generate a vector of uniformly distributed random unsigned ints
// -----------------------------------------------------------------------------

std::vector<unsigned int> get_random_unsigned_vect(int64_t nValues,
						   int seed);

std::vector<arma_cube> make_cube_vector(int64_t nLons,
				    int64_t nLats,
				    int64_t nAlts,
				    int64_t nComps);

arma_cube dot_product(std::vector<arma_cube> vec1,
		      std::vector<arma_cube> vec2);

std::vector<arma_cube> cross_product(std::vector<arma_cube> vec1,
				     std::vector<arma_cube> vec2);

std::vector<precision_t> make_vector_from_fvec(arma_vec in_fvec);
arma_vec make_fvec_from_vector(std::vector<precision_t> in_vector);

std::string tostr(int64_t num_to_convert, int64_t zero_padding_len);

precision_t str_to_num(std::string input);

json read_json(std::string json_file);
bool write_json(std::string json_file, json json_output);

bool compare(precision_t value1, precision_t value2);

// -----------------------------------------------------------------------------
// calculate mean of vector
// -----------------------------------------------------------------------------

precision_t mean(std::vector<precision_t> values);

// -----------------------------------------------------------------------------
// calculate standard deviation of vector
// -----------------------------------------------------------------------------

precision_t standard_deviation(std::vector<precision_t> values);

//-------------------------------------------------------------
// Get min, mean, and max of an arma_cube
//-------------------------------------------------------------

std::vector<precision_t> get_min_mean_max(const arma_cube &value);

//-------------------------------------------------------------
// Find the name of given species in neutrals and ions. Throw exception if not found
//-------------------------------------------------------------

const arma_cube& find_species_density(const std::string &name,
                                      Neutrals &neutrals,
                                      Ions &ions);

//-------------------------------------------------------------
// Get min, mean, and max of either a neutral or ion species
//-------------------------------------------------------------

std::vector<precision_t> get_min_mean_max_density(const std::string &name,
                                                  Neutrals &neutrals,
                                                  Ions &ions);

#endif  // INCLUDE_TOOLS_H_
