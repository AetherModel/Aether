// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_TOOLS_H_
#define INCLUDE_TOOLS_H_

// Structure for a 2x2 matrix
struct mat_2x2{
    arma_mat A11; 
    arma_mat A12; 
    arma_mat A21; 
    arma_mat A22;
};

// ----------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------

void display_vector(arma_vec vec);

// -----------------------------------------------------------------------------
// synchronize a variable across all processors
// -----------------------------------------------------------------------------

bool sync_across_all_procs(bool value);

// ----------------------------------------------------------------------------
// Calculate the average value across all processors
// ----------------------------------------------------------------------------

precision_t sync_mean_across_all_procs(precision_t value);

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

//-------------------------------------------------------------
// Checks whether two arma vectors are approximately equal
//-------------------------------------------------------------
bool is_approx_equal(arma_vec &vec1, arma_vec &vec2, precision_t tol);

//-------------------------------------------------------------
// Overload col vector function with row vec
//-------------------------------------------------------------
bool is_approx_equal(Row<precision_t> &vec1, Row<precision_t> &vec2, precision_t tol);

//-------------------------------------------------------------
// Checks whether a vector is constant (all values the same)
// Method uses variance as a factor
//-------------------------------------------------------------
bool is_approx_constant(arma_vec &vec, precision_t tol);

// --------------------------------------------------------------------------
// Convert spherical vector (velocities) to reference (contravariant) vector
// Units of the velocities and transformation laws must be the same
// u and v are spherical velocities
// u1 and u2 are contravariant velocities
// --------------------------------------------------------------------------
void sphvect2ref(arma_mat& u, arma_mat& v, arma_mat& u1, arma_mat& u2, mat_2x2 &A_inv_mat);

// --------------------------------------------------------------------------
// Convert spherical vector (velocities) to reference (contravariant) vector
// Units of the velocities and transformation laws must be the same
// u and v are spherical velocities
// u1 and u2 are contravariant velocities
// --------------------------------------------------------------------------
void refvect2sph(arma_mat &u1, arma_mat &u2, arma_mat &u, arma_mat &v, mat_2x2 &A_mat);

#endif  // INCLUDE_TOOLS_H_
