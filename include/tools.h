// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_TOOLS_H_
#define INCLUDE_TOOLS_H_

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

json read_json(std::string json_file);
bool write_json(std::string json_file, json json_output);

#endif  // INCLUDE_TOOLS_H_
