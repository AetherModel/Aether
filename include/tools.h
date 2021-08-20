// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_TOOLS_H_
#define INCLUDE_TOOLS_H_

std::vector<fcube> make_cube_vector(int64_t nLons,
				    int64_t nLats,
				    int64_t nAlts,
				    int64_t nComps);

fcube dot_product(std::vector<fcube> vec1,
		  std::vector<fcube> vec2);

std::vector<fcube> cross_product(std::vector<fcube> vec1,
				 std::vector<fcube> vec2);

std::vector<float> make_vector_from_fvec(fvec in_fvec);
fvec make_fvec_from_vector(std::vector<float> in_vector);

#endif  // INCLUDE_TOOLS_H_
