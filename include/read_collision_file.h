// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_COLLISION_H_
#define INCLUDE_COLLISION_H_

#include "../include/aether.h"

void read_collision_file(Neutrals &neutrals,
			 Ions &ions);

void parse_nu_in_table(std::vector<std::vector<std::string>> csv,
		       Neutrals &neutrals,
		       Ions &ions);

void parse_resonant_nu_in_table(std::vector<std::vector<std::string>> csv,
				Neutrals &neutrals,
				Ions &ions);

void parse_bst_in_table(std::vector<std::vector<std::string>> csv,
				Neutrals &neutrals,
				Ions &ions);

void check_collision_frequncies(Ions ions,
				Neutrals neutrals);

#endif  // INCLUDE_IONS_H_

