// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_READ_INDICES_H_
#define INCLUDE_READ_INDICES_H_

#include <vector>
#include <string>

#include "indices.h"
#include "report.h"

index_file_output_struct read_f107_file(std::string f107_file,
					Indices indices,
					Report &report);
index_file_output_struct read_omni_file(std::string omni_file,
					Indices indices,
					Report &report);
int pair_omniweb_index(std::string, Indices index);
float get_omniweb_missing_value(std::string);

int read_and_store_indices(Indices &indices, Inputs args, Report &report);


#endif  // INCLUDE_READ_INDICES_H_
