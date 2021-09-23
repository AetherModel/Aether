// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_AURORA_H_
#define INCLUDE_AURORA_H_

void read_aurora(Neutrals &neutrals,  
		 Ions &ions,
		 Inputs args,
		 Report &report);

fvec calculate_fang(float eflux,  // in ergs/cm2/s
                    float avee,   // in keV
                    float Ebin,   // eV
                    fvec rhoH, 
                    std::vector<float> Ci, 
                    float dE,     // eV
                    fvec H, 
                    Report &report);

void calc_aurora(Grid grid, 
		 Neutrals &neutrals, 
		 Ions &ions, 
		 Inputs args,
		 Report &report);

#endif  // INCLUDE_AURORA_H_

