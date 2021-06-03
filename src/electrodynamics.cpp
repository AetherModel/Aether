
#include <iostream>
#include <string>
#include <netcdf>

#include "../include/aether.h"
#include <vector>
using namespace std;
//using namespace netCDF;
//using namespace netCDF::exceptions;

// g++ electrodynamics.cpp -o -std=c++11

//Initialize Electrodynamics

Electrodynamics::Electrodynamics(Inputs input, Report &report){
    read_netcdf_electrodynamics_file(input.get_electrodynamics_file(), report);
}

fcube Electrodynamics::get_potential(fcube magLat, fcube magLocalTime, Report &report){
    fcube pot(magLat.n_rows, magLat.n_cols, magLat.n_slices);
    pot.zeros();
    int time_pos = static_cast<int>(time_index);
    fmat e_potentials = input_electrodynamics[0].potential[time_pos];
    for (int i = 0; i < magLat.n_slices; ++i){
        set_grid(magLat.slice(i), magLocalTime.slice(i), report);
        pot.slice(i) = get_values(e_potentials, magLat.n_rows, magLat.n_cols);
    }
    return pot;
}

fmat Electrodynamics::get_eflux(fcube magLat, fcube magLocalTime, Report &report){
    int i = magLat.n_slices-1;//figure out
    set_grid(magLat.slice(i), magLocalTime.slice(i), report);
    int time_pos = static_cast<int>(time_index);
    fmat e_e_flux = input_electrodynamics[0].energy_flux[time_pos];
    return get_values(e_e_flux, magLat.n_rows, magLat.n_cols);
}

fmat Electrodynamics::get_avee(fcube magLat, fcube magLocalTime, Report &report){
    int i = magLat.n_slices-1;//figure out
    set_grid(magLat.slice(i), magLocalTime.slice(i), report);
    int time_pos = static_cast<int>(time_index);
    fmat e_avee = input_electrodynamics[0].average_energy[time_pos];
    return get_values(e_avee, magLat.n_rows, magLat.n_cols);
}

fmat Electrodynamics::get_values(fmat matToInterpolateOn, int rows, int cols){
    fmat slice(rows, cols);
    slice.zeros();
    for (int r = 0; r < rows; ++r){
        for (int c = 0; c < cols; ++c){
            float h_pos = input_electrodynamics[0].mlts_indices(r, c);
            float v_pos = input_electrodynamics[0].lats_indices(r, c);
            int r_start, c_start;
            if (h_pos < 0){
                c_start = 0;
            }
            else if (h_pos >= matToInterpolateOn.n_cols-1){
                c_start = matToInterpolateOn.n_cols-2;
            }
            else {
                c_start = h_pos;
            }
            if (v_pos < 0){
                r_start = 0;
            }
            else if (v_pos >= matToInterpolateOn.n_rows-1){
                r_start = matToInterpolateOn.n_rows-2;
            }
            else {
                r_start = v_pos;
            }
            float first_row_slope = matToInterpolateOn(r_start, c_start + 1) - matToInterpolateOn(r_start, c_start);
            float second_row_slope = matToInterpolateOn(r_start + 1, c_start + 1) - matToInterpolateOn(r_start + 1, c_start);
            float first_row_val = matToInterpolateOn(r_start, c_start) + first_row_slope * (h_pos - c_start);
            float second_row_val = matToInterpolateOn(r_start + 1, c_start) + second_row_slope * (h_pos - c_start);
            //time to vertically linear interpolate these two values
            float vertical_slope = second_row_val - first_row_val;
            float final_value = first_row_val + vertical_slope * (v_pos - r_start); 
            slice(r, c) = final_value;
        }
    }
    return slice;
}

//average energy and eflux energy as well call to get_values

std::tuple<fcube, fmat, fmat> Electrodynamics::get_electrodynamics(fcube magLat, fcube magLocalTime, Report &report){
    fcube pot = get_potential(magLat, magLocalTime, report);
    fmat eflux = get_eflux(magLat, magLocalTime, report);
    fmat avee = get_avee(magLat, magLocalTime, report);
    return std::make_tuple(pot, eflux, avee);
}

void Electrodynamics::set_time(double time, Report &report){
    std::string function = "Electrodynamics::set_time";
    static int iFunction = -1;
    //report.enter(function, iFunction);

    time_needed = time;

    int64_t iLow, iMid, iHigh, N;
    double interpolation_index, x, dt;
    double intime = time;
    std::vector<double> times = input_electrodynamics[0].times;
    // Check to see if the time is below the bottom time in the vector:
    iLow = 0;
    iHigh = times.size()-1;
    if (intime <= times[iLow]) {
        interpolation_index = 0.0;
        std::cout << "Warning: current time below first available potential-vector time, using first time" << "\n";
    }
    else if (intime >= times[iHigh]) {
        interpolation_index = iHigh;
        std::cout << "Warning: current time above last available potential-vector time, using last time" << "\n";
    }

    // At this point, we know that it is somewhere between the highest
    // and lowest values:
    else{
        iMid = (iHigh+iLow)/2;

        while (iHigh-iLow > 1) {
            // Break if iMid <= time < iMid+1
            if (times[iMid] == intime) break;
            if (times[iMid] <= intime &&
            times[iMid+1] > intime) break;
            // Upper Half:
            if (times[iMid] < intime) {
            iLow = iMid;
            iMid = (iHigh+iLow)/2;
            } else {
            iHigh = iMid;
            iMid = (iHigh+iLow)/2;
            }
        }

        // At this point, time should be between iMid and iMid+1:

        dt = (times[iMid+1] - times[iMid]);
        x = (intime - times[iMid]) / dt;

        interpolation_index = iMid + x;
    }
    
    time_index = interpolation_index;
    //report.exit(function);
}

void Electrodynamics::set_grid(fmat lats, fmat mlts, Report &report){
    std::string function = "Electrodynamics::set_grid";
    static int iFunction = -1;
    //report.enter(function, iFunction);

    lats_needed = lats;
    mlts_needed = mlts;

    //uses first input_electrodynamics struct
    fvec lat_search = input_electrodynamics[0].mlats;
    fvec mlt_search = input_electrodynamics[0].mlts;

    input_electrodynamics[0].lats_indices = get_interpolation_indices(lats, lat_search);
    input_electrodynamics[0].mlts_indices = get_interpolation_indices(mlts, mlt_search);
    // If we have read in a file, we need to create interpolation indices here:

    // This is a bit more complicated, since we need to loop through
    // all of the points in lats and use the 1d interpolation scheme to find
    // the index, store, and move to next point; then repeat with mlts.

    
    //report.exit(function);
}

fmat Electrodynamics::get_interpolation_indices(fmat vals, fvec search){
    fmat res(vals.n_rows, vals.n_cols, fill::zeros);

    for (int i = 0; i < vals.n_rows; ++i){
        for (int j = 0; j < vals.n_cols; ++j){

            float in = vals(i, j);
            int64_t iLow, iMid, iHigh, N;
            double interpolation_index, x, dx;

            
            // Check to see if the time is below the bottom time in the vector:
            iLow = 0;
            if (in <= search(iLow)) {
                interpolation_index = 0.0;
            }
            // Check to see if the time is above the top time in the vector:
            iHigh = search.n_rows-1;
            if (in >= search(iHigh)) {
                interpolation_index = iHigh;
            }

            // At this point, we know that it is somewhere between the highest
            // and lowest values:

            iMid = (iHigh+iLow)/2;

            while (iHigh-iLow > 1) {
                // Break if iMid <= time < iMid+1
                if (search[iMid] == in) break;
                if (search[iMid] <= in &&
                search[iMid+1] > in) break;
                // Upper Half:
                if (search[iMid] < in) {
                iLow = iMid;
                iMid = (iHigh+iLow)/2;
                } else {
                iHigh = iMid;
                iMid = (iHigh+iLow)/2;
                }
            }
            // At this point, time should be between iMid and iMid+1:

            dx = (search[iMid+1] - search[iMid]);
            x = (in - search[iMid]) / dx;

            interpolation_index = iMid + x;
            res(i, j) = interpolation_index;
        }
    }
    return res;
}


void Electrodynamics::set_imf_bx(float value){
    imf_bx_needed = value;
}

void Electrodynamics::set_imf_by(float value){
    imf_by_needed = value;
}

void Electrodynamics::set_imf_bz(float value){
    imf_bz_needed = value;
}

void Electrodynamics::set_sw_v(float value){
    sw_v_needed = value;
}

void Electrodynamics::set_sw_n(float value){
    sw_n_needed = value;
}

void Electrodynamics::set_hp(float value){
    hp_needed = value;
}

void Electrodynamics::set_au(float value){
    au_needed = value;
}

void Electrodynamics::set_al(float value){
    al_needed = value;
}

void Electrodynamics::set_ae(float value){
    ae_needed = value;
}

void Electrodynamics::set_kp(float value){
    kp_needed = value;
}


