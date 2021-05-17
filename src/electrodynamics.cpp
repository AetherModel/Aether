
#include <iostream>
#include <string>
#include <netcdf>
#include "read_netcdf_electrodynamics.cpp"
#include "../include/electrodynamics.h"
#include "../include/report.h"
#include <vector>
using namespace std;
//using namespace netCDF;
//using namespace netCDF::exceptions;

// g++ electrodynamics.cpp -o -std=c++11


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
    if (intime <= times[iLow]) {
        interpolation_index = 0.0;
    }

    // Check to see if the time is above the top time in the vector:
    iHigh = sizeof(times)-1;
    if (intime >= times[iHigh]) {
        interpolation_index = iHigh;
    }

    // At this point, we know that it is somewhere between the highest
    // and lowest values:

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

    // If we have read in a file, we need to create interpolation indices here:

    // Please look at the function Indices::get_index to find the
    // general technique for doing linear interpolation.  We want to
    // break this in two, where the first step is finding the index and linear
    // coefficient to use (iMid and x, where the function can return iMid+x)
    // and the end part where is actually uses iMid+x to do the interpolation
    
    //report.exit(function);
}

fmat get_interpolation_indices(fmat vals, fvec search){
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


