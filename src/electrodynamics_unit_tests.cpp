
#include <iostream>
#include <string>
#include <netcdf>
#include "../include/electrodynamics.h"
#include "../include/report.h"
#include <vector>
#include <cassert>
using namespace std;

//g++ electrodynamics_unit_tests.cpp -o e_test.exe -std=c++11

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


//integer indices
bool test1(){
    fmat vals = { {1, 3},
          {2, 4} };
    fvec search = {1,2,3,4,5};
    fmat res = get_interpolation_indices(vals, search);

    //results verification
    for (int i = 0; i < vals.n_rows; ++i){
        for (int j = 0; j < vals.n_cols; ++j){
            int ind = res(i, j);
            if (vals(i, j) != search(ind)){
                return false;
            }
        }
    }
    return true;
}

//empty
bool test2(){
    fmat vals = {};
    fvec search = {1,2,3,4,5};
    fmat res = get_interpolation_indices(vals, search);

    //results verification
    for (int i = 0; i < vals.n_rows; ++i){
        for (int j = 0; j < vals.n_cols; ++j){
            int ind = res(i, j);
            if (vals(i, j) != search(ind)){
                return false;
            }             
        }
    }
    return true;
}

//decimal indices
bool test3(){
    fmat vals = {{0, 3},{4, 5}, {7, 9}};
    fvec search = {1,3,4,8};
    fmat res = get_interpolation_indices(vals, search);

    fmat answer = {{-0.5, 1}, {2, 2.25}, {2.75, 3.25}};
    //results verification
    for (int i = 0; i < vals.n_rows; ++i){
        for (int j = 0; j < vals.n_cols; ++j){
            if (res(i, j) != answer(i, j)){
                return false;
            }
        }
    }
    return true;
}

//decimal indices
bool test4(){
    fmat vals = {{0, 3},{4, 5}};
    fvec search = {2,4,8};
    fmat res = get_interpolation_indices(vals, search);

    fmat answer = {{-1, 0.5}, {1, 1.25}};
    //results verification
    for (int i = 0; i < vals.n_rows; ++i){
        for (int j = 0; j < vals.n_cols; ++j){
            if (res(i, j) != answer(i, j)){
                return false;
            }
        }
    }
    return true;
}

int main(){
    assert(test1());
    assert(test2());
    assert(test3());
    assert(test4());
}

