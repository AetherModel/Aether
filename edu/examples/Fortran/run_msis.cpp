
// These codes are used to call MSIS from C++
// Developed by A. Ridley at University of Michigan, April, 2023

// To compile:
// gfortran -c msis_constants.F90
// gfortran -c msis_gfn.F90
// gfortran -c msis_tfn.F90
// gfortran -c msis_dfn.F90
// gfortran -c msis_calc.F90
// gfortran -c msis_gtd8d.F90
// gfortran -c msis_init.F90
// gfortran -c msis_utils.F90
// gfortran -c msis2.1_test.F90
// 
// This is for the fortran version of main:
// gfortran -o msis2.1 msis*.o
// 
// This is for the C++ version of main:
// g++ -c run_msis.cpp
// gfortran -c call_msis.f90
// 
// g++ -o msis_c -lgfortran run_msis.o msis_dfn.o msis_init.o msis_calc.o msis_gfn.o msis_tfn.o msis_constants.o msis_gtd8d.o msis_utils.o call_msis.o

#include <iostream>

// Define the two fortran codes we have to call:

extern "C" void init_msis(void);
extern "C" void call_msis_f(int *iYear,
			    int *iDay,
			    float *second,
			    float *gLonDeg,
			    float *gLatDeg,
			    float *altKm,
			    float *f107,
			    float *f107a,
			    float *ap,
			    float[], float[]);

using namespace std;

int main() {

  // Call init msis in fortran:
  init_msis();
  
  // This is the first test case in the msis2.1_test_in.txt file:
  
  // 70178  64800    0.2   50.0   55.0  21.67  153.3  146.5   35.0
  int iYear = 70;
  int iDay = 178;
  float second = 64800.0;
  float altKm = 0.2;
  float gLatDeg = 50.0;
  float gLonDeg = 55.0;
  float f107 = 146.5;
  float f107a = 153.3;
  float ap = 35.0;

  float density_back[10];
  float temperature_back[2];

  // Call msis in fortran:
  call_msis_f(&iYear, &iDay, &second, &gLonDeg, &gLatDeg, &altKm,
	      &f107, &f107a, &ap, density_back, temperature_back);

  // Report the densities and temperatures returned:
  for (int i = 0; i < 10; i++)
    std::cout << "c++ (density) " << i << ": " << density_back[i] << "\n";
  std::cout << "c++ (temp) : "
	    << temperature_back[0] << " "
	    << temperature_back[1] << "\n";

  // Answers should be:
  /* 
     inputs:
     70178  64800    0.2   50.0   55.0  21.67  153.3  146.5   35.0
     outputs:
     0.1255E+15 0.9999E-37 0.1884E+20 0.5052E+19 0.2252E+18
     0.1160E-02 0.9999E-37 0.9999E-37 0.9999E-37 0.9999E-37  
     294.10
  */
  
  return 0;
  
}
