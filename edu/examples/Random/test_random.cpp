// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <chrono>
#include <random>

/*  

This is to demonstrate and test random numbers.

g++ test_random.cpp
a.out

Produces:

seed : 1971
Roll # 0
  seed : 33126597
  Input mean and test mean : 10 and 9.97498
  Input std and test std : 2 and 1.99518
Roll # 1
  seed : 560451206
  Input mean and test mean : 10 and 9.98721
  Input std and test std : 2 and 2.0023
Roll # 2
  seed : 640143500
  Input mean and test mean : 10 and 9.9607
  Input std and test std : 2 and 1.99034
etc.

 */


std::vector<double> get_normal_random_vect(double mean,
					   double std,
					   int64_t nValues,
					   int seed) {
  std::default_random_engine generator (seed); 
  std::normal_distribution<double> distribution(mean, std);
  std::vector<double> values(nValues);
  for (int64_t iVal = 0; iVal < nValues; iVal++)
    values[iVal] = distribution(generator);
  return values;
}

std::vector<unsigned int> get_random_unsigned_vect(int64_t nValues,
						   int seed) {
  std::default_random_engine get_random(seed);
  std::vector<unsigned int> values(nValues);
  for (int64_t iVal = 0; iVal < nValues; iVal++)
    values[iVal] = get_random();
  return values;
}


int main() {

  int seed = int(std::chrono::system_clock::now().time_since_epoch().count());
  seed = 1971;
  std::cout << "seed : " << seed << "\n";

  int nTimes = 10;
  
  std::vector<unsigned int> list_of_seeds = get_random_unsigned_vect(nTimes,
								     seed);
  
  double random_mean = 10.0;
  double random_std = 2.0;
  int nVals, iVal;

  for (int iTime = 0; iTime < nTimes; iTime++) {

    seed = list_of_seeds[iTime];
    
    std::vector<double> values = get_normal_random_vect(random_mean,
							random_std,
							nVals,
							seed);
    double mean_test = 0.0;
    double std_test = 0.0;

    for (iVal = 0; iVal < nVals; iVal++)
      mean_test = mean_test + values[iVal];
    mean_test = mean_test / nVals;

    for (iVal = 0; iVal < nVals; iVal++)
      std_test = std_test +
	(mean_test - values[iVal]) * (mean_test - values[iVal]);
    std_test = sqrt(std_test / nVals);

    std::cout << "Roll # " << iTime << "\n";
    std::cout << "  seed : " << seed << "\n";

    std::cout << "  Input mean and test mean : "
	      << random_mean << " and " << mean_test << "\n";
  
    std::cout << "  Input std and test std : "
	      << random_std << " and " << std_test << "\n";
  
  }
  
  return 0;

}
