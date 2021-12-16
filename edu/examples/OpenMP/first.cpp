
#include "omp.h"
#include <iostream>

// g++ -fopenmp first.cpp
// setenv OMP_NUM_THREADS 8
// a.out

int main() {
#pragma omp parallel
  {
  int64_t nProcs = omp_get_num_threads();
  int64_t ID = omp_get_thread_num();
  std::cout << "Proc " << ID << " of " << nProcs << "\n";
  }
}
