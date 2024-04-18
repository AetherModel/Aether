
#include "omp.h"
#include <iostream>

// g++ -fopenmp pi.cpp
// setenv OMP_NUM_THREADS 8
// a.out

static long nSteps = 10000000000;
double step;

int main() {

  double pi;
  step = 1.0 / (double) nSteps;

  #pragma omp parallel
  {
    int64_t nProcs = omp_get_num_threads();
    int64_t iProc = omp_get_thread_num();
    if (iProc == 0) std::cout << "Proc " << iProc << " of " << nProcs << "\n";

    int64_t i;
    double x, sum;

    sum = 0;
    for (i = iProc; i < nSteps; i += nProcs) {
      x = (i + 0.5) * step;
      sum += 4.0 / (1.0 + x * x);
    }
    #pragma omp critical
      pi += sum * step;

  }

  std::cout << "pi : " << pi << "\n";
  printf("pi : %lf\n", pi);
  
}
