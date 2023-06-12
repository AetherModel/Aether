
// gfortran -c print_hi.f90 -o print_hi.o
// g++ -c main_v2.cpp -o main_v2.o
// g++ main_v2.o print_hi.o -o main -lgfortran

#include <iostream>

extern "C" void print_hi(void);
extern "C" void print_double(int *i, float *x, float *y);
extern "C" void pass_arrays(int[], int[]);
extern "C" void get_array(int[]);

using namespace std;

int main() {
  int j;
  float x, y;
  j = 4;
  x = 3.14;
  y = 0.0;
  print_hi();
  print_double(&j, &x, &y);
  cout << "float returned from Fortran: " << y << endl;

  int array_to[10];
  int array_back[10];
  for (int i = 0; i < 10; i++)
    array_to[i] = i+1;

  pass_arrays(array_to, array_back);
  for (int i = 0; i < 10; i++) {
    std::cout << "c++ " << i << ": " << array_back[i] << "\n";
    array_back[i] = 0;
  }
  
  get_array(array_back);
  for (int i = 0; i < 10; i++)
    std::cout << "c++ (saved fortran) " << i << ": " << array_back[i] << "\n";
  
  return 0;
}
