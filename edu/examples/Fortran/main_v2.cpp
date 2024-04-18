
// gfortran -c print_hi.f90 -o print_hi.o
// g++ -c main_v2.cpp -o main_v2.o
// g++ main_v2.o print_hi.o -o main -lgfortran

#include <iostream>
#include <cstring> 

extern "C" void print_hi(void);
extern "C" void print_double(int *i, float *x, float *y);
extern "C" void pass_arrays(int[], int[]);
extern "C" void get_array(int[]);
extern "C" void test_passing_string(int[]);

using namespace std;

const int iLength_ = 100;

int* copy_string_to_int(string inString) {
  const int length = inString.length(); 
  // declaring character array (+1 for null terminator) 
  int* outArray = new int[iLength_]; 
  for (int i = 0; i < length; i++) {
    outArray[i] = inString[i];
  }
  return outArray;
}

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
  
  string testString = "this is a test: UA/input/file.csv";
  int* testArray = copy_string_to_int(testString);

  test_passing_string(testArray);

  return 0;
}
