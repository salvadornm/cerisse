#include "clad/Differentiator/Differentiator.h"
#include <iostream>


inline double f(double x, double y) { return x * y; }

int main() {
  // Call clad to generate the derivative of f wrt x.
  auto f_dx = clad::differentiate(f, "x");
  // Execute the generated derivative function.
  std::cout << f_dx.execute(/*x=*/3, /*y=*/4) << std::endl;
  // Dump the generated derivative code to standard output.
  f_dx.dump();

  exit(0);

  }