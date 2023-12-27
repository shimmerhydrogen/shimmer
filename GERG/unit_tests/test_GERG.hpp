#ifndef __test_GERG_H
#define __test_GERG_H

#include <iostream>
#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"

void callFevalsqrt()
{
  // Call MATLAB sqrt function on array

  using namespace matlab::engine;
  using namespace matlab::data;

  // Start MATLAB engine synchronously
  std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();

  // Create  MATLAB data array factory
  ArrayFactory factory;

  // Define a four-element array
  TypedArray<double> const argArray =
      factory.createArray({ 1,4 }, { -2.0, 2.0, 6.0, 8.0 });

  // Call MATLAB function
  TypedArray<std::complex<double>> const results =
      matlabPtr->feval(u"sqrt", argArray);

  // Display results
  int i = 0;
  for (auto r : results) {
    double a = argArray[i++];
    double realPart = r.real();
    double imgPart = r.imag();
    std::cout << "Square root of " << a << " is " <<
                 realPart << " + " << imgPart << "i" << std::endl;
  }
}

int test_GERG(int , char **)
{
  callFevalsqrt();
  return EXIT_SUCCESS;
}

#endif
