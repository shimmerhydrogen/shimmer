#ifndef __test_reducing_parameters_H
#define __test_reducing_parameters_H

#include <iostream>

#include "x_mock.hpp"
#include "MATLAB_GERG_functions.hpp"

namespace GERG_test
{
  int test_reducing_parameters(int , char **)
  {
    const auto x = x_mock();

    const auto reducing_parameters = GERG::reducing_parameters(x);

    std::cout<< "CALL test_reducing_parameters"<< std::endl;


    return EXIT_SUCCESS;
  }
}

#endif // __test_reducing_parameters_H
