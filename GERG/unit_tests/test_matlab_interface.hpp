#ifndef __test_matlab_interface_H
#define __test_matlab_interface_H

#include <iostream>
#include "Matlab_interface.hpp"
#include "test_utilities.hpp"

namespace GERG_test
{
  int test_matlab_interface(int , char **)
  {
    const std::string matlab_shimmer_folder = GERG::Matlab_interface::matlab_shimmer_directory_path();
    const auto& matlab = GERG::Matlab_interface::get_instance(
                           {
                             matlab_shimmer_folder
                           });

    ASSERT_TRUE(matlab.is_directory_in_matlab_path(matlab_shimmer_folder));

    return EXIT_SUCCESS;
  }
}

#endif // __test_matlab_interface_H
