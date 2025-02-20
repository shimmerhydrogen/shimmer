#ifndef __test_filter_components_H
#define __test_filter_components_H

#include "Eigen/Eigen"
#include "shimmer_gerg_utilities.hpp"

#include "test_utilities.hpp"

namespace shimmer_gerg
{
  namespace test
  {
    // *********************************************************
    int test_filter_components(int , char **)
    {
      using namespace shimmer_gerg::utilities;

      constexpr double tolerance = std::numeric_limits<double>::epsilon();
      Eigen::ArrayXd x(21);
      x<< 0.8, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

      const auto filter = shimmer_gerg::utilities::filter_components(x,
                                                                     tolerance);


      const std::vector<unsigned int> expected_filter = { 0, 2 };
      ASSERT_EQ(expected_filter,
                filter);

      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_H
