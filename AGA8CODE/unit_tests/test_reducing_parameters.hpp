#ifndef __test_reducing_parameters_H
#define __test_reducing_parameters_H

#include "Eigen/Eigen"
#include "shimmer_gerg_functions.hpp"
#include "shimmer_gerg_utilities.hpp"

#include "test_gerg_mock.hpp"
#include "test_utilities.hpp"

namespace shimmer_gerg
{
  namespace test
  {
    // *********************************************************
    int test_reducing_parameters(int , char **)
    {
      using namespace shimmer_gerg::utilities;

      const auto tolerance = gerg_mock::tolerance();
      const auto comps = gerg_mock::comps();
      const auto mol_frac = gerg_mock::mol_frac();
      const auto x = gerg_mock::x();

      shimmer_gerg::gerg_functions::setup_GERG();

      const auto reducing_parameters = shimmer_gerg::gerg_functions::reducing_parameters(x,
                                                                                         tolerance);
      const auto expected_reducing_parameters = gerg_mock::reducing_parameters();

      ASSERT_DOUBLE_EQ_TOL(expected_reducing_parameters.Tr,
                           reducing_parameters.Tr,
                           tolerance);
      ASSERT_DOUBLE_EQ_TOL(expected_reducing_parameters.Dr,
                           reducing_parameters.Dr,
                           tolerance);

      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_H
