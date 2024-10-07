#ifndef __test_filter_collection_H
#define __test_filter_collection_H

#include "Eigen/Eigen"
#include "shimmer_teqp_utilities.hpp"

#include "teqp/models/GERG/GERG.hpp"
#include "test_utilities.hpp"

#include "test_gerg_mock.hpp"

namespace shimmer_teqp
{
  namespace test
  {
    // *********************************************************
    int test_filter_collection(int , char **)
    {
      using namespace shimmer_teqp::utilities;

      const auto x = gerg_mock::x();
      const auto filter = shimmer_teqp::utilities::filter_components(x,
                                                                     gerg_mock::tolerance());

      ASSERT_EQ(gerg_mock::filter(),
                filter);

      const auto& names = teqp::GERG2008::component_names;

      const auto comps = shimmer_teqp::utilities::filter_collection(filter,
                                                                    names);

      ASSERT_EQ(gerg_mock::comps(),
                comps);

      const auto mol_frac = shimmer_teqp::utilities::filter_collection(filter,
                                                                       x);

      ASSERT_VECTOR_DOUBLE_EQ(gerg_mock::mol_frac(),
                              mol_frac);

      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_H
