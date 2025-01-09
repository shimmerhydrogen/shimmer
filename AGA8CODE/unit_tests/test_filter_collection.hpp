#ifndef __test_filter_collection_H
#define __test_filter_collection_H

#include "Eigen/Eigen"
#include "shimmer_gerg_utilities.hpp"

#include "test_utilities.hpp"

#include "test_gerg_mock.hpp"

namespace shimmer_gerg
{
  namespace test
  {
    // *********************************************************
    int test_filter_collection(int , char **)
    {
      using namespace shimmer_gerg::utilities;

      const auto x = gerg_mock::x();
      const auto filter = shimmer_gerg::utilities::filter_components(x,
                                                                     gerg_mock::tolerance());

      ASSERT_EQ(gerg_mock::filter(),
                filter);

      const auto& names = shimmer_gerg::utilities::component_names;

      const auto comps = shimmer_gerg::utilities::filter_collection(filter,
                                                                    names);

      ASSERT_EQ(gerg_mock::comps(),
                comps);

      const auto mol_frac = shimmer_gerg::utilities::filter_collection(filter,
                                                                       x);

      ASSERT_VECTOR_DOUBLE_EQ(gerg_mock::mol_frac(),
                              mol_frac);

      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_H
