#ifndef __test_filter_collection_H
#define __test_filter_collection_H

#include "Eigen/Eigen"
#include "shimmer_teqp_utilities.hpp"

#include "teqp/models/GERG/GERG.hpp"
#include "test_utilities.hpp"

namespace shimmer_teqp
{
  namespace test
  {
    // *********************************************************
    int test_filter_collection(int , char **)
    {
      using namespace shimmer_teqp::utilities;

      const std::vector<unsigned int> filter = { 0, 2 };

      Eigen::ArrayXd x(21);
      x<< 0.8, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      const auto& names = teqp::GERG2008::component_names;

      const auto comps = shimmer_teqp::utilities::filter_collection(filter,
                                                                    names);

      const std::vector<std::string> expected_comps = { "methane", "carbondioxide" };
      ASSERT_EQ(expected_comps, comps);

      const auto mol_frac = shimmer_teqp::utilities::filter_collection(filter,
                                                                       x);


      Eigen::ArrayXd expected_mol_frac(2);
      expected_mol_frac<< 0.8, 0.2;
      ASSERT_VECTOR_DOUBLE_EQ(expected_mol_frac, expected_mol_frac);

      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_H
