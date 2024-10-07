#ifndef __test_reducing_parameters_H
#define __test_reducing_parameters_H

#include "Eigen/Eigen"
#include "shimmer_teqp_utilities.hpp"

#include "teqp/models/GERG/GERG.hpp"
#include "test_utilities.hpp"

#include <iostream>

namespace shimmer_teqp
{
  namespace test
  {
    // *********************************************************
    int test_reducing_parameters(int , char **)
    {
      using namespace shimmer_teqp::utilities;

      constexpr double tolerance = 1.0e-14;
      Eigen::ArrayXd x(21);
      x<< 0.8, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

      const auto filter = shimmer_teqp::utilities::filter_components(x,
                                                                     tolerance);

      const auto& names = teqp::GERG2008::component_names;
      const auto comps = shimmer_teqp::utilities::filter_collection(filter,
                                                                    names);
      const auto mol_frac = shimmer_teqp::utilities::filter_collection(filter,
                                                                       x);

      auto model = teqp::GERG2008::GERG2008ResidualModel(comps);

      ASSERT_DOUBLE_EQ_TOL(2.082743529542123e+02,
                           model.red.get_Tr(mol_frac),
                           tolerance);
      ASSERT_DOUBLE_EQ_TOL(1.022348540995680e+04,
                           model.red.get_rhor(mol_frac),
                           tolerance);

      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_H
