#ifndef __test_reducing_parameters_H
#define __test_reducing_parameters_H

#include "Eigen/Eigen"
#include "shimmer_teqp_utilities.hpp"

#include "teqp/models/GERG/GERG.hpp"
#include "test_gerg_mock.hpp"
#include "test_utilities.hpp"

namespace shimmer_teqp
{
  namespace test
  {
    // *********************************************************
    int test_reducing_parameters(int , char **)
    {
      using namespace shimmer_teqp::utilities;

      const auto tolerance = gerg_mock::tolerance();
      const auto comps = gerg_mock::comps();
      const auto mol_frac = gerg_mock::mol_frac();

      auto model = teqp::GERG2008::GERG2008ResidualModel(comps);

      const auto expected_reducing_parameters = gerg_mock::reducing_parameters();

      ASSERT_DOUBLE_EQ_TOL(expected_reducing_parameters.Tr,
                           model.red.get_Tr(mol_frac),
                           tolerance);
      ASSERT_DOUBLE_EQ_TOL(expected_reducing_parameters.Dr,
                           model.red.get_rhor(mol_frac),
                           tolerance);

      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_H
