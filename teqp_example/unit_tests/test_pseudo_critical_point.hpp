#ifndef __test_pseudo_critical_point_H
#define __test_pseudo_critical_point_H

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
    int test_pseudo_critical_point(int , char **)
    {
      using namespace shimmer_teqp::utilities;

      const auto tolerance = gerg_mock::tolerance();
      const auto comps = gerg_mock::comps();
      const auto mol_frac = gerg_mock::mol_frac();

      auto model = teqp::GERG2008::GERG2008ResidualModel(comps);

      const auto expected_pseudo_critical_point = gerg_mock::pseudo_critical_point();

      const auto pseudo_critical_point = shimmer_teqp::utilities::pseudo_critical_point(mol_frac,
                                                                                        model.red.get_Tcvec(),
                                                                                        model.red.get_vcvec(),
                                                                                        gerg_mock::tolerance());

      ASSERT_DOUBLE_EQ_TOL(expected_pseudo_critical_point.Tcx,
                           pseudo_critical_point.Tcx,
                           tolerance);
      ASSERT_DOUBLE_EQ_TOL(expected_pseudo_critical_point.Vcx,
                           pseudo_critical_point.Vcx,
                           tolerance);
      ASSERT_DOUBLE_EQ_TOL(expected_pseudo_critical_point.Dcx,
                           pseudo_critical_point.Dcx,
                           tolerance);

      return EXIT_SUCCESS;
    }
    // *********************************************************
  }
}

#endif // __test_temp_H
