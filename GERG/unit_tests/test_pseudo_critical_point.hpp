#ifndef __test_pseudo_critical_point_H
#define __test_pseudo_critical_point_H

#include <iostream>

#include "Eigen/Eigen"
#include "GERG_functions.hpp"
#include "MATLAB_GERG_functions.hpp"
#include "test_utilities.hpp"

namespace GERG_test
{
  // *********************************************************
  Eigen::MatrixXd pseudo_critical_point_x_mock()
  {
    Eigen::MatrixXd x(3, 21);

    x.row(0)<< 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    x.row(1)<< 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    x.row(2)<< 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    return x;
  }
  // *********************************************************
  GERG::Pseudo_critical_point<Eigen::MatrixXd> pseudo_critical_point_result_mock()
  {
    GERG::Pseudo_critical_point<Eigen::MatrixXd> result;

    result.Tcx.resize(3, 1);
    result.Tcx.col(0)<< 1.905640000000000e+02, 1.905640000000000e+02, 1.905640000000000e+02;

    result.Dcx.resize(3, 1);
    result.Dcx.col(0)<< 1.013934271900000e+01, 1.013934271900000e+01, 1.013934271900000e+01;

    result.Vcx.resize(3, 1);
    result.Vcx.col(0)<< 9.862572236818776e-02, 9.862572236818776e-02, 9.862572236818776e-02;

    return result;
  }
  // *********************************************************
  int test_pseudo_critical_point(int , char **)
  {
    const auto x = pseudo_critical_point_x_mock();
    const auto result = GERG::pseudo_critical_point(x, 3);
    const auto expected_result = pseudo_critical_point_result_mock();

    ASSERT_MATRIX_DOUBLE_EQ(expected_result.Tcx,
                            result.Tcx);
    ASSERT_MATRIX_DOUBLE_EQ(expected_result.Dcx,
                            result.Dcx);
    ASSERT_MATRIX_DOUBLE_EQ(expected_result.Vcx,
                            result.Vcx);

    return EXIT_SUCCESS;
  }
  // *********************************************************
}

#endif // __test_pseudo_critical_point_H
