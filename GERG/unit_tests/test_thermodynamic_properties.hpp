#ifndef __test_thermodynamic_properties_H
#define __test_thermodynamic_properties_H

#include <iostream>

#include "Eigen/Eigen"
#include "MATLAB_GERG_functions.hpp"
#include "test_utilities.hpp"

namespace GERG_test
{
  // *********************************************************
  Eigen::MatrixXd thermodynamic_properties_P_mock()
  {
    Eigen::MatrixXd P(3, 1);
    P.col(0)<< 5.101325000000000e+03, 5.101325000000000e+03, 5.101325000000000e+03;

    return P;
  }
  // *********************************************************
  Eigen::MatrixXd thermodynamic_properties_T_mock()
  {
    Eigen::MatrixXd T(3, 1);
    T.col(0)<< 2.931500000000000e+02, 2.931500000000000e+02, 2.931500000000000e+02;

    return T;
  }
  // *********************************************************
  Eigen::MatrixXd thermodynamic_properties_x_mock()
  {
    Eigen::MatrixXd x(3, 21);

    x.row(0)<< 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    x.row(1)<< 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    x.row(2)<< 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    return x;
  }
  // *********************************************************
  GERG::Pseudo_critical_point<Eigen::MatrixXd> thermodynamic_properties_pseudo_critical_point_mock()
  {
    GERG::Pseudo_critical_point<Eigen::MatrixXd> pseudo_critical_point;

    pseudo_critical_point.Tcx.resize(3, 1);
    pseudo_critical_point.Tcx.col(0)<< 1.905640000000000e+02, 1.905640000000000e+02, 1.905640000000000e+02;

    pseudo_critical_point.Dcx.resize(3, 1);
    pseudo_critical_point.Dcx.col(0)<< 1.013934271900000e+01, 1.013934271900000e+01, 1.013934271900000e+01;

    pseudo_critical_point.Vcx.resize(3, 1);
    pseudo_critical_point.Vcx.col(0)<< 9.862572236818776e-02, 9.862572236818776e-02, 9.862572236818776e-02;

    return pseudo_critical_point;
  }
  // *********************************************************
  GERG::Reducing_parameters<Eigen::MatrixXd> thermodynamic_properties_reducing_parameters_mock()
  {
    GERG::Reducing_parameters<Eigen::MatrixXd> reducing_parameters;

    reducing_parameters.Tr.resize(3, 1);
    reducing_parameters.Tr.col(0)<< 1.905640000000000e+02, 1.905640000000000e+02, 1.905640000000000e+02;

    reducing_parameters.Dr.resize(3, 1);
    reducing_parameters.Dr.col(0)<< 1.013934271900000e+01, 1.013934271900000e+01, 1.013934271900000e+01;

    return reducing_parameters;
  }
  // *********************************************************
  GERG::Thermodynamic_properties_parameters thermodynamic_properties_parameters_mock()
  {
    GERG::Thermodynamic_properties_parameters parameters;

    parameters.Type = GERG::Thermodynamic_properties_parameters::Types::Gas_phase;

    return parameters;
  }
  // *********************************************************
  GERG::Thermodynamic_properties<Eigen::MatrixXd> thermodynamic_properties_result_mock()
  {
    GERG::Thermodynamic_properties<Eigen::MatrixXd> result;

    result.D.resize(3, 1);
    result.D.col(0)<< 2.299809852349306e+00, 2.299809852349306e+00, 2.299809852349306e+00;

    result.P1.resize(3, 1);
    result.P1.col(0)<< 5.101324999999993e+03, 5.101324999999993e+03, 5.101324999999993e+03;

    result.Z.resize(3, 1);
    result.Z.col(0)<< 9.100525880503002e-01, 9.100525880503002e-01, 9.100525880503002e-01;

    result.gamma.resize(3, 1);
    result.gamma.col(0)<< 1.303742642844443e+00, 1.303742642844443e+00, 1.303742642844443e+00;

    return result;
  }
  // *********************************************************
  int test_thermodynamic_properties(int , char **)
  {
    const auto P = thermodynamic_properties_P_mock();
    const auto T = thermodynamic_properties_T_mock();
    const auto x = thermodynamic_properties_x_mock();
    const auto reducing_parameters = thermodynamic_properties_reducing_parameters_mock();
    const auto pseudo_critical_point = thermodynamic_properties_pseudo_critical_point_mock();
    const auto parameters = thermodynamic_properties_parameters_mock();

    const auto result = GERG::thermodynamic_properties(P, T,
                                                       x, 3,
                                                       reducing_parameters,
                                                       pseudo_critical_point,
                                                       parameters);

    const auto expected_result = thermodynamic_properties_result_mock();

    ASSERT_MATRIX_DOUBLE_EQ(expected_result.D,
                            result.D);
    ASSERT_MATRIX_DOUBLE_EQ(expected_result.P1,
                            result.P1);
    ASSERT_MATRIX_DOUBLE_EQ(expected_result.Z,
                            result.Z);
    ASSERT_MATRIX_DOUBLE_EQ(expected_result.gamma,
                            result.gamma);

    return EXIT_SUCCESS;
  }
  // *********************************************************
}

#endif // __test_thermodynamic_properties_H
