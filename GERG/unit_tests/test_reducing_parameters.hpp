#ifndef __test_reducing_parameters_H
#define __test_reducing_parameters_H

#include <iostream>

#include "Eigen/Eigen"
#include "MATLAB_GERG_functions.hpp"
#include "test_utilities.hpp"

namespace GERG_test
{
  Eigen::MatrixXd reducing_parameters_x_mock()
  {
    Eigen::MatrixXd x(3, 21);

    x.row(0)<< 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    x.row(1)<< 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    x.row(2)<< 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    return x;
  }

  GERG::Reducing_parameters<Eigen::MatrixXd> reducing_parameters_result_mock()
  {
    GERG::Reducing_parameters<Eigen::MatrixXd> result;

    result.Tr.resize(3, 1);
    result.Tr.col(0)<< 1.905640000000000e+02, 1.905640000000000e+02, 1.905640000000000e+02;

    result.Dr.resize(3, 1);
    result.Dr.col(0)<< 1.013934271900000e+01, 1.013934271900000e+01, 1.013934271900000e+01;

    return result;
  }

  int test_reducing_parameters(int , char **)
  {
    const auto x = reducing_parameters_x_mock();
    const auto result = GERG::reducing_parameters(x);
    const auto expected_result = reducing_parameters_result_mock();

    ASSERT_EQ(1.0, 2.0);

    ASSERT_MATRIX_DOUBLE_EQ(expected_result.Tr,
                            result.Tr);
    ASSERT_MATRIX_DOUBLE_EQ(expected_result.Dr,
                            result.Dr);

    std::cout<< "CALL test_reducing_parameters"<< std::endl;


    return EXIT_SUCCESS;
  }
}

#endif // __test_reducing_parameters_H
