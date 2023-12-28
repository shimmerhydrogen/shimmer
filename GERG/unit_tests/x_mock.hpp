#ifndef __x_mock_H
#define __x_mock_H

#include "Eigen/Eigen"

namespace GERG_test
{
  Eigen::MatrixXd x_mock()
  {
    Eigen::MatrixXd x(3, 21);

    return x;
  }
}

#endif // __x_mock_H
