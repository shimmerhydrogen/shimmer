#ifndef __TEST_UTILITIES_HPP__
#define __TEST_UTILITIES_HPP__

#include <iostream>

#define MIN_DOUBLE_TOL 2.0 * std::numeric_limits<double>::epsilon()

#define ASSERT_MSG(exp, msg) \
  if (!(exp)) { std::cerr<< msg<< std::endl; \
  abort(); }

#define ASSERT_TRUE_MSG(value, msg) \
  ASSERT_MSG(value, msg)

#define ASSERT_TRUE(value) \
  ASSERT_TRUE_MSG(value, "Expected true value")

#define ASSERT_EQ_MSG(expected, actual, msg) \
  ASSERT_MSG(expected == actual, msg)

#define ASSERT_EQ(expected, actual) \
  std::stringstream msg; \
  msg<< "Expected "<< expected<< " equal to "<< actual<< std::endl; \
  ASSERT_EQ_MSG(expected, actual, msg.str())

#define ASSERT_LE_MSG(expected, actual, msg) \
  ASSERT_MSG(expected < actual || expected == actual, msg)

#define ASSERT_LE(expected, actual) \
  std::stringstream msg; \
  msg<< "Expected "<< expected<< " lower or equal to "<< actual<< std::endl; \
  ASSERT_LE_MSG(expected, actual, msg.str())

#define ASSERT_DOUBLE_EQ_TOL_MSG(expected, actual, tolerance, msg) \
  ASSERT_LE_MSG(std::abs(expected - actual), tolerance * (std::abs(expected) <= tolerance ? 1.0 : std::abs(expected)), msg)

#define ASSERT_DOUBLE_EQ_TOL(expected, actual, tolerance) \
  std::stringstream msg; \
  msg.precision(16); \
  msg<< std::scientific<< "Expected "<< expected<< " equal to "<< actual<< " respect tol "<< tolerance<< std::endl; \
  ASSERT_DOUBLE_EQ_TOL_MSG(expected, actual, tolerance, msg.str())

#define ASSERT_DOUBLE_EQ(expected, actual) \
  ASSERT_DOUBLE_EQ_TOL(expected, actual, MIN_DOUBLE_TOL)

#define ASSERT_VECTOR_DOUBLE_EQ_TOL(expected, actual, tolerance) \
  ASSERT_EQ_MSG(expected.size(), actual.size(), "Array sizes differ") \
  for (unsigned int idx = 0; idx < expected.size(); idx++) { \
  std::stringstream msg; \
  msg.precision(16); \
  msg<< std::scientific<< "Expected "<< expected[idx]<< " equal to "<< actual[idx]<< " respect tol "<< tolerance<< "at index: " << idx<< std::endl; \
  ASSERT_DOUBLE_EQ_TOL_MSG(expected[idx], actual[idx], tolerance, msg.str()) }

#define ASSERT_VECTOR_DOUBLE_EQ(expected, actual) \
  ASSERT_VECTOR_DOUBLE_EQ_TOL(expected, actual, MIN_DOUBLE_TOL)

#define ASSERT_MATRIX_DOUBLE_EQ_TOL(expected, actual, tolerance) \
  ASSERT_EQ_MSG(expected.rows(), actual.rows(), "Matrix rows differ. Expect " + std::to_string(expected.rows()) + " actual " + std::to_string(actual.rows())) \
  ASSERT_EQ_MSG(expected.cols(), actual.cols(), "Matrix cols differ. Expect " + std::to_string(expected.rows()) + " actual " + std::to_string(actual.rows())) \
  for (unsigned int r = 0; r < expected.rows(); r++) { \
  for (unsigned int c = 0; c < expected.cols(); c++) { \
  std::stringstream msg; \
  msg.precision(16); \
  msg<< std::scientific<< "Expected "<< expected(r, c)<< " equal to "<< actual(r, c)<< " respect tol "<< tolerance<< " at position: (" << r<< ","<< c<< ")"<< std::endl; \
  ASSERT_DOUBLE_EQ_TOL_MSG(expected(r, c), actual(r, c), tolerance, msg.str()) } }

#define ASSERT_MATRIX_DOUBLE_EQ(expected, actual) \
  ASSERT_MATRIX_DOUBLE_EQ_TOL(expected, actual, MIN_DOUBLE_TOL)

#endif // __TEST_UTILITIES_HPP__
