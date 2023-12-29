#include <iostream>
#include "test_matlab_interface.hpp"
#include "test_pseudo_critical_point.hpp"
#include "test_reducing_parameters.hpp"
#include "test_thermodynamic_properties.hpp"

int main(int argc, char ** argv)
{
  if (argc > 1)
  {
    const std::string testName = argv[1];

    if (testName.compare("test_matlab_interface") == 0)
      return GERG_test::test_matlab_interface(argc, argv);

    if (testName.compare("test_reducing_parameters") == 0)
      return GERG_test::test_reducing_parameters(argc, argv);

    if (testName.compare("test_pseudo_critical_point") == 0)
      return GERG_test::test_pseudo_critical_point(argc, argv);

    if (testName.compare("test_thermodynamic_properties") == 0)
      return GERG_test::test_thermodynamic_properties(argc, argv);
  }

  return EXIT_SUCCESS;
}
