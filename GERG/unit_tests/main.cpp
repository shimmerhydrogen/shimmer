#include <iostream>
#include "test_GERG.hpp"
#include "test_reducing_parameters.hpp"

int main(int argc, char ** argv)
{
  GERG_test::test_reducing_parameters(argc, argv);

  if (argc > 1)
  {
    const std::string testName = argv[1];

    if (testName.compare("test_GERG") == 0)
      return test_GERG(argc, argv);

    if (testName.compare("test_reducing_parameters") == 0)
      return GERG_test::test_reducing_parameters(argc, argv);
  }

  return EXIT_SUCCESS;
}
