#include <iostream>
#include "test_GERG.hpp"
#include "test_GERG_2.hpp"

int main(int argc, char ** argv)
{
  test_GERG(argc, argv);

  if (argc > 1)
  {
    const std::string testName = argv[1];

    if (testName.compare("test_GERG") == 0)
      return test_GERG(argc, argv);

    if (testName.compare("test_GERG_2") == 0)
      return test_GERG_2(argc, argv);
  }

  return EXIT_SUCCESS;
}
