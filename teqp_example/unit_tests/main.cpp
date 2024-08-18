#include <iostream>

#include "test_temp.hpp"

int main(int argc, char ** argv)
{
  if (argc > 1)
  {
    const std::string testName = argv[1];

    if (testName.compare("test_temp") == 0)
      return shimmer_teqp::test::test_temp(argc, argv);
  }

  return EXIT_SUCCESS;
}
