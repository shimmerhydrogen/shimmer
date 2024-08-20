#include <iostream>

#include "test_filter_components.hpp"
#include "test_temp_2.hpp"

int main(int argc, char ** argv)
{
  const std::map<std::string, std::function<int(int , char **)>> tests =
  {
  { "test_filter_components", shimmer_teqp::test::test_filter_components },
  { "test_temp_2", shimmer_teqp::test::test_temp_2 },
};

  if (argc > 1)
  {
    const std::string test_name = argv[1];
    const auto& test_pos = tests.find(test_name);

    if (test_pos == tests.end())
    {
      std::cerr<< "Test "<< test_name<< " not found"<< std::endl;
      return EXIT_FAILURE;
    }

    return test_pos->second(argc, argv);
  }
  else
  {
    for (const auto& test : tests)
    {
      std::cout<< "EXECUTING test "<< test.first<< "..."<< std::endl;
      const auto test_result = test.second(argc, argv);
      if (test_result == EXIT_SUCCESS)
      {
        std::cout<< "\tSUCCESS"<< std::endl;
        continue;
      }

      std::cerr<< "\tFAILED "<< test_result<< std::endl;
    }
  }

  return EXIT_SUCCESS;
}
