#include <iostream>

#include "test_filter_components.hpp"
#include "test_filter_collection.hpp"
#include "test_temp_2.hpp"

int main(int argc, char ** argv)
{
  const std::map<std::string, std::function<int(int , char **)>> tests =
  {
  { "test_filter_collection", shimmer_teqp::test::test_filter_collection },
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
    const unsigned int num_tests = tests.size();
    unsigned int test_number = 1;
    unsigned int num_tests_success = 0;

    std::cout<< "Executing "<< num_tests<< " tests ..."<< std::endl;
    std::cout<< std::endl;

    for (const auto& test : tests)
    {
      std::cout<< "\tStart test "<< test_number<< ": "<< test.first<< std::endl;

      const double start_time = clock();
      const auto test_result = test.second(argc, argv);
      const double end_time = clock();

      if (test_result == EXIT_SUCCESS)
        num_tests_success++;

      const std::string test_result_message =
          (test_result == EXIT_SUCCESS) ?
            "Passed" :
            "Failed - " + std::to_string(test_result);
      const double test_time = (end_time - start_time) / static_cast<double>(CLOCKS_PER_SEC);

      std::cout.precision(2);
      std::cout<< test_number<< "/"<< num_tests<< " ";
      std::cout<< "Test #"<< test_number<< ": "<< test.first<< " ... ";
      std::cout<< test_result_message<< "\t";
      std::cout<< std::scientific<< test_time<< " sec"<< std::endl;

      test_number++;
    }

    std::cout<< std::endl;
    std::cout.precision(1);
    std::cout<< std::fixed<< "Passed "<< 100.0 * (num_tests_success / num_tests)<< " tests"<< " ";
    std::cout<< (num_tests - num_tests_success)<< " test failed out of "<< num_tests<< std::endl;
  }

  return EXIT_SUCCESS;
}
