#include <iostream>

#include "test_filter_components.hpp"
#include "test_filter_collection.hpp"
#include "test_reducing_parameters.hpp"
#include "test_temp_2.hpp"

int main(int argc, char ** argv)
{
  const std::map<std::string, std::function<int(int , char **)>> tests =
  {
  { "test_filter_collection", shimmer_teqp::test::test_filter_collection },
  { "test_filter_components", shimmer_teqp::test::test_filter_components },
  { "test_reducing_parameters", shimmer_teqp::test::test_reducing_parameters },
  { "test_temp_2", shimmer_teqp::test::test_temp_2 },
};

  std::list<std::string> tests_to_run;
  if (argc > 1)
  {
    const unsigned int num_tests_to_run = argc - 1;
    for (unsigned int t = 0; t < num_tests_to_run; t++)
      tests_to_run.push_back(argv[t + 1]);
  }
  else
  {
    for (const auto& test : tests)
      tests_to_run.push_back(test.first);
  }

  const unsigned int num_tests = tests_to_run.size();

  if (num_tests == 0)
    return EXIT_SUCCESS;

  unsigned int test_number = 1;
  unsigned int num_tests_success = 0;

  std::cout<< "Executing "<< num_tests<< " tests ..."<< std::endl;
  std::cout<< std::endl;

  const double start_total_time = clock();
  for (const auto& test_name : tests_to_run)
  {
    const auto test_pos = tests.find(test_name);
    if (test_pos == tests.end())
      continue;

    std::cout<< "    Start test "<< test_number<< ": "<< test_name<< std::endl;

    const double start_time = clock();
    const auto test_result = test_pos->second(argc, argv);
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
    std::cout<< "Test #"<< test_number<< ": "<< test_name<< " ...\t";
    std::cout<< test_result_message<< "\t";
    std::cout<< std::scientific<< test_time<< " sec"<< std::endl;

    test_number++;
  }
  const double end_total_time = clock();
  const double total_tests_time = (end_total_time - start_total_time) / static_cast<double>(CLOCKS_PER_SEC);

  std::cout<< std::endl;
  std::cout.precision(1);
  std::cout<< std::fixed<< "Passed "<< 100.0 * (num_tests_success / num_tests)<< " tests"<< " ";
  std::cout<< (num_tests - num_tests_success)<< " test failed out of "<< num_tests<< std::endl;
  std::cout.precision(1);
  std::cout<< std::scientific<< "Total Test time (real) = "<< std::scientific<< total_tests_time<< " sec"<< std::endl;


  return EXIT_SUCCESS;
}
