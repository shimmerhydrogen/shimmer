cmake_minimum_required(VERSION 3.9)

project(GERGUnitTest LANGUAGES CXX C)

# Set C++ Standard
###########################################################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## MATLAB
find_package(Matlab REQUIRED)

## Threads
find_package(Threads REQUIRED)

find_package(Eigen3 CONFIG REQUIRED)
list(APPEND GERG_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE Eigen3::Eigen)

## GERG
if (NOT TARGET GERG::GERG)
	find_package(GERG CONFIG REQUIRED)
	list(APPEND GERG_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE GERG::GERG)
else()
	get_target_property(GERG_LIBRARY_HEADERS GERG::GERG PUBLIC_HEADER)
	get_target_property(GERG_LIBRARY_NAME GERG::GERG EXPORT_NAME)
	get_target_property(GERG_LIBRARY_LINKED_LIBRARIES GERG::GERG LINK_LIBRARIES)
	get_target_property(GERG_LIBRARY_INCLUDE_DIRECTORIES GERG::GERG INCLUDE_DIRECTORIES)
	get_target_property(GERG_LIBRARY_SOURCES GERG::GERG SOURCES)

	list(APPEND GERG_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE ${GERG_LIBRARY_NAME} ${GERG_LIBRARY_LINKED_LIBRARIES})
	list(APPEND GERG_LIBRARY_TEST_HEADERS ${GERG_LIBRARY_HEADERS})
	list(APPEND GERG_LIBRARY_TEST_INCLUDE ${GERG_LIBRARY_INCLUDE_DIRECTORIES})
	list(APPEND GERG_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE GERG::GERG)
endif()

# Insert Sources
###############################################################################################
list(APPEND GERG_LIBRARY_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

list(APPEND GERG_LIBRARY_TEST_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities.hpp)

create_test_sourcelist(GERG_LIBRARY_TEST_TESTS
    ${CMAKE_CURRENT_SOURCE_DIR}/test_matlab_interface.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_reducing_parameters.hpp
  )

# Create library test executable
###############################################################################################
add_executable(${PROJECT_NAME}
    ${GERG_LIBRARY_TEST_SOURCES}
    ${GERG_LIBRARY_TEST_HEADERS}
    ${GERG_LIBRARY_TEST_TESTS}
)
target_link_libraries(${PROJECT_NAME} ${GERG_LIBRARY_TEST_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${GERG_LIBRARY_TEST_INCLUDE})

enable_testing()

foreach (test ${GERG_LIBRARY_TEST_TESTS})
  get_filename_component(TName ${test} NAME_WE)
  add_test(NAME ${TName} COMMAND ${PROJECT_NAME} ${TName})
endforeach ()
