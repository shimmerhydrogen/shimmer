cmake_minimum_required(VERSION 3.20)

project(hello LANGUAGES CXX C VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 20)

set(TEQP_DIR "/home/geoscore/Desktop/GEO++/shimmer/teqp")

list(APPEND CMAKE_PREFIX_PATH "${TEQP_DIR}/build/install")

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(autodiff REQUIRED)
find_package(nlohmann_json_schema_validator REQUIRED)

add_executable(hello main.cpp)
target_link_libraries(hello PRIVATE teqpcpp PRIVATE Eigen3::Eigen PRIVATE autodiff::autodiff PRIVATE nlohmann_json_schema_validator)
target_link_directories(hello PRIVATE "${TEQP_DIR}/build")
target_link_directories(hello PRIVATE "${TEQP_DIR}/build/install/lib")
target_include_directories(hello PRIVATE "${TEQP_DIR}/include")
target_include_directories(hello PRIVATE "${TEQP_DIR}/install/include")
#target_include_directories(hello PRIVATE "${TEQP_DIR}/externals/nlohmann_json")
#target_include_directories(hello PRIVATE "${TEQP_DIR}/externals/autodiff")
