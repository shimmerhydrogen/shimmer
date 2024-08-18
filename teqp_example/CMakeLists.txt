cmake_minimum_required(VERSION 3.20)

project(SHMTEQP LANGUAGES CXX C VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create library variables
###############################################################################################
set(SHMTEQP_INSTALL_DIR ${PROJECT_NAME})
set(SHMTEQP_INSTALL_BINARY_DIR bin/)
set(SHMTEQP_INSTALL_INCLUDE_DIR include/)
set(SHMTEQP_INSTALL_LIB_DIR lib/)
set(SHMTEQP_INSTALL_ARCHIVE_DIR lib/)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
Message(STATUS "SHMTEQP Build configuration: " ${CMAKE_BUILD_TYPE})

if(NOT TEQP_DIR)
  set(TEQP_DIR "/home/geoscore/Desktop/GEO++/shimmer/teqp" CACHE STRING "Find Teqp directory" FORCE) # TODO: remove this and put error message!
endif()
Message(STATUS "SHMTEQP teqp folder: " ${TEQP_DIR})

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${SHMTEQP_INSTALL_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
Message(STATUS "SHMTEQP Library will be installed in: " ${CMAKE_INSTALL_PREFIX})

option(BUILD_SHMTEQP_UNIT_TESTS "Build SHMTEQP unitary tests" OFF)

list(APPEND CMAKE_PREFIX_PATH "${TEQP_DIR}/build/install")

# Add dependencies
###########################################################################################
list(APPEND SHMTEQP_LIBRARY_LINKED_LIBRARIES teqpcpp)
list(APPEND SHMTEQP_LIBRARY_LINKED_DIRECTORIES "${TEQP_DIR}/build")
list(APPEND SHMTEQP_LIBRARY_LINKED_DIRECTORIES "${TEQP_DIR}/build/install/lib")
list(APPEND SHMTEQP_LIBRARY_INCLUDE "${TEQP_DIR}/include")
list(APPEND SHMTEQP_LIBRARY_INCLUDE "${TEQP_DIR}/install/include")

find_package(Boost REQUIRED)

find_package(Eigen3 REQUIRED)
list(APPEND SHMTEQP_LIBRARY_LINKED_LIBRARIES Eigen3::Eigen)

find_package(autodiff REQUIRED)
list(APPEND SHMTEQP_LIBRARY_LINKED_LIBRARIES autodiff::autodiff)

find_package(nlohmann_json_schema_validator REQUIRED)
list(APPEND SHMTEQP_LIBRARY_LINKED_LIBRARIES nlohmann_json_schema_validator)

# Insert Sources
###############################################################################################
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
list(APPEND SHMTEQP_LIBRARY_SOURCES ${src_sources})
list(APPEND SHMTEQP_LIBRARY_HEADERS ${src_headers})
list(APPEND SHMTEQP_LIBRARY_INCLUDE ${src_includes})

# Create public headers
###############################################################################################
set(SHMTEQP_LIBRARY_ADDITIONAL_INCLUDE_PATH ${CMAKE_CURRENT_BINARY_DIR}/SHMTEQP_additional_include)

foreach(header ${SHMTEQP_LIBRARY_HEADERS})
        get_filename_component(headerDirectory ${header} DIRECTORY)
        get_filename_component(headerFileName ${header} NAME_WE)
        get_filename_component(headerExtension ${header} EXT)

  if(${headerExtension} MATCHES ".hpp.in")
    set(headerExtension ".hpp")
    set(headerDirectory ${SHMTEQP_LIBRARY_ADDITIONAL_INCLUDE_PATH})

    configure_file(${header} ${headerDirectory}/${headerFileName}${headerExtension} @ONLY)
  elseif(NOT ${headerExtension} MATCHES ".hpp" AND NOT ${headerExtension} MATCHES ".h")
    message(SEND_ERROR "Header extension not recognized for file ${headerFileName}${headerExtension}")
  endif()

  list(APPEND SHMTEQP_LIBRARY_PUBLIC_HEADERS ${headerDirectory}/${headerFileName}${headerExtension})
endforeach()

list(APPEND SHMTEQP_LIBRARY_INCLUDE ${SHMTEQP_LIBRARY_ADDITIONAL_INCLUDE_PATH})

# Create library
###############################################################################################
add_library(${PROJECT_NAME} STATIC ${SHMTEQP_LIBRARY_SOURCES} ${SHMTEQP_LIBRARY_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
add_library(SHMTEQP::SHMTEQP ALIAS ${PROJECT_NAME})

target_link_directories(${PROJECT_NAME} PRIVATE ${SHMTEQP_LIBRARY_LINKED_DIRECTORIES})
target_link_libraries(${PROJECT_NAME} ${SHMTEQP_LIBRARY_LINKED_LIBRARIES})

target_include_directories(${PROJECT_NAME} PRIVATE ${SHMTEQP_LIBRARY_INCLUDE})
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${SHMTEQP_INSTALL_INCLUDE_DIR}>)
target_compile_definitions(${PROJECT_NAME} PUBLIC ${SHMTEQP_compile_definitions})

target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

# Create install target
###############################################################################################
set_target_properties(${PROJECT_NAME} PROPERTIES 
PUBLIC_HEADER "${SHMTEQP_LIBRARY_PUBLIC_HEADERS}"
EXPORT_NAME SHMTEQP)

install(
	TARGETS ${PROJECT_NAME} 
        EXPORT SHMTEQPTargets
    RUNTIME DESTINATION ${SHMTEQP_INSTALL_BINARY_DIR} # runtime executable
    LIBRARY DESTINATION ${SHMTEQP_INSTALL_LIB_DIR} # dynamic libraries
    ARCHIVE DESTINATION ${SHMTEQP_INSTALL_ARCHIVE_DIR} # static libraries
    PUBLIC_HEADER DESTINATION ${SHMTEQP_INSTALL_INCLUDE_DIR} # headers
)


include(CMakePackageConfigHelpers)

configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/SHMTEQPConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/SHMTEQPConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/${SHMTEQP_INSTALL_LIB_DIR}/cmake/${SHMTEQP_INSTALL_DIR}
        PATH_VARS SHMTEQP_INSTALL_INCLUDE_DIR
  )

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/SHMTEQPConfigVersion.cmake
        VERSION ${SHMTEQP_VERSION}
	COMPATIBILITY AnyNewerVersion
)

export(
	TARGETS ${PROJECT_NAME} 
        NAMESPACE SHMTEQP::
        FILE ${CMAKE_CURRENT_BINARY_DIR}/SHMTEQPTargets.cmake
)

install(
        EXPORT SHMTEQPTargets
        NAMESPACE SHMTEQP::
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${SHMTEQP_INSTALL_LIB_DIR}/cmake/${SHMTEQP_INSTALL_DIR}
)

install(
	FILES
                ${CMAKE_CURRENT_BINARY_DIR}/SHMTEQPConfig.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/SHMTEQPConfigVersion.cmake
	DESTINATION
                ${CMAKE_INSTALL_PREFIX}/${SHMTEQP_INSTALL_LIB_DIR}/cmake/${SHMTEQP_INSTALL_DIR}
)

# Create UnitTests
###############################################################################################
if (BUILD_SHMTEQP_UNIT_TESTS)
    include(CTest)
    add_subdirectory(unit_tests)
endif()

#add_executable(hello main.cpp shimmer_teqp_utilities.hpp)
#target_link_libraries(hello PRIVATE teqpcpp PRIVATE Eigen3::Eigen PRIVATE autodiff::autodiff PRIVATE nlohmann_json_schema_validator)
#target_link_directories(hello PRIVATE "${TEQP_DIR}/build")
#target_link_directories(hello PRIVATE "${TEQP_DIR}/build/install/lib")
#target_include_directories(hello PRIVATE "${TEQP_DIR}/include")
#target_include_directories(hello PRIVATE "${TEQP_DIR}/install/include")
