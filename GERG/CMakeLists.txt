cmake_minimum_required(VERSION 3.9)

project(GERG LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
###########################################################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create library variables
###############################################################################################
set(GERG_INSTALL_DIR ${PROJECT_NAME})
set(GERG_INSTALL_BINARY_DIR bin/)
set(GERG_INSTALL_INCLUDE_DIR include/)
set(GERG_INSTALL_LIB_DIR lib/)
set(GERG_INSTALL_ARCHIVE_DIR lib/)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
Message(STATUS "GERG Build configuration: " ${CMAKE_BUILD_TYPE})

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${GERG_INSTALL_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
Message(STATUS "GERG Library will be installed in: " ${CMAKE_INSTALL_PREFIX})

option(BUILD_GERG_UNIT_TESTS "Build GERG unitary tests" OFF)

# Add dependencies
###########################################################################################
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

## MATLAB
find_package(Matlab REQUIRED)
list(APPEND GERG_LIBRARY_LINKED_LIBRARIES Matlab::MatlabEngine)
list(APPEND GERG_LIBRARY_LINKED_LIBRARIES Matlab::MatlabDataArray)

## Threads
#find_package(Threads REQUIRED)
#list(APPEND GERG_LIBRARY_LINKED_LIBRARIES Threads::Threads)

# Insert Sources
###############################################################################################
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
list(APPEND GERG_LIBRARY_SOURCES ${gerg_src_sources})
list(APPEND GERG_LIBRARY_HEADERS ${gerg_src_headers})
list(APPEND GERG_LIBRARY_INCLUDE ${gerg_src_includes})

# Create public headers
###############################################################################################
set(GERG_LIBRARY_ADDITIONAL_INCLUDE_PATH ${CMAKE_CURRENT_BINARY_DIR}/GERG_additional_include)

foreach(header ${GERG_LIBRARY_HEADERS})
        get_filename_component(headerDirectory ${header} DIRECTORY)
        get_filename_component(headerFileName ${header} NAME_WE)
        get_filename_component(headerExtension ${header} EXT)

  if(${headerExtension} MATCHES ".hpp.in")
    set(headerExtension ".hpp")
    set(headerDirectory ${GERG_LIBRARY_ADDITIONAL_INCLUDE_PATH})

    configure_file(${header} ${headerDirectory}/${headerFileName}${headerExtension} @ONLY)
  elseif(NOT ${headerExtension} MATCHES ".hpp" AND NOT ${headerExtension} MATCHES ".h")
    message(SEND_ERROR "Header extension not recognized for file ${headerFileName}${headerExtension}")
  endif()

  list(APPEND GERG_LIBRARY_PUBLIC_HEADERS ${headerDirectory}/${headerFileName}${headerExtension})
endforeach()

list(APPEND GERG_LIBRARY_INCLUDE ${GERG_LIBRARY_ADDITIONAL_INCLUDE_PATH})

# Create library
###############################################################################################
add_library(${PROJECT_NAME} STATIC ${GERG_LIBRARY_SOURCES} ${GERG_LIBRARY_HEADERS})
add_library(GERG::GERG ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} ${GERG_LIBRARY_LINKED_LIBRARIES})

target_include_directories(${PROJECT_NAME} PRIVATE ${GERG_LIBRARY_INCLUDE})
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${GERG_INSTALL_INCLUDE_DIR}>)
target_compile_definitions(${PROJECT_NAME} PUBLIC ${GERG_compile_definitions})

target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

# Create install target
###############################################################################################
set_target_properties(${PROJECT_NAME} PROPERTIES 
PUBLIC_HEADER "${GERG_LIBRARY_PUBLIC_HEADERS}" 
EXPORT_NAME GERG)

install(
	TARGETS ${PROJECT_NAME} 
	EXPORT GERGTargets
    RUNTIME DESTINATION ${GERG_INSTALL_BINARY_DIR} # runtime executable
    LIBRARY DESTINATION ${GERG_INSTALL_LIB_DIR} # dynamic libraries
    ARCHIVE DESTINATION ${GERG_INSTALL_ARCHIVE_DIR} # static libraries
    PUBLIC_HEADER DESTINATION ${GERG_INSTALL_INCLUDE_DIR} # headers
)


include(CMakePackageConfigHelpers)

configure_package_config_file(
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/GERGConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/GERGConfig.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/${GERG_INSTALL_LIB_DIR}/cmake/${GERG_INSTALL_DIR}
	PATH_VARS GERG_INSTALL_INCLUDE_DIR
  )

write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/GERGConfigVersion.cmake
	VERSION ${GERG_VERSION}
	COMPATIBILITY AnyNewerVersion
)

export(
	TARGETS ${PROJECT_NAME} 
	NAMESPACE GERG::
	FILE ${CMAKE_CURRENT_BINARY_DIR}/GERGTargets.cmake
)

install(
	EXPORT GERGTargets 
	NAMESPACE GERG:: 
	DESTINATION ${CMAKE_INSTALL_PREFIX}/${GERG_INSTALL_LIB_DIR}/cmake/${GERG_INSTALL_DIR}
)

install(
	FILES
		${CMAKE_CURRENT_BINARY_DIR}/GERGConfig.cmake
		${CMAKE_CURRENT_BINARY_DIR}/GERGConfigVersion.cmake
	DESTINATION
		${CMAKE_INSTALL_PREFIX}/${GERG_INSTALL_LIB_DIR}/cmake/${GERG_INSTALL_DIR}
)

# Create UnitTests
###############################################################################################
if (BUILD_GERG_UNIT_TESTS)
    include(CTest)
    add_subdirectory(unit_tests)
endif()
