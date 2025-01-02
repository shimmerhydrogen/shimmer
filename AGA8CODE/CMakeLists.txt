cmake_minimum_required(VERSION 3.20)

project(SHMGERG LANGUAGES CXX C VERSION 1.0.0)

# Create library variables
###############################################################################################
set(SHMGERG_INSTALL_DIR ${PROJECT_NAME})
set(SHMGERG_INSTALL_BINARY_DIR bin/)
set(SHMGERG_INSTALL_INCLUDE_DIR include/)
set(SHMGERG_INSTALL_LIB_DIR lib/)
set(SHMGERG_INSTALL_ARCHIVE_DIR lib/)
set(SHMGERG_INSTALL_DESTINATION lib/cmake/${SHMGERG_INSTALL_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
Message(STATUS "SHMGERG Build configuration: " ${CMAKE_BUILD_TYPE})


# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${SHMGERG_INSTALL_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
Message(STATUS "SHMGERG Library will be installed in: " ${CMAKE_INSTALL_PREFIX})

option(BUILD_SHMGERG_UNIT_TESTS "Build SHMGERG unitary tests" OFF)

# Add dependencies
###########################################################################################
find_package(Eigen3 REQUIRED)
list(APPEND SHMGERG_LIBRARY_LINKED_LIBRARIES Eigen3::Eigen)

# Insert Sources
###############################################################################################
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
list(APPEND SHMGERG_LIBRARY_SOURCES ${src_sources})
list(APPEND SHMGERG_LIBRARY_HEADERS ${src_headers})
list(APPEND SHMGERG_LIBRARY_INCLUDE ${src_includes})

MESSAGE(WARNING ${SHMGERG_LIBRARY_SOURCES})

# Create public headers
###############################################################################################
set(SHMGERG_LIBRARY_ADDITIONAL_INCLUDE_PATH ${CMAKE_CURRENT_BINARY_DIR}/SHMGERG_additional_include)

foreach(header ${SHMGERG_LIBRARY_HEADERS})
        get_filename_component(headerDirectory ${header} DIRECTORY)
        get_filename_component(headerFileName ${header} NAME_WE)
        get_filename_component(headerExtension ${header} EXT)

  if(${headerExtension} MATCHES ".hpp.in")
    set(headerExtension ".hpp")
    set(headerDirectory ${SHMGERG_LIBRARY_ADDITIONAL_INCLUDE_PATH})

    configure_file(${header} ${headerDirectory}/${headerFileName}${headerExtension} @ONLY)
  elseif(NOT ${headerExtension} MATCHES ".hpp" AND NOT ${headerExtension} MATCHES ".h")
    message(SEND_ERROR "Header extension not recognized for file ${headerFileName}${headerExtension}")
  endif()

  list(APPEND SHMGERG_LIBRARY_PUBLIC_HEADERS ${headerDirectory}/${headerFileName}${headerExtension})
endforeach()

list(APPEND SHMGERG_LIBRARY_INCLUDE ${SHMGERG_LIBRARY_ADDITIONAL_INCLUDE_PATH})

# Create library
###############################################################################################
add_library(${PROJECT_NAME} STATIC ${SHMGERG_LIBRARY_SOURCES} ${SHMGERG_LIBRARY_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
add_library(SHMGERG::SHMGERG ALIAS ${PROJECT_NAME})

target_link_directories(${PROJECT_NAME} PRIVATE ${SHMGERG_LIBRARY_LINKED_DIRECTORIES})
target_link_libraries(${PROJECT_NAME} ${SHMGERG_LIBRARY_LINKED_LIBRARIES})

target_include_directories(${PROJECT_NAME} PRIVATE ${SHMGERG_LIBRARY_INCLUDE})
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${SHMGERG_INSTALL_INCLUDE_DIR}>)
target_compile_definitions(${PROJECT_NAME} PUBLIC ${SHMGERG_compile_definitions})

target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

# Create install target
###############################################################################################
set_target_properties(${PROJECT_NAME} PROPERTIES 
	PUBLIC_HEADER "${SHMGERG_LIBRARY_PUBLIC_HEADERS}"
	EXPORT_NAME SHMGERG
	LINKER_LANGUAGE CXX
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF)

install(
	TARGETS ${PROJECT_NAME} 
        EXPORT SHMGERGTargets
    RUNTIME DESTINATION ${SHMGERG_INSTALL_BINARY_DIR} # runtime executable
    LIBRARY DESTINATION ${SHMGERG_INSTALL_LIB_DIR} # dynamic libraries
    ARCHIVE DESTINATION ${SHMGERG_INSTALL_ARCHIVE_DIR} # static libraries
    PUBLIC_HEADER DESTINATION ${SHMGERG_INSTALL_INCLUDE_DIR} # headers
)


include(CMakePackageConfigHelpers)

configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/SHMGERGConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/SHMGERGConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/${SHMGERG_INSTALL_DESTINATION}
        PATH_VARS SHMGERG_INSTALL_INCLUDE_DIR
  )

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/SHMGERGConfigVersion.cmake
        VERSION ${SHMGERG_VERSION}
	COMPATIBILITY AnyNewerVersion
)

export(
	TARGETS ${PROJECT_NAME} 
        NAMESPACE SHMGERG::
        FILE ${CMAKE_CURRENT_BINARY_DIR}/SHMGERGTargets.cmake
)

install(
        EXPORT SHMGERGTargets
        NAMESPACE SHMGERG::
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${SHMGERG_INSTALL_DESTINATION}
)

install(
	FILES
                ${CMAKE_CURRENT_BINARY_DIR}/SHMGERGConfig.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/SHMGERGConfigVersion.cmake
	DESTINATION
                ${CMAKE_INSTALL_PREFIX}/${SHMGERG_INSTALL_DESTINATION}
)

# Create UnitTests
###############################################################################################
if (BUILD_SHMGERG_UNIT_TESTS)
    include(CTest)
    add_subdirectory(unit_tests)
endif()
