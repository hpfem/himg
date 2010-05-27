set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cxx_flag_overrides.cmake)
project(himg)

cmake_minimum_required(VERSION 2.6)

#---------------------------------------------------------
# begin of configuration section
# The values below are sane defaults to produce a build on all platforms,
# if you want to override them, set them in CMake.vars

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# build with debugging symbols
set(DEBUG YES)
set(RELEASE NO)

set(WITH_JPEG NO)

# allow to override the default values in CMake.vars
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
	include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

# end of configuration section
#---------------------------------------------------------

find_package(HERMES_COMMON REQUIRED)
find_package(HERMES2D REQUIRED)
find_package(UMFPACK REQUIRED)
find_package(LIBPNG REQUIRED)
if(WITH_JPEG)
    message("JPEG: No implemented yet")
    find_package(LIBJPEG REQUIRED)
endif(WITH_JPEG)

add_subdirectory(src)
