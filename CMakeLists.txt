cmake_minimum_required(VERSION 3.18)

project(bcm3)

# Boost settings
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()
set(Boost_USE_MULTITHREADED ON)
if(CMAKE_HOST_WIN32)
    set(Boost_USE_STATIC_LIBS ON)
    set(Boost_USE_STATIC_RUNTIME OFF)
endif(CMAKE_HOST_WIN32)
find_package(Boost REQUIRED COMPONENTS program_options date_time chrono filesystem OPTIONAL_COMPONENTS unit_test_framework)

# NetCDF
find_package(netCDF REQUIRED)

# Some required global sections
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 14)

# Warning handling
if(CMAKE_HOST_WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
endif(CMAKE_HOST_WIN32)	

# Optimization settings
if(CMAKE_HOST_WIN32)
    add_definitions(/arch:AVX2)
endif(CMAKE_HOST_WIN32)
if(CMAKE_HOST_UNIX)
    set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
endif(CMAKE_HOST_UNIX)

add_subdirectory(dependencies/cluster-1.52a/src)
add_subdirectory(dependencies/cvode-5.3.0)
add_subdirectory(dependencies/eigen-3.4-rc1)
add_subdirectory(dependencies/HungarianAlgorithm-master)
add_subdirectory(dependencies/hungarian2)
add_subdirectory(src/bcminf)
add_subdirectory(src/bcmopt)
add_subdirectory(src/bcmrbridge)
add_subdirectory(src/likelihoods)
add_subdirectory(src/odecommon)
add_subdirectory(src/pharmaco)
add_subdirectory(src/sampler)
add_subdirectory(src/stats)
add_subdirectory(src/utils)

# Look for an installed libsbml; otherwise use the included version
find_package(libsbml)
if(NOT libsbml_FOUND)
    message("Using the included libsbml")
    if(CMAKE_HOST_WIN32)
        add_compile_definitions("LIBLAX_STATIC=1;LIBSBML_STATIC=1")
    endif(CMAKE_HOST_WIN32)
    add_subdirectory(dependencies/expat-2.6.0)
    add_subdirectory(dependencies/libsbml-5.20.2)
    set(LIBSBML_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/dependencies/libsbml-5.20.2)
    set(LIBSBML_LIBRARY libsbml)
else()
    message("Using system libsbml")
    set(LIBSBML_LIBRARY libsbml)
endif()

add_subdirectory(src/sbml)
add_subdirectory(src/cellpop)
add_subdirectory(src/fISA)
add_subdirectory(src/testing)

if (Boost_unit_test_framework_FOUND)
    enable_testing()
    add_subdirectory(tests/stats)
else()
    message("Boost unit_test_framework is not available; disabling unit tests")
endif()
