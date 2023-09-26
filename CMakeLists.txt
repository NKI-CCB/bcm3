cmake_minimum_required(VERSION 3.16)
project(bcm3)

# Boost settings
set(Boost_USE_MULTITHREADED		ON)
if(CMAKE_HOST_WIN32)
    set(Boost_USE_STATIC_LIBS		ON)
    set(Boost_USE_STATIC_RUNTIME	OFF)
endif(CMAKE_HOST_WIN32)
find_package(Boost REQUIRED COMPONENTS program_options filesystem system date_time chrono unit_test_framework)



# bthijssen: I tried to make these work with find_package, but didn't manage to find a way that works
# on Windows with either prebuilt netCDF/libsbml or custom ones. So for now it is required to manually
# specify the location of these libraries
include("external_dependency_locations.txt")

# Some required global sections
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)

find_package(Torch REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

# Optimization settings
if(CMAKE_HOST_WIN32)
    add_definitions(/arch:AVX)
endif(CMAKE_HOST_WIN32)
if(CMAKE_HOST_UNIX)
    set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
    set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif(CMAKE_HOST_UNIX)

add_subdirectory(dependencies/cluster-1.52a/src)
add_subdirectory(dependencies/cvode-5.3.0)
add_subdirectory(dependencies/eigen-3.4-rc1)
add_subdirectory(dependencies/HungarianAlgorithm-master)
add_subdirectory(src/bcminf)
add_subdirectory(src/bcmopt)
add_subdirectory(src/bcmrbridge)
add_subdirectory(src/cellpop)
add_subdirectory(src/fISA)
add_subdirectory(src/likelihoods)
add_subdirectory(src/odecommon)
add_subdirectory(src/sampler)
add_subdirectory(src/sbml)
add_subdirectory(src/stats)
add_subdirectory(src/utils)

enable_testing()

add_subdirectory(tests/stats)