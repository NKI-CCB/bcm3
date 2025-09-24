add_library(utils)

target_sources(utils
  PRIVATE
    Clustering.h
    Clustering.cpp
    Correlation.h
    Correlation.cpp
    CSVParser.h
    CSVParser.cpp
    EigenPartialPivLUSomewhatSparse.h
    Logger.h
    Logger.cpp
    MathFunctions.h
    MathFunctions.cpp
    NetCDFBundler.h
    NetCDFBundler.cpp
    NetCDFDataFile.h
    NetCDFDataFile.cpp
    Platform.h
    ProbabilityDistributions.h
    ProbabilityDistributions.cpp
    RNG.h
    RNG.cpp
    Spinlock.h
    SummaryStats.h
    SummaryStats.cpp
    TaskManager.h
    TaskManager.cpp
    Timer.h
    Typedefs.h
    Utils.h
    Utils.cpp
    VectorUtils.h
    VectorUtils.cpp
    Version.h
)

target_include_directories(utils
  PUBLIC
    .
  PRIVATE
    ${NETCDF_INCLUDE_DIR}
)

target_link_libraries(utils
  PUBLIC
    Eigen
    Boost::headers
    ${Boost_LIBRARIES}
    Threads::Threads
  PRIVATE
    netCDF::netcdf
    clusterc
)

set_target_properties(utils PROPERTIES DEBUG_POSTFIX d)
target_precompile_headers(utils PRIVATE Utils.h)
