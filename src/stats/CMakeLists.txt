add_library(stats)

target_sources(stats
  PRIVATE
    Clustering.h
    Clustering.cpp
    SummaryStats.h
    SummaryStats.cpp
)

target_include_directories(stats
  PUBLIC
    ../utils
)

target_link_libraries(stats
  PUBLIC
    utils
)

set_target_properties(stats PROPERTIES DEBUG_POSTFIX d)
target_precompile_headers(stats PRIVATE ../utils/Utils.h)