add_library(stats)

target_sources(stats
  PRIVATE
    checks.h
    checks.cpp
    mvn.h
    mvn.cpp
    mvt.h
    mvt.cpp
)

target_include_directories(stats
  PUBLIC
    .
)

target_link_libraries(stats
  PUBLIC
    utils
)

set_target_properties(stats PROPERTIES DEBUG_POSTFIX d)
target_precompile_headers(stats PRIVATE ../utils/Utils.h)