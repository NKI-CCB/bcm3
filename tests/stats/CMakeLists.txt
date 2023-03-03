add_executable(test_stats)

target_sources(test_stats
  PRIVATE
    mvn.cpp
)

target_include_directories(test_stats PRIVATE ${BOOST_INCLUDE_DIRS})

target_include_directories(test_stats
  PUBLIC
    ../../src/stats
)

target_link_libraries(test_stats
  PUBLIC
    stats
)

set_target_properties(test_stats PROPERTIES DEBUG_POSTFIX d)
target_precompile_headers(test_stats PRIVATE ../../src/utils/Utils.h)

add_test(test_stats test_stats)
