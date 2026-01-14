add_library(hungarian_matching)

target_sources(hungarian_matching
  PRIVATE
    matching.h
    matching.cpp
)

target_include_directories(hungarian_matching
  PUBLIC
    .
  PRIVATE
    ../../src/utils
)

target_link_libraries(hungarian_matching
  PRIVATE
    utils
)

set_target_properties(hungarian_matching PROPERTIES DEBUG_POSTFIX d)
set_target_properties(hungarian_matching PROPERTIES FOLDER Dependencies)
