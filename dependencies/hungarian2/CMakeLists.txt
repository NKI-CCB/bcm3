add_library(hungarian_matching2)

target_sources(hungarian_matching2
  PRIVATE
    hungarian.h
    hungarian.cpp
)

target_include_directories(hungarian_matching2
  PUBLIC
    .
  PRIVATE
    ../../src/utils
)

target_link_libraries(hungarian_matching2
  PRIVATE
    utils
)

set_target_properties(hungarian_matching2 PROPERTIES DEBUG_POSTFIX d)
set_target_properties(hungarian_matching2 PROPERTIES FOLDER Dependencies)
