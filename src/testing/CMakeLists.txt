add_executable(testing)

target_sources(testing
  PRIVATE
    main.cpp
)

target_link_libraries(testing
  PRIVATE
    utils
    hungarian
    hungarian2
)

set_target_properties(testing PROPERTIES DEBUG_POSTFIX d)
set_target_properties(testing PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../../../bin)