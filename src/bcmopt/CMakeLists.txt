add_executable(bcmopt)

target_sources(bcmopt
  PRIVATE
    main.cpp
)

target_include_directories(bcmopt
  PRIVATE
    ../likelihoods
)

target_link_libraries(bcmopt
  PRIVATE
    likelihoods
)

set_target_properties(bcmopt PROPERTIES DEBUG_POSTFIX d)
set_target_properties(bcmopt PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../../../bin)