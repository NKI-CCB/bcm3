add_executable(bcminf)

target_sources(bcminf
  PRIVATE
    main.cpp
)

target_include_directories(bcminf
  PRIVATE
    ../likelihoods
)

target_link_libraries(bcminf
  PRIVATE
    likelihoods
)

set_target_properties(bcminf PROPERTIES DEBUG_POSTFIX d)
set_target_properties(bcminf PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../../../bin)