add_library(bcmrbridge SHARED)

target_sources(bcmrbridge
  PRIVATE
    interface.h
    interface.cpp
    interface_cellpop.cpp
    interface_fISA.cpp
    interface_incucyte.cpp
    interface_ODE.cpp
    interface_PK.cpp
    interface_popPK.cpp
    definitions.def
)

target_include_directories(bcmrbridge
  PRIVATE
    ../likelihoods
)

target_link_libraries(bcmrbridge
  PRIVATE
    likelihoods
)

set_target_properties(bcmrbridge PROPERTIES DEBUG_POSTFIX d)
set_target_properties(bcmrbridge PROPERTIES PREFIX "")
set_target_properties(bcmrbridge PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../../../bin)
set_target_properties(bcmrbridge PROPERTIES LIBRARY_OUTPUT_DIRECTORY ../../../bin)
