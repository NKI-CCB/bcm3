add_library(bcmrbridge SHARED)

target_sources(bcmrbridge
  PRIVATE
    interface.h
    interface.cpp
    interface_incucyte.cpp
    interface_ODE.cpp
    interface_pharmaco_population.cpp
    interface_pharmaco_single.cpp
    interface_PK.cpp
    interface_popPK.cpp
    definitions.def
)

if(INCLUDE_CELLPOP)
    target_sources(bcmrbridge PRIVATE interface_cellpop.cpp definitions_cellpop.def)
endif(INCLUDE_CELLPOP)

if(INCLUDE_fISA)
    target_sources(bcmrbridge PRIVATE interface_fISA.cpp definitions_fISA.def)
endif(INCLUDE_fISA)

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
