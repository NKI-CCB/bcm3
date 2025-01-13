add_library(pharmaco)

target_sources(pharmaco
  PRIVATE
    PharmacokineticModel.h
    PharmacokineticModel.cpp
    PharmacoLikelihoodSingle.h
    PharmacoLikelihoodSingle.cpp
)

target_include_directories(pharmaco
  PUBLIC
    .
)

target_link_libraries(pharmaco
  PUBLIC
    utils
    ${CMAKE_DL_LIBS}
  PRIVATE
    sampler
)

target_precompile_headers(pharmaco PRIVATE ../utils/Utils.h)
set_target_properties(pharmaco PROPERTIES DEBUG_POSTFIX d)
