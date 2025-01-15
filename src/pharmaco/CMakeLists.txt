add_library(pharmaco)

target_sources(pharmaco
  PRIVATE
    DrugConstants.h
    PharmacokineticModel.h
    PharmacokineticModel.cpp
    PharmacoLikelihoodPopulation.h
    PharmacoLikelihoodPopulation.cpp
    PharmacoLikelihoodSingle.h
    PharmacoLikelihoodSingle.cpp
    PharmacoPatient.h
    PharmacoPatient.cpp
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
