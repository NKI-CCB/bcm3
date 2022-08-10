add_library(fISA)

target_sources(fISA
  PRIVATE
    fISAExperiment.h
    fISAExperiment.cpp
    fISAExperimentDrugRange.h
    fISAExperimentDrugRange.cpp
    fISAExperimentIncucyteSequential.h
    fISAExperimentIncucyteSequential.cpp
    fISAExperimentSingleCondition.h
    fISAExperimentSingleCondition.cpp
    fISALikelihood.h
    fISALikelihood.cpp
    SignalingNetwork.h
    SignalingNetwork.cpp
)

target_include_directories(fISA
  PUBLIC
    .
)

target_link_libraries(fISA
  PRIVATE
    sampler
    sbml
)

target_precompile_headers(fISA PRIVATE ../utils/Utils.h)
set_target_properties(fISA PROPERTIES DEBUG_POSTFIX d)