add_library(dynamicISA)

target_sources(dynamicISA
  PRIVATE
    dynamicISALikelihood.h
    dynamicISALikelihood.cpp
    dynamicISAExperiment.h
    dynamicISAExperiment.cpp
    SignalingLink.h
    SignalingLink.cpp
    SignalingModel.h
    SignalingModel.cpp
    SignalingMolecule.h
    SignalingMolecule.cpp
)

target_include_directories(dynamicISA
  PUBLIC
    ../utils
    ../odecommon
    ../sampler
    ../sbml
  SYSTEM PUBLIC
    ${SBML_INCLUDE_DIR}
)

target_link_libraries(dynamicISA
  PUBLIC
    utils
    sampler
    sbml
    odecommon
)

target_precompile_headers(dynamicISA PRIVATE ../utils/Utils.h)
set_target_properties(dynamicISA PROPERTIES DEBUG_POSTFIX d)