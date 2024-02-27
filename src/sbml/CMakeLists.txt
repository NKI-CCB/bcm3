add_library(sbml)

target_sources(sbml
  PRIVATE
    SBMLAssignmentRule.h
    SBMLAssignmentRule.cpp
    SBMLModel.h
    SBMLModel.cpp
    SBMLModelParameters.h
    SBMLModelParameters.cpp
    SBMLRatelaws.h
    SBMLRatelaws.cpp
    SBMLReaction.h
    SBMLReaction.cpp
    SBMLSpecies.h
    SBMLSpecies.cpp
)

target_include_directories(sbml
  PUBLIC
    .
    ../../dependencies/libsbml-5.20.2/
  PRIVATE
)

target_link_libraries(sbml
  PUBLIC
    utils
  PRIVATE
    sampler
    odecommon
    libsbml
)

set_target_properties(sbml PROPERTIES DEBUG_POSTFIX d)
target_precompile_headers(sbml PRIVATE ../utils/Utils.h)
