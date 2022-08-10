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


target_link_libraries(sbml
  PUBLIC
    Eigen
    Boost::headers
    debug ${SBML_LIBRARIES_DEBUG}
    optimized ${SBML_LIBRARIES}
)

set_target_properties(sbml PROPERTIES DEBUG_POSTFIX d)
target_include_directories(sbml PRIVATE ../utils ../sampler ../odecommon ../../dependencies/cvode-5.3.0/include)
target_include_directories(sbml SYSTEM PUBLIC ${Boost_INCLUDE_DIR} ${SBML_INCLUDE_DIR})
target_precompile_headers(sbml PRIVATE ../utils/Utils.h)
