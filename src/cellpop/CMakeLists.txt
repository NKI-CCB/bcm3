add_library(cellpop)

target_sources(cellpop
  PRIVATE
    Cell.h
    Cell.cpp
    CellPopulationLikelihood.h
    CellPopulationLikelihood.cpp
    DataLikelihoodBase.h
    DataLikelihoodBase.cpp
    DataLikelihoodDuration.h
    DataLikelihoodDuration.cpp
    DataLikelihoodTimeCourse.h
    DataLikelihoodTimeCourse.cpp
    DataLikelihoodTimePoints.h
    DataLikelihoodTimePoints.cpp
    Experiment.h
    Experiment.cpp
    TreatmentTrajectory.h
    TreatmentTrajectory.cpp
    VariabilityDescription.h
    VariabilityDescription.cpp
)

target_include_directories(cellpop
  PUBLIC
    .
)

target_link_libraries(cellpop
  PUBLIC
    utils
    sbml
    ${CMAKE_DL_LIBS}
    hungarian_matching
  PRIVATE
    odecommon
    sampler
)

target_precompile_headers(cellpop PRIVATE ../utils/Utils.h)
set_target_properties(cellpop PROPERTIES DEBUG_POSTFIX d)
