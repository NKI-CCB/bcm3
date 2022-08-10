add_library(likelihoods)

target_sources(likelihoods
  PRIVATE
    LikelihoodFactory.h
    LikelihoodFactory.cpp
    LikelihoodPharmacokineticTrajectory.h
    LikelihoodPharmacokineticTrajectory.cpp
    LikelihoodPopPKTrajectory.h
    LikelihoodPopPKTrajectory.cpp
	TestLikelihoodBanana.h
	TestLikelihoodBanana.cpp
	TestLikelihoodCircular.h
	TestLikelihoodCircular.cpp
)

target_include_directories(likelihoods
  PUBLIC
    ../dynamicISA
    ../fISA
    ../cellpop
)

target_link_libraries(likelihoods
  PUBLIC
    dynamicISA
    fISA
    cellpop
    utils
)

set_target_properties(likelihoods PROPERTIES DEBUG_POSTFIX d)