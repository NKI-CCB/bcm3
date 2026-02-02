add_library(likelihoods)

target_sources(likelihoods
  PRIVATE
    LikelihoodDLL.h
    LikelihoodDLL.cpp
    LikelihoodDummy.h
    LikelihoodDummy.cpp
    LikelihoodFactory.h
    LikelihoodFactory.cpp
    LikelihoodIncucytePopulation.h
    LikelihoodIncucytePopulation.cpp
    LikelihoodMitosisTimeEstimation.h
    LikelihoodMitosisTimeEstimation.cpp
    LikelihoodODE.h
    LikelihoodODE.cpp
    LikelihoodPharmacokineticTrajectory.h
    LikelihoodPharmacokineticTrajectory.cpp
    LikelihoodPopPKTrajectory.h
    LikelihoodPopPKTrajectory.cpp
	TestLikelihoodBanana.h
	TestLikelihoodBanana.cpp
	TestLikelihoodCircular.h
	TestLikelihoodCircular.cpp
    TestLikelihoodMultimodalGaussians.h
    TestLikelihoodMultimodalGaussians.cpp
    TestLikelihoodTruncatedT.h
    TestLikelihoodTruncatedT.cpp
)

target_link_libraries(likelihoods
  PUBLIC
    utils
    sampler
    odecommon
    pharmaco
)

if(INCLUDE_CELLPOP)
    target_compile_definitions(likelihoods PRIVATE BCM_INCLUDE_CELLPOP=1)
    target_link_libraries(likelihoods PUBLIC cellpop)
endif(INCLUDE_CELLPOP)

if(INCLUDE_fISA)
    target_compile_definitions(likelihoods PRIVATE BCM_INCLUDE_fISA=1)
    target_link_libraries(likelihoods PUBLIC fISA)
endif(INCLUDE_fISA)

set_target_properties(likelihoods PROPERTIES DEBUG_POSTFIX d)
target_precompile_headers(likelihoods PRIVATE ../utils/Utils.h)
