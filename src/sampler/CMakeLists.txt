add_library(sampler)

target_sources(sampler
  PRIVATE
    DataReference.h
    DataReference.cpp
    GMM.h
    GMM.cpp
    Likelihood.h
    Likelihood.cpp
    Prior.h
    Prior.cpp
    PriorIndependence.h
    PriorIndependence.cpp
    ProgressIndicator.h
    ProgressIndicator.cpp
    ProgressIndicatorConsole.h
    ProgressIndicatorConsole.cpp
    SampleHandlerNetCDF.h
    SampleHandlerNetCDF.cpp
    SampleHandlerStoreMaxAPosteriori.h
    SampleHandlerStoreMaxAPosteriori.cpp
    SampleHandlerTSV.h
    SampleHandlerTSV.cpp
    SamplerPT.h
    SamplerPT.cpp
    SamplerPTChain.h
    SamplerPTChain.cpp
    UnivariateMarginal.h
    UnivariateMarginal.cpp
    VariableSet.h
    VariableSet.cpp
)

target_include_directories(sampler
  PUBLIC
    ../stats
    ../utils
)

target_link_libraries(sampler
  PUBLIC
    stats
    utils
)

set_target_properties(sampler PROPERTIES DEBUG_POSTFIX d)
target_precompile_headers(sampler PRIVATE ../utils/Utils.h)