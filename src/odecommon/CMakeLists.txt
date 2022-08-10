add_library(odecommon)

target_sources(odecommon
  PRIVATE
    CVODESolver.h
    CVODESolver.cpp
    CVODESolverDelay.h
    CVODESolverDelay.cpp
    LinearAlgebraSelector.h
    nvector_serial_eigen.h
    nvector_serial_eigen.cpp
    ODESolver.h
    ODESolver.cpp
    ODESolverDP5.h
    ODESolverDP5.cpp
    sunlinsol_dense_eigen.h
    sunlinsol_dense_eigen.cpp
    sunlinsol_sparse_eigen.h
    sunlinsol_sparse_eigen.cpp
    sunmatrix_dense_eigen.h
    sunmatrix_dense_eigen.cpp
    sunmatrix_sparse_eigen.h
    sunmatrix_sparse_eigen.cpp
)

target_include_directories(odecommon
  PUBLIC
    .
     ../../dependencies/cvode-5.3.0/include
)

target_link_libraries(odecommon
  PUBLIC
    utils
    sundials_cvode_static
    sundials_nvecserial_static
)

set_target_properties(odecommon PROPERTIES DEBUG_POSTFIX d)
target_precompile_headers(odecommon PRIVATE ../utils/Utils.h)
