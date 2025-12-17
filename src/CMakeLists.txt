aux_source_directory(. source)
file(GLOB header *.h)

add_library(cms SHARED ${source})
target_link_libraries(cms
  PUBLIC 
  fmt::fmt
  spdlog::spdlog
  Spectra
  SuiteSparse::CHOLMOD
  Eigen3_WITH_BLAS
)