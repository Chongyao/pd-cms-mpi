add_executable(main main.cpp)
target_link_libraries(main
  PRIVATE 
  fmt::fmt
  spdlog::spdlog
  Spectra
  SuiteSparse::CHOLMOD
  Eigen3_WITH_BLAS
)
