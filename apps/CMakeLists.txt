add_executable(main main.cpp)
target_link_libraries(main
  cms
)



set(STEPS
  1_env_check
  2_partition_matrix
  3_mpi_loader
)

foreach(STEP IN LISTS STEPS)
    message(${STEP})
    add_executable(step${STEP} step${STEP}.cpp)
    target_link_libraries(step${STEP} PRIVATE cms MPI::MPI_CXX ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endforeach()

