add_executable(main main.cpp)
target_link_libraries(main
  cms
)


# === Step 1 测试程序 ===
add_executable(step1_check step1_env_check.cpp)

target_link_libraries(step1_check PRIVATE
  MPI::MPI_CXX
  cms
  ${BLAS_LIBRARIES}      # <--- 新增：链接 BLAS (MKL)
  ${LAPACK_LIBRARIES}    # <--- 新增：通常最好也加上 LAPACK
)

# 如果开启了 OpenMP，链接它
if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(step1_check PRIVATE OpenMP::OpenMP_CXX)
endif()