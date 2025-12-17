# CPMAddPackage(
#   NAME Eigen3
#   GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
#   GIT_TAG "3.4.0"
#   PATCH_COMMAND
#   git restore . && git apply "${PROJECT_SOURCE_DIR}/cmake/eigen3_skip_build_demo.patch"
#   OPTIONS
#   "BUILD_TESTING OFF"
#   "EIGEN_BUILD_DOC OFF"
#   "EIGEN_BUILD_PKGCONFIG OFF"

# )

# add_library(Eigen3_WITH_BLAS INTERFACE)
# target_link_libraries(Eigen3_WITH_BLAS
#   INTERFACE Eigen3::Eigen
#   ${BLAS_LIBRARIES}
#   ${LAPACK_LIBRARIES}
# )
# target_compile_definitions(Eigen3_WITH_BLAS
#   INTERFACE
#   EIGEN_USE_MKL_ALL
#   EIGEN_USE_BLAS
#   EIGEN_USE_LAPACKE
#   INTERFACE
#   "UI"
# )

CPMAddPackage(
    NAME Eigen3
    GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
    GIT_TAG "3.4.0"
    DOWNLOAD_ONLY YES  # <--- 关键：只下载，不add_subdirectory
)

# 如果下载成功，手动创建一个干净的 Eigen3::Eigen 目标
if(Eigen3_ADDED)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE "${Eigen3_SOURCE_DIR}")
    
    # 打印一下，确认路径正确
    message(STATUS "Eigen3 defined manually from: ${Eigen3_SOURCE_DIR}")
endif()

# --------------- 你的剩余逻辑保持不变 ---------------

add_library(Eigen3_WITH_BLAS INTERFACE)
target_link_libraries(Eigen3_WITH_BLAS
    INTERFACE Eigen3::Eigen
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)
target_compile_definitions(Eigen3_WITH_BLAS
    INTERFACE
    EIGEN_USE_MKL_ALL
    EIGEN_USE_BLAS
    EIGEN_USE_LAPACKE
    EIGEN_MKL_NO_DIRECT_CALL  # <--- 添加这一行
)
