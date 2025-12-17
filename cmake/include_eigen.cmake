CPMAddPackage(
    NAME Eigen3
    GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
    GIT_TAG "3.4.0"
    DOWNLOAD_ONLY YES
)

if(Eigen3_ADDED)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE "${Eigen3_SOURCE_DIR}")
    message(STATUS "Eigen3 defined manually from: ${Eigen3_SOURCE_DIR}")
endif()


add_library(Eigen3_WITH_BLAS INTERFACE)
target_link_libraries(Eigen3_WITH_BLAS
    INTERFACE Eigen3::Eigen
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)
target_compile_definitions(Eigen3_WITH_BLAS
    INTERFACE
    EIGEN_USE_BLAS
    EIGEN_USE_LAPACKE
)

if (${BLA_VENDOR} STREQUAL "Intel10_64lp")
    target_compile_definitions(Eigen3_WITH_BLAS
        INTERFACE
        EIGEN_USE_INTEL_MKL
        EIGEN_MKL_NO_DIRECT_CALL
    )
endif()
