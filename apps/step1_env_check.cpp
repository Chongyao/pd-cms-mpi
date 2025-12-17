// apps/step1_env_check.cpp
#include <iostream>
#include <mpi.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Spectra/SymEigsSolver.h>
#include <spdlog/spdlog.h>

// 简单的宏检查 OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    // 1. 初始化 MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        spdlog::info("Step 1 Environment Check");
        spdlog::info("MPI Initialized. Size: {}", size);

        // 2. 检查 OpenMP
#ifdef _OPENMP
        spdlog::info("OpenMP detected. Max threads: {}", omp_get_max_threads());
#else
        spdlog::warn("OpenMP NOT detected.");
#endif

        // 3. 检查 Eigen
        Eigen::MatrixXd m(2, 2);
        m(0, 0) = 3;
        m(1, 0) = 2.5;
        m(0, 1) = -1;
        m(1, 1) = m(1, 0) + m(0, 1);
        spdlog::info("Eigen check: \n{}", (std::stringstream() << m).str());

        // 4. 检查 Spectra (简单编译测试)
        // 定义一个简单的 10x10 对角矩阵
        int n = 10;
        Eigen::MatrixXd M = Eigen::MatrixXd::Random(n, n);
        Eigen::MatrixXd A = M + M.transpose(); // 对称矩阵

        // Spectra 求解器: 求解 A 的 3 个最大特征值
        Spectra::DenseSymMatProd<double> op(A);
        Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigs(op, 3, 6);
        
        eigs.init();
        int nconv = eigs.compute(Spectra::SortRule::LargestAlge);

        if (eigs.info() == Spectra::CompInfo::Successful) {
            Eigen::VectorXd evalues = eigs.eigenvalues();
            spdlog::info("Spectra check: Found {} eigenvalues.", evalues.size());
            spdlog::info("Top eigenvalue: {}", evalues(0));
        } else {
            spdlog::error("Spectra computation failed.");
        }
    }else{
        spdlog::info("MPI Rank {} ready.", rank);
    }

    MPI_Finalize();
    return 0;
}