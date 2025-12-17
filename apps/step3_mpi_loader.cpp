// apps/step3_mpi_loader.cpp
#include <mpi.h>
#include <spdlog/spdlog.h>
#include <fmt/core.h>
#include "save_spm.hpp" // 包含 zcy::io::read_spm

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 默认数据目录为当前目录，也可以通过参数传入
    std::string data_dir = ".";
    if (argc > 1) {
        data_dir = argv[1];
    }

    try {
        // 1. 构造文件名 (严格对应 step2 的输出)
        // 注意：这里假设 step2 生成的是二进制格式，尽管后缀是 .csc
        std::string K_path = fmt::format("{}/K_block_{}.csc", data_dir, rank);
        std::string M_path = fmt::format("{}/M_block_{}.csc", data_dir, rank);

        // 2. 定义本地矩阵
        Eigen::SparseMatrix<double> K_local, M_local;

        // 3. 读取数据
        if (zcy::io::read_spm(K_path.c_str(), K_local) != 0) {
            throw std::runtime_error(fmt::format("Failed to read K from {}", K_path));
        }
        if (zcy::io::read_spm(M_path.c_str(), M_local) != 0) {
            throw std::runtime_error(fmt::format("Failed to read M from {}", M_path));
        }

        // 4. 验证与打印 (为了不刷屏，让 Rank 0 打印汇总，或者有序打印)
        // 这里简单地每个 Rank 打印一句日志
        spdlog::info("Rank {:02d}: Loaded K({}, {}) nnz={} | M({}, {}) nnz={}", 
                     rank, 
                     K_local.rows(), K_local.cols(), K_local.nonZeros(),
                     M_local.rows(), M_local.cols(), M_local.nonZeros());

        // 可以在这里加一个简单的数值检查，比如检查对角线元素不为0 (可选)
        
    } catch (const std::exception& e) {
        spdlog::error("Rank {:02d} Error: {}", rank, e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Finalize();
    return 0;
}