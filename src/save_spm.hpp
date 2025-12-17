#pragma once
#include <Eigen/SparseCore>
#include <cstdint>
#include <fstream>
#include <iostream>
namespace zcy {
namespace io {
    template <int Major>
    int write_spm(const char* path, const Eigen::SparseMatrix<double, Major>& A)
    {
        std::ofstream ofs(path, std::ios::binary);
        if (ofs.fail()) {
            return __LINE__;
        }
        const int64_t mat_size[4] = { A.rows(), A.cols(), A.nonZeros(), Major };
        ofs.write((const char*)&mat_size[0], 4 * sizeof(int64_t));
        // write value, Innerindex, OuterIndex
        ofs.write((const char*)A.valuePtr(), A.nonZeros() * sizeof(double));
        ofs.write((const char*)A.innerIndexPtr(), A.nonZeros() * sizeof(int));
        if (Major == Eigen::ColMajor) {
            ofs.write((const char*)A.outerIndexPtr(), (A.cols() + 1) * sizeof(int));
        } else {
            ofs.write((const char*)A.outerIndexPtr(), (A.rows() + 1) * sizeof(int));
        }
        ofs.close();
        return 0;
    }

    template <int Major>
    int write_large_spm(const char* path, const Eigen::SparseMatrix<double, Major, std::ptrdiff_t>& A)
    {
        std::ofstream ofs(path, std::ios::binary);
        if (ofs.fail()) {
            return __LINE__;
        }
        const std::ptrdiff_t mat_size[4] = { A.rows(), A.cols(), A.nonZeros(), Major };
        ofs.write((const char*)&mat_size[0], 4 * sizeof(std::ptrdiff_t));
        // write value, Innerindex, OuterIndex
        ofs.write((const char*)A.valuePtr(), A.nonZeros() * sizeof(double));
        ofs.write((const char*)A.innerIndexPtr(), A.nonZeros() * sizeof(std::ptrdiff_t));
        if (Major == Eigen::ColMajor) {
            ofs.write((const char*)A.outerIndexPtr(), (A.cols() + 1) * sizeof(std::ptrdiff_t));
        } else {
            ofs.write((const char*)A.outerIndexPtr(), (A.rows() + 1) * sizeof(std::ptrdiff_t));
        }
        ofs.close();
        return 0;
    }


    template <int Major>
    int _read_spm(const char* path, Eigen::SparseMatrix<double, Major>& A)
    {
        std::ifstream ifs(path, std::ios::binary);
        if (ifs.fail()) {
            std::cerr << "Open " << path << "failed.\n";
            return 1;
        }
        int64_t mat_size[4];
        ifs.read((char*)&mat_size[0], 4 * sizeof(int64_t));
        A.resize(mat_size[0], mat_size[1]);
        Eigen::VectorXd val;
        Eigen::VectorXi idx, offset;
        val.resize(mat_size[2]);
        idx.resize(mat_size[2]);
        if (Major != mat_size[3]) {
            return 2;
        }
        if (Major == Eigen::ColMajor) {
            offset.resize(mat_size[1] + 1);
        } else {
            offset.resize(mat_size[0] + 1);
        }
        ifs.read((char*)val.data(), val.size() * sizeof(double));
        ifs.read((char*)idx.data(), idx.size() * sizeof(int));
        ifs.read((char*)offset.data(), offset.size() * sizeof(int));
        A = Eigen::Map<Eigen::SparseMatrix<double, Major>>(mat_size[0], mat_size[1], mat_size[2], offset.data(), idx.data(), val.data());

        ifs.close();
        return 0;
    }
    template <int Major>
    constexpr int EigenTheOtherMajor()
    {
        if constexpr (Major == Eigen::RowMajor) {
            return Eigen::ColMajor;
        } else {
            return Eigen::RowMajor;
        }
    }
    template <int Major>
    int read_spm(const char* path, Eigen::SparseMatrix<double, Major>& A)
    {
        int rtn = _read_spm<Major>(path, A);
        if (rtn == 0 || rtn == 1) {
            return rtn;
        }
        if (rtn == 2) {
            Eigen::SparseMatrix<double, EigenTheOtherMajor<Major>()> _A;
            rtn = _read_spm<EigenTheOtherMajor<Major>()>(path, _A);
            A = _A;
        }
        return rtn;
    }

}
}
