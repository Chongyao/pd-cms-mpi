#include "cms_io.hpp"
#include "cms_partition.hpp"
#include <fmt/core.h>
#include <fstream>
#include <spdlog/spdlog.h>
namespace zcy {
using namespace std;
void write_partition(const std::string& out_path, const cms_partition& part)
{
    // ofstream part_ofs(fmt::format("{}/{}.partition", out_path, part_name));
    ofstream part_ofs(out_path);
    part_ofs << part.order.size() << " " << part.offset.size() - 2 << "\n";
    for (auto o : part.order) {
        part_ofs << o << "\n";
    }
    for (auto o : part.offset) {
        part_ofs << o << "\n";
    }
    part_ofs.close();
}
cms_partition read_partition(const std::string& filename)
{
    std::ifstream ifs(filename, std::ios::in);
    if (!ifs) {
        std::cerr << "Failed to open file for reading: " << filename << '\n';
        throw std::runtime_error("Failed to open partition file");
    }

    size_t nrows, npart;
    ifs >> nrows >> npart;
    if (!ifs) {
        std::cerr << "Failed to read nrows or npart from: " << filename << '\n';
        throw std::runtime_error("Failed to read nrows or npart from partition file");
    }

    vector<size_t> order(nrows), offset(npart + 2);
    for (size_t i = 0; i < nrows; ++i) {
        ifs >> order[i];
    }
    for (size_t i = 0; i < npart + 2; ++i) {
        ifs >> offset[i];
    }

    return cms_partition(order, offset);
}
tuple<Eigen::MatrixXd, Eigen::VectorXd> readSubstructureEV(const std::string& filename, int nrows, int& nconv)
{
    // new implementation: read (nconv, eigenvalues[nconv], eigenvectors[nrows*nconv]) and return
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        std::cerr << "Failed to open file for reading: " << filename << '\n';
        nconv = 0;
        return { Eigen::MatrixXd(), Eigen::VectorXd() };
    }
    std::cout << "open file: " << filename << "\n";

    ifs.read(reinterpret_cast<char*>(&nconv), sizeof(int));
    if (!ifs) {
        std::cerr << "Failed to read nconv from: " << filename << '\n';
        nconv = 0;
        return { Eigen::MatrixXd(), Eigen::VectorXd() };
    }
    if (nconv < 0) {
        std::cerr << "Invalid nconv in file: " << filename << '\n';
        nconv = 0;
        return { Eigen::MatrixXd(), Eigen::VectorXd() };
    }

    Eigen::VectorXd eigenvalues = Eigen::VectorXd::Zero(nconv);
    if (nconv > 0) {
        ifs.read(reinterpret_cast<char*>(eigenvalues.data()), sizeof(double) * static_cast<size_t>(nconv));
        if (!ifs) {
            std::cerr << "Failed to read eigenvalues from: " << filename << '\n';
            nconv = 0;
            return { Eigen::MatrixXd(), Eigen::VectorXd() };
        }
    }

    const size_t count = static_cast<size_t>(nrows) * static_cast<size_t>(nconv);
    Eigen::MatrixXd eigenvectors;
    if (count > 0) {
        eigenvectors.resize(nrows, nconv);
        ifs.read(reinterpret_cast<char*>(eigenvectors.data()), sizeof(double) * count);
        if (!ifs) {
            std::cerr << "Failed to read eigenvectors from: " << filename << '\n';
            nconv = 0;
            return { Eigen::MatrixXd(), Eigen::VectorXd() };
        }
    } else {
        eigenvectors.resize(nrows, 0);
    }

    ifs.close();
    return { eigenvectors, eigenvalues };
}
void writeSubstructureEV(const std::string& filename, const Eigen::MatrixXd& eigenvectors, const Eigen::VectorXd& eigenvalues)
{
    ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        std::cerr << "Failed to open file for writing: " << filename << "\n";
        return;
    }
    int nconv = static_cast<int>(eigenvalues.size());
    ofs.write(reinterpret_cast<const char*>(&nconv), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(eigenvalues.data()), nconv * sizeof(double));
    ofs.write(reinterpret_cast<const char*>(eigenvectors.data()), eigenvectors.size() * sizeof(double));
    ofs.close();
}

std::string getPartitionName(const std::string& path)
{
    const std::string suf = ".partition";
    size_t posSlash = path.find_last_of('/');
    size_t start = (posSlash == std::string::npos) ? 0 : posSlash + 1;

    size_t end = path.rfind(suf);
    if (end == std::string::npos || end < start) {
        end = path.size();
    }
    return path.substr(start, end - start);
}
std::vector<cms_partition> read_cms_partitions(const std::vector<std::string>& part_paths, size_t num_partition)
{
    if (part_paths.size() != num_partition) {
        std::cerr << "Error: --num-partitions=" << num_partition
                  << " but got " << part_paths.size() << " --part-path values.\n";
        throw std::runtime_error("partition size does not match");
    }

    vector<cms_partition> partitions;
    partitions.reserve(num_partition);
    std::vector<std::vector<size_t>> orders(num_partition), offsets(num_partition);

    for (size_t i = 0; i < num_partition; ++i) {
        partitions.push_back(read_partition(part_paths[i]));
    }

    return partitions;
}


}
