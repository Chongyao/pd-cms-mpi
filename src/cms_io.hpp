#pragma once
#include "cms_partition.hpp"
#include <tuple>
namespace zcy {
// io of partition
void write_partition(const std::string& path, const cms_partition& part);
cms_partition read_partition(const std::string& filename);

std::string getPartitionName(const std::string& path);

// io of eigen vectors of substructures
std::tuple<Eigen::MatrixXd, Eigen::VectorXd> readSubstructureEV(const std::string& filename, int nrows, int& nconv);
void writeSubstructureEV(const std::string& filename, const Eigen::MatrixXd& eigenvectors, const Eigen::VectorXd& eigenvalues);
cms_res readCmsResults(const std::string& input_path, const std::string& partition_name, const cms_partition& part);
void writeCmsResults(const std::string& output_path, const std::string& partition_name, const cms_partition& part, const cms_res& result);

// io of cms partitions and results together
std::vector<cms_partition> read_cms_partitions(const std::vector<std::string>& part_paths, size_t num_partition, const std::string& output_path);

std::tuple<std::vector<cms_partition>, std::vector<cms_res>> read_cms_partitions_and_res(const std::vector<std::string>& part_paths, size_t num_partition, const std::string& output_path);
void write_cms_partitions_and_res(const std::vector<cms_partition>& parts, const std::vector<cms_res>& ress, const std::string& output_path, const std::vector<std::string>& part_paths);
}
