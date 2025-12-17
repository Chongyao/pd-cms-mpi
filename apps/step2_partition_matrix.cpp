#include <string>
#include <spdlog/spdlog.h>
#include "Eigen_ext.hpp"
#include "save_spm.hpp"
#include "cms_partition.hpp"
#include "cms_io.hpp"
using namespace zcy;
int main(int argc, char** argv){
    if(argc != 5){
        spdlog::info("Usage: {} K_path M_path partition_file_path", argv[0]);
        return 1;
    }
    const std::string K_path = argv[1];
    const std::string M_path = argv[2];
    const std::string partition_file_path = argv[3];
    const std::string output_path = argv[4];

    Eigen::SparseMatrix<double> K, M;
    zcy::io::read_spm(K_path.c_str(), K);
    zcy::io::read_spm(M_path.c_str(), M);
    spdlog::info("Read K: {} x {}, nnz = {}", K.rows(), K.cols(), K.nonZeros());
    spdlog::info("Read M: {} x {}, nnz = {}", M.rows(), M.cols(), M.nonZeros());

    std::vector<std::string> partitions;
    {
      // split partitio_file_path by ','
      size_t start = 0;
      size_t end = partition_file_path.find(',');
      while (end != std::string::npos) {
          partitions.push_back(partition_file_path.substr(start, end - start));
          start = end + 1;
          end = partition_file_path.find(',', start);      
      }
      partitions.push_back(partition_file_path.substr(start));
    }
    size_t rank_id = 0;
    for(auto &part_file : partitions){
      spdlog::info("Partition file: {}", part_file);
      cms_partition p = read_partition(part_file);
      const auto &order = p.order;
      const auto &offset = p.offset;
      spdlog::info("Read partition: num_block = {}", p.num_block());

      for(size_t sub_id = 0; sub_id < p.num_block(); ++sub_id){
        std::span<const size_t> substructure(order.data() + offset[sub_id],
                                           offset[sub_id + 1] - offset[sub_id]);
        Eigen::SparseMatrix<double> Kii = Eigen::indexingSparseMatrix(K, substructure);
        Eigen::SparseMatrix<double> Mii = Eigen::indexingSparseMatrix(M, substructure);

        std::string Kii_path = fmt::format("{}/K_block_{}.csc", output_path, rank_id);
        std::string Mii_path = fmt::format("{}/M_block_{}.csc", output_path, rank_id);
        zcy::io::write_spm(Kii_path.c_str(), Kii);
        zcy::io::write_spm(Mii_path.c_str(), Mii);
        ++rank_id;
        spdlog::info("Write Kii: {} x {}, nnz = {} to {}", Kii.rows(), Kii.cols(), Kii.nonZeros(), Kii_path);
        spdlog::info("Write Mii: {} x {}, nnz = {} to {}", Mii.rows(), Mii.cols(), Mii.nonZeros(), Mii_path);
      }
    }
    



    return 0;
}
