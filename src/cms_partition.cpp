#include "cms_partition.hpp"
namespace zcy {
cms_partition::cms_partition(const std::tuple<std::vector<size_t>, std::vector<size_t>>& order_offset)
    : order(std::get<0>(order_offset))
    , offset(std::get<1>(order_offset))
{
}
cms_partition::cms_partition(const std::vector<size_t>& order, const std::vector<size_t>& offset)
    : order(order)
    , offset(offset)
{
}

size_t cms_partition::num_block() const
{
    return offset.size() - 2;
}

cms_res::cms_res(const std::vector<Eigen::MatrixXd>& phi, const std::vector<Eigen::VectorXd>& evalue)
    : phi(phi)
    , evalue(evalue)
{
}
double cms_res::memory_cost_in_mb() const
{

    double cost = 0.0;
    for (const auto& phi : phi) {
        cost += phi.size() * sizeof(double) / (1024.0 * 1024.0);
    }
    for (const auto& evalue : evalue) {
        cost += evalue.size() * sizeof(double) / (1024.0 * 1024.0);
    }
    for (const auto& M_blk : M_red_blk) {
        cost += M_blk.size() * sizeof(double) / (1024.0 * 1024.0);
    }
    cost += K_red_blk.size() * sizeof(double) / (1024.0 * 1024.0);
    return cost;
}
double cms_res::memory_cost_of_bb_in_mb() const
{
    double cost = 0.0;
    for (const auto& M_blk : M_red_blk) {
        cost += M_blk.size() * sizeof(double) / (1024.0 * 1024.0);
    }
    cost += K_red_blk.size() * sizeof(double) / (1024.0 * 1024.0);
    return cost;
}
double cms_res::memory_cost_of_eig_substructure() const
{
    double cost = 0;
    for (const auto& phi : phi) {
        cost += phi.size() * sizeof(double) / (1024.0 * 1024.0);
    }
    for (const auto& evalue : evalue) {
        cost += evalue.size() * sizeof(double) / (1024.0 * 1024.0);
    }
    return cost;
}
}
