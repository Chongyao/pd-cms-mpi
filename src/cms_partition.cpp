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



}
