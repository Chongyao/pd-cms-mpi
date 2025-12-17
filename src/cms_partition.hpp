#pragma once
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cstddef>
#include <span>
#include <iostream>
#include <span>
#include <unordered_set>
#include <vector>

namespace zcy {
struct cms_partition {
    cms_partition() = default;
    cms_partition(const std::tuple<std::vector<size_t>, std::vector<size_t>>& order_offset);
    cms_partition(const std::vector<size_t>& order, const std::vector<size_t>& offset);

    std::vector<size_t> order;
    std::vector<size_t> offset;

    std::vector<std::vector<size_t>> blkwise_bd;

    void set_blockwise_boundary(const Eigen::SparseMatrix<double>& K)
    {
        size_t num_blocks = offset.size() - 2;

        std::vector<int> tag(order.size(), -1);
        for (size_t pi = 0; pi < num_blocks; ++pi) {
            for (size_t i = offset[pi]; i < offset[pi + 1]; ++i) {
                tag[order[i]] = pi;
            }
        }

        // loop through boundary
        blkwise_bd.resize(num_blocks);
        std::span<const size_t> bd_span(order.begin() + offset[num_blocks], order.end());
        for (size_t i = 0; i < bd_span.size(); ++i) {
            size_t bd_dof = bd_span[i];
            for (Eigen::SparseMatrix<double>::InnerIterator it(K, bd_dof); it; ++it) {
                size_t nz_dof = it.index();
                if (tag[nz_dof] != -1) {
                    blkwise_bd[tag[nz_dof]].push_back(bd_dof);
                }
            }
        }

        for (auto& bd : blkwise_bd) {
            std::sort(bd.begin(), bd.end());
        }
    }

    size_t num_block() const;
};

}
