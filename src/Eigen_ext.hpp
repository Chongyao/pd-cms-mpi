#pragma once
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <span>
#include <unordered_map>
namespace Eigen {
template <int dim>
struct VectorCmp {
    using VectorDimd = Eigen::Matrix<double, dim, 1>;
    bool operator()(const VectorDimd& a, const VectorDimd& b) const
    {
        for (size_t i = 0; i < dim; ++i) {
            if (a(i) != b(i))
                return a(i) < b(i);
        }
    }
};
template <int size, typename T_value>
using stdmapVec = std::map<Eigen::Matrix<double, size, 1>, T_value, VectorCmp<size>>;

template <typename Derived>
Eigen::Reshaped<const Derived, Eigen::Dynamic, 1> reshapeAsVector(const Eigen::MatrixBase<Derived>& m)
{
    return Eigen::Reshaped<const Derived, Eigen::Dynamic, 1>(m.derived(), m.size(), 1);
}

template <typename Derived>
Eigen::Map<Eigen::VectorXd> mapAsVector(Eigen::MatrixBase<Derived>& m)
{
    return Eigen::Map<Eigen::VectorXd>(m.derived().data(), m.size());
}

template <typename Derived>
Eigen::Map<const Eigen::VectorXd> mapAsConstVector(const Eigen::MatrixBase<Derived>& m)
{
    return Eigen::Map<const Eigen::VectorXd>(m.derived().data(), m.size());
}

template <typename T>
Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> mapAsConstVector(const std::span<T>& m)
{
    return Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>(m.data(), m.size());
}
template <typename T>
Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> mapAsConstVector(const std::vector<T>& m)
{
    return Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>(m.data(), m.size());
}

template <typename T>
Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> mapAsVector(std::vector<T>& m)
{
    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(m.data(), m.size());
}

/**
 * @brief Extracts a sparse sub-block from a larger sparse matrix.
 *
 * @param K_global The original large sparse matrix.
 * @param global_rows The global row indices to extract.
 * @param global_cols The global column indices to extract.
 * @return A new sparse matrix containing only the specified elements,
 * with indices re-numbered from 0.
 */
template <int Major, int outMajor = Eigen::ColMajor>
Eigen::SparseMatrix<double, Major> extractSparseBlock(
    const Eigen::SparseMatrix<double, Major>& K_global,
    std::span<const size_t> global_rows,
    std::span<const size_t> global_cols)
{
    // 1. Create maps from global indices to the new local indices (0, 1, 2...).
    // This allows for O(1) average lookup time.
    std::unordered_map<size_t, size_t> row_map;
    for (size_t i = 0; i < global_rows.size(); ++i) {
        row_map[global_rows[i]] = i;
    }

    std::unordered_map<size_t, size_t> col_map;
    for (size_t i = 0; i < global_cols.size(); ++i) {
        col_map[global_cols[i]] = i;
    }

    // 2. Iterate through the non-zero elements of the global matrix
    //    and build a triplet list for the sub-block.
    std::vector<Eigen::Triplet<double>> triplets;
    for (int k = 0; k < K_global.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K_global, k); it; ++it) {
            auto row_it = row_map.find(it.row());
            if (row_it != row_map.end()) { // If the row is in our selection
                auto col_it = col_map.find(it.col());
                if (col_it != col_map.end()) { // If the column is also in our selection
                    // Add the element to our triplet list with the new LOCAL indices.
                    triplets.emplace_back(row_it->second, col_it->second, it.value());
                }
            }
        }
    }

    // 3. Create the new sparse matrix from the triplets.
    Eigen::SparseMatrix<double, outMajor> K_sub_block(global_rows.size(), global_cols.size());
    K_sub_block.setFromTriplets(triplets.begin(), triplets.end());
    K_sub_block.makeCompressed();
    return K_sub_block;
}
template <int Major, int dim = 1>
std::vector<Eigen::Triplet<double>> indexingSparseMatrixTriplets(const Eigen::SparseMatrix<double, Major>& A, const std::span<const size_t> v)
{
    // Map original indices to submatrix positions
    std::unordered_map<int, int> index_map;
    for (size_t k = 0; k < v.size(); ++k) {
        for (size_t i = 0; i < dim; ++i) {
            index_map[(v[k] * dim) + i] = k * dim + i; // v(k) is the original index, k is the new index
        }
    }

    // Collect valid elements using triplets
    std::vector<Triplet<double>> triplets;
    triplets.reserve(A.nonZeros());
    for (int col = 0; col < A.outerSize(); ++col) {
        for (SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            int row = it.row();
            auto row_iter = index_map.find(row);
            auto col_iter = index_map.find(col);

            // Check if both row and column are in v
            if (row_iter != index_map.end() && col_iter != index_map.end()) {
                triplets.emplace_back(row_iter->second, col_iter->second, it.value());
            }
        }
    }
    return triplets;
}
template <int Major, int dim = 1>
Eigen::SparseMatrix<double, Major> indexingSparseMatrix(const Eigen::SparseMatrix<double, Major>& A, const std::span<const size_t> v)
{

    // Get triplets for the submatrix
    auto triplets = indexingSparseMatrixTriplets<Major, dim>(A, v);
    // Build the submatrix
    SparseMatrix<double> submatrix(v.size() * dim, v.size() * dim);
    submatrix.setFromTriplets(triplets.begin(), triplets.end());
    submatrix.makeCompressed();
    return submatrix;
}
template <int Major, int dim = 1>
Eigen::SparseMatrix<double, Major> indexingSparseMatrix(const Eigen::SparseMatrix<double, Major>& A, const std::vector<size_t>& v)
{
    return indexingSparseMatrix<Major, dim>(A, std::span<const size_t> { v.data(), v.size() });
}
template <int Major, int dim = 1>
std::vector<Eigen::Triplet<double>> indexingSparseMatrixTriplets(const Eigen::SparseMatrix<double, Major>& A, const std::vector<size_t>& v)
{
    return indexingSparseMatrixTriplets<Major, dim>(A, std::span<const size_t> { v.data(), v.size() });
}

template <typename T = Eigen::SparseMatrix<double>>
T combine_topLeft_and_bottomRight(const T& topLeft, const T& bottomRight)
{
    using Scalar = typename T::Scalar;
    using StorageIndex = typename T::StorageIndex;
    std::vector<Scalar> vals(topLeft.nonZeros() + bottomRight.nonZeros());
    std::copy(topLeft.valuePtr(), topLeft.valuePtr() + topLeft.nonZeros(), vals.data());
    std::copy(bottomRight.valuePtr(), bottomRight.valuePtr() + bottomRight.nonZeros(), vals.data() + topLeft.nonZeros());

    std::vector<StorageIndex> idx(topLeft.nonZeros() + bottomRight.nonZeros());
    std::copy(topLeft.innerIndexPtr(), topLeft.innerIndexPtr() + topLeft.nonZeros(), idx.data());
    std::copy(bottomRight.innerIndexPtr(), bottomRight.innerIndexPtr() + bottomRight.nonZeros(), idx.data() + topLeft.nonZeros());
    Eigen::Map<Eigen::VectorX<StorageIndex>>(idx.data() + topLeft.nonZeros(), bottomRight.nonZeros()).array() += topLeft.innerSize();

    std::vector<StorageIndex> ptr(topLeft.outerSize() + bottomRight.outerSize() + 1);
    std::copy(topLeft.outerIndexPtr(), topLeft.outerIndexPtr() + topLeft.outerSize(), ptr.data());
    std::copy(bottomRight.outerIndexPtr(), bottomRight.outerIndexPtr() + bottomRight.outerSize() + 1, ptr.data() + topLeft.outerSize());
    Eigen::Map<Eigen::VectorX<StorageIndex>>(ptr.data() + topLeft.outerSize(), bottomRight.outerSize() + 1).array() += topLeft.outerIndexPtr()[topLeft.outerSize()];

    Eigen::Map<T> result(topLeft.rows() + bottomRight.rows(), topLeft.cols() + bottomRight.cols(), vals.size(), ptr.data(), idx.data(), vals.data());

    return result;
}
}
