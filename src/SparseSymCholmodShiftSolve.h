#ifndef SPECTRA_SPARSE_SYM_SHIFT_SOLVE_H
#define SPECTRA_SPARSE_SYM_SHIFT_SOLVE_H

#include "timer.hpp"
#include <Eigen/CholmodSupport>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <iostream>
#include <stdexcept>

namespace Spectra {

///
/// \ingroup MatOp
///
/// This class defines the shift-solve operation on a sparse real symmetric matrix \f$A\f$,
/// i.e., calculating \f$y=(A-\sigma I)^{-1}x\f$ for any real \f$\sigma\f$ and
/// vector \f$x\f$. It is mainly used in the SymEigsShiftSolver eigen solver.
///
/// \tparam Scalar_      The element type of the matrix, for example,
///                      `float`, `double`, and `long double`.
/// \tparam Uplo         Either `Eigen::Lower` or `Eigen::Upper`, indicating which
///                      triangular part of the matrix is used.
/// \tparam Flags        Either `Eigen::ColMajor` or `Eigen::RowMajor`, indicating
///                      the storage format of the input matrix.
/// \tparam StorageIndex The type of the indices for the sparse matrix.
///
template <typename Scalar_, int Uplo = Eigen::Lower, int Flags = Eigen::ColMajor, typename StorageIndex = int>
class SparseSymCholmodShiftSolve {
public:
    ///
    /// Element type of the matrix.
    ///
    using Scalar = Scalar_;

private:
    using Index = Eigen::Index;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MapConstVec = Eigen::Map<const Vector>;
    using MapVec = Eigen::Map<Vector>;
    using SparseMatrix = Eigen::SparseMatrix<Scalar, Flags, StorageIndex>;
    using ConstGenericSparseMatrix = const Eigen::Ref<const SparseMatrix>;

    ConstGenericSparseMatrix m_mat;
    const Index m_n;
    Eigen::CholmodSupernodalLLT<SparseMatrix> m_solver;
    mutable size_t cnt { 0 };
    // Eigen::CholmodSimplicialLDLT<SparseMatrix> m_solver;
    // Eigen::SparseLU<SparseMatrix> m_solver;
    // Eigen::SimplicialLDLT<SparseMatrix> m_solver;

public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat An **Eigen** sparse matrix object, whose type can be
    /// `Eigen::SparseMatrix<Scalar, ...>` or its mapped version
    /// `Eigen::Map<Eigen::SparseMatrix<Scalar, ...> >`.
    ///
    template <typename Derived>
    SparseSymCholmodShiftSolve(const Eigen::SparseMatrixBase<Derived>& mat)
        : m_mat(mat)
        , m_n(mat.rows())
    {
        static_assert(
            static_cast<int>(Derived::PlainObject::IsRowMajor) == static_cast<int>(SparseMatrix::IsRowMajor),
            "SparseSymShiftSolve: the \"Flags\" template parameter does not match the input matrix (Eigen::ColMajor/Eigen::RowMajor)");

        if (mat.rows() != mat.cols()) {
            throw std::invalid_argument("SparseSymShiftSolve: matrix must be square");
        }

        // {
        //     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig;
        //     eig.compute(Eigen::MatrixXd(mat));
        //     std::cout << eig.eigenvalues().topRows(10) << std::endl;
        // }
        zcy::TIMING::TIC("analyze");
        m_solver.analyzePattern(mat);
        zcy::TIMING::TOC("analyze", false);
        if (m_solver.info() != Eigen::Success) {
            throw std::invalid_argument("SparseSymShiftSolve: analyze pattern fail");
        }
    }

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_n; }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_n; }

    ///
    /// Set the real shift \f$\sigma\f$.
    ///
    void set_shift(const Scalar& sigma)
    {
        SparseMatrix mat = m_mat.template selfadjointView<Uplo>();

        SparseMatrix identity(m_n, m_n);
        identity.setIdentity();
        mat = mat - sigma * identity;

        // mat.diagonal().array() -= sigma;

        // m_solver.isSymmetric(true);
        zcy::TIMING::TIC("factorize");
        m_solver.factorize(mat);
        zcy::TIMING::TOC("factorize", false);
        if (m_solver.info() != Eigen::Success) {
            throw std::invalid_argument("SparseSymShiftSolve: factorization failed with the given shift");
        }
    }

    ///
    /// Perform the shift-solve operation \f$y=(A-\sigma I)^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(A - sigma * I) * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in, m_n);
        MapVec y(y_out, m_n);
        y.noalias() = m_solver.solve(x);
        ++cnt;
        // std::cout << cnt++ << "\n";
    }
};

} // namespace Spectra

#endif // SPECTRA_SPARSE_SYM_SHIFT_SOLVE_H
