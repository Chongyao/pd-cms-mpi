#include "eig.hpp"
#include "SparseSymCholmodShiftSolve.h"
#include "sym_shift_invert_with_cholmod.hpp"
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/Util/GEigsMode.h>
#include <iostream>
#include <stdexcept>

void solve_minimal_eig(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, double shift, int maxits, double tol)
{
    // using T_A_op = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using T_A_op = Spectra::SymShiftInvert<double, Eigen::Dense, Eigen::Dense>;
    // using T_B_op = Spectra::SparseSymMatProd<double>;
    using T_B_op = Spectra::DenseSymMatProd<double>;
    T_A_op A_op(A, B);
    T_B_op B_op(B);
    size_t dof = A.rows();
    Spectra::SymGEigsShiftSolver<T_A_op, T_B_op, Spectra::GEigsMode::ShiftInvert> eigs(A_op, B_op, num_eig, std::min(num_eig * 2, dof - 1), shift);
    // Spectra::SymGEigsShiftSolver<T_A_op, T_B_op, Spectra::GEigsMode::Buckling> eigs(A_op, B_op, num_eig, std::min(num_eig * 2, dof - 1), shift);

    eigs.init();
    int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 100000, tol, Spectra::SortRule::SmallestMagn);
    vals = eigs.eigenvalues();
    vecs = eigs.eigenvectors();
}
void solve_minimal_eig(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& diag_B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, double shift)
{
    if (A.rows() != A.cols() || A.rows() != diag_B.size()) {
        throw std::invalid_argument("SparseSymShiftSolve: matrix A must be square, and dof should be the same as diag_B");
    }
    const size_t dof = A.rows();

    Eigen::VectorXd B_squared_root = diag_B.array().sqrt().inverse();

    // compute B^{1/2} A B^{1/2}
    Eigen::SparseMatrix<double> A_galerkin = A;
#pragma omp parallel for
    for (int k = 0; k < A_galerkin.outerSize(); ++k) {
        auto B_squared_root_col = B_squared_root(k);
        for (Eigen::SparseMatrix<double>::InnerIterator it(A_galerkin, k); it; ++it) {
            it.valueRef() *= B_squared_root(it.row()) * B_squared_root_col;
            if (it.row() == k) {
                it.valueRef() += shift;
            }
        }
    }

    Spectra::SparseSymCholmodShiftSolve<double> A_op(A_galerkin);
    Spectra::SymEigsShiftSolver<Spectra::SparseSymCholmodShiftSolve<double>> eigs(A_op, num_eig, std::min(num_eig * 2, dof - 1), 0);
    eigs.init();

    int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 100000, 1e-6, Spectra::SortRule::SmallestMagn);

    if (eigs.info() == Spectra::CompInfo::Successful) {
        // std::cout << nconv << "eigen values and vectors have been found found by spectra." << "\n";
        vals = eigs.eigenvalues();
        vals.array() -= shift;
        vecs = eigs.eigenvectors();

// scale vecs back
#pragma omp parallel for
        for (size_t i = 0; i < vecs.rows(); ++i) {
            vecs.row(i) *= B_squared_root(i);
        }

    } else {
        std::cerr << "[EPS Error]: Eigen Value Problem Solving failed.\n";
        throw;
    }
}

void solve_maximal_eig(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& diag_B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, size_t maxit, double tol)
{
    if (A.rows() != A.cols() || A.rows() != diag_B.size()) {
        throw std::invalid_argument("SparseSymShiftSolve: matrix A must be square, and dof should be the same as diag_B");
    }
    const size_t dof = A.rows();

    Eigen::VectorXd B_squared_root = diag_B.array().sqrt().inverse();

    // compute B^{1/2} A B^{1/2}
    Eigen::SparseMatrix<double> A_galerkin = A;
#pragma omp parallel for
    for (int k = 0; k < A_galerkin.outerSize(); ++k) {
        auto B_squared_root_col = B_squared_root(k);
        for (Eigen::SparseMatrix<double>::InnerIterator it(A_galerkin, k); it; ++it) {
            it.valueRef() *= B_squared_root(it.row()) * B_squared_root_col;
        }
    }

    using T_op = Spectra::SparseSymMatProd<double>;
    T_op A_op(A_galerkin);
    Spectra::SymEigsSolver<T_op> eigs(A_op, num_eig, std::min(num_eig * 2, dof - 1));
    eigs.init();
    eigs.compute();

    int nconv = eigs.compute(Spectra::SortRule::LargestMagn, maxit, tol, Spectra::SortRule::SmallestMagn);

    if (eigs.info() == Spectra::CompInfo::Successful) {
        // std::cout << nconv << "eigen values and vectors have been found found by spectra." << std::endl;
        vals = eigs.eigenvalues();
        vecs = eigs.eigenvectors();

// scale vecs back
#pragma omp parallel for
        for (size_t i = 0; i < vecs.rows(); ++i) {
            vecs.row(i) *= B_squared_root(i);
        }

    } else {
        std::cerr << "[EPS Error]: Eigen Value Problem Solving failed.\n";
        throw;
    }
}

void solve_maximal_eig_sparseB(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, size_t maxit, double tol)
{
    using T_A_op = Spectra::SparseSymMatProd<double>;
    using T_B_op = Spectra::SparseCholesky<double>;
    T_A_op A_op(A);
    T_B_op B_op(B);

    const size_t dof = A.rows();
    Spectra::SymGEigsSolver<T_A_op, T_B_op, Spectra::GEigsMode::Cholesky> eigs(A_op, B_op, num_eig, std::min(num_eig * 2, dof - 1));
    eigs.init();
    int nconv = eigs.compute(Spectra::SortRule::LargestMagn, maxit, tol, Spectra::SortRule::SmallestMagn);
    vals = eigs.eigenvalues();
    vecs = eigs.eigenvectors();
}

void solve_minimal_eig_sparseB(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, double shift, int maxits, double tol, double* t_fac, size_t* nsolve)
{
    // using T_A_op = Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using T_A_op = Spectra::SymShiftInvertWithCholmod<double, Eigen::Sparse, Eigen::Sparse>;
    using T_B_op = Spectra::SparseSymMatProd<double>;
    T_A_op A_op(A, B);
    T_B_op B_op(B);
    size_t dof = A.rows();
    Spectra::SymGEigsShiftSolver<T_A_op, T_B_op, Spectra::GEigsMode::ShiftInvert> eigs(A_op, B_op, num_eig, std::min(num_eig * 2, dof - 1), shift);
    eigs.init();
    int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 100000, tol, Spectra::SortRule::SmallestMagn);
    vals = eigs.eigenvalues();
    vecs = eigs.eigenvectors();
    if (t_fac) {
        *t_fac = A_op.fac_time_cost;
    }
    if (nsolve) {
        *nsolve = A_op.cnt;
    }
}
void solve_maximal_eig(const Eigen::SparseMatrix<double>& A, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, size_t maxit, double tol)
{

    const size_t dof = A.rows();

    // compute B^{1/2} A B^{1/2}

    using T_op = Spectra::SparseSymMatProd<double>;
    T_op A_op(A);
    Spectra::SymEigsSolver<T_op> eigs(A_op, num_eig, std::min(num_eig * 2, dof - 1));
    eigs.init();
    int nconv = eigs.compute(Spectra::SortRule::LargestMagn, maxit, tol, Spectra::SortRule::SmallestMagn);

    if (eigs.info() == Spectra::CompInfo::Successful) {
        // std::cout << nconv << "eigen values and vectors have been found found by spectra." << "\n";
        vals = eigs.eigenvalues();
        vecs = eigs.eigenvectors();

    } else {
        std::cerr << "[EPS Error]: Eigen Value Problem Solving failed.\n";
        // std::cerr << A << "\n";
        throw;
    }
}
double solve_maximal_eig(const Eigen::SparseMatrix<double>& A, size_t maxit, double tol)
{

    Eigen::MatrixXd vec(A.rows(), 1);
    Eigen::VectorXd val(1);
    solve_maximal_eig(A, 1, val, vec, maxit, tol);
    return val(0);
}
void solve_minimal_eig(const Eigen::SparseMatrix<double>& A, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, double shift)
{
    const size_t dof = A.rows();

    Eigen::SparseMatrix<double> A_shift = A;
    A_shift.diagonal() += shift * Eigen::VectorXd::Ones(A.rows());

    Spectra::SparseSymCholmodShiftSolve<double> A_op(A_shift);
    Spectra::SymEigsShiftSolver<Spectra::SparseSymCholmodShiftSolve<double>> eigs(A_op, num_eig, std::min(num_eig * 2, dof - 1), 0);
    eigs.init();

    int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 100000, 1e-6, Spectra::SortRule::SmallestMagn);

    if (eigs.info() == Spectra::CompInfo::Successful) {
        // std::cout << nconv << "eigen values and vectors have been found found by spectra." << "\n";
        vals = eigs.eigenvalues();
        vals.array() -= shift;
        vecs = eigs.eigenvectors();

    } else {
        std::cerr << "[EPS Error]: Eigen Value Problem Solving failed.\n";
        throw;
    }
}
std::pair<double, double> get_specturm_span(const Eigen::SparseMatrix<double>& A, size_t offset)
{
    std::pair<double, double> low_up;
    Eigen::MatrixXd vec(A.rows(), 1);
    Eigen::VectorXd val(1);
    solve_maximal_eig(A, 1, val, vec);
    low_up.second = val(0);
    if (offset != 0) {
        vec.resize(A.rows(), offset + 1);
        val.resize(offset + 1);
    }
    solve_minimal_eig(A, offset + 1, val, vec);
    low_up.first = val(offset);
    return low_up;
}
