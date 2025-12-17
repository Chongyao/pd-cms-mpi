#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <optional>
/* Solve generalized eigen problem:
   Ax = wBx.

   Here we assume A is real,symmetric, sparse matrix, B is a positive diagonal matrix.
   The above eigen problem could be simplified as simplier eigen problem:
   B^{-1/2} A B^{-1/2} y = w y

   where y = B^{1/2}x,
      => x = B^{-1/2}y.

   To find the minial eigen values, we could use shift-invert method to solve
   it. The method needs a inverse solver, we use cholesky solver from cholmod to
   accelerate it.
 */
void solve_minimal_eig(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& diag_B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, double shift = 1e-4);
void solve_minimal_eig(const Eigen::SparseMatrix<double>& A, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, double shift = 1e-4);

void solve_maximal_eig(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& diag_B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, size_t maxit = 100000, double tol = 1e-10);
void solve_maximal_eig(const Eigen::SparseMatrix<double>& A, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, size_t maxit = 100000, double tol = 1e-10);
double solve_maximal_eig(const Eigen::SparseMatrix<double>& A, size_t maxit = 100000, double tol = 1e-10);

void solve_minimal_eig(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, double shift, int maxit, double tol);
void solve_minimal_eig_sparseB(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, double shift, int maxit, double tol, double* t_fac = nullptr, size_t* nsolve = nullptr);

void solve_maximal_eig_sparseB(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, size_t maxit = 1000, double tol = 1e-10);

std::pair<double, double> get_specturm_span(const Eigen::SparseMatrix<double>& A, size_t offset = 0);
// void solve_minimal_eig(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, size_t num_eig, Eigen::VectorXd& vals, Eigen::MatrixXd& vecs, double shift, int maxits, double tol);
