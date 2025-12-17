#include "eig.hpp"
#include <string>
#include "save_spm.hpp"     
#include "timer.hpp"
int main(int argc, char** argv){


    // #pragma omp parallel for
    // for(int i=0; i<100; i++){
    //     // 仅作并行测试
    //     std::cout << i << "\n";
    // }
    const std::string K_path = argv[1]; // e.g., "./data/K_block_0.csc"
    const std::string M_path = argv[2]; // e.g., "./data/M
    int nev = std::stoi(argv[3]); // e.g., 100  
    Eigen::SparseMatrix<double> K, M;
    zcy::io::read_spm(K_path.c_str(), K);
    zcy::io::read_spm(M_path.c_str(), M);

    Eigen::MatrixXd eig_vecs;
    Eigen::VectorXd eig_vals;
    zcy::TIMING::TIC("Solve Eigenvalues");
    solve_minimal_eig_sparseB(K, M, nev, eig_vals, eig_vecs, 0, 10000, 1e-8);
    zcy::TIMING::TOC("Solve Eigenvalues", true);
    return 0;
}