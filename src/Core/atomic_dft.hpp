//
// Created by F.Moitzi on 10.01.2023.
//

#ifndef LSMS_SRC_CORE_ATOMIC_DFT_HPP_
#define LSMS_SRC_CORE_ATOMIC_DFT_HPP_

#include <vector>
#include <memory>

#include "SingleSite/AtomData.hpp"
#include "XCLDA.hpp"

namespace lsms {

constexpr int BROYDEN_MAX_STEP = 6;
constexpr int BROYDEN_ITER_RESET = 25;
constexpr double BROYDEN_W0 = 0.01;

class AtomicDFT {

 private:

  int max_iter;
  int max_eig_iter;
  double e_eig_tol;
  double e_tol;
  int iprint;
  double alpha;

  std::unique_ptr<lsms::XCBase> xc;

 public:

  explicit AtomicDFT(
      std::vector<int> functional = {0, 0, 0},
      int max_iter = 100,
      int max_eig_iter = 100,
      double e_tol = 1.0e-10,
      double e_eig_tol = 1.0e-12, int iprint = 0, double alpha = 0.05
  );

  std::tuple<std::vector<double>, double> solve(int Z,
               const std::vector<double> &r_mesh,
               double h,
               int N,
               Matrix<double> &density,
               std::vector<double> &tot_potential
  );

};

void generate_starting_potential(std::vector<double> &potential,
                                 const std::vector<double> &r_mesh,
                                 int core_charge,
                                 std::size_t end
);

double generate_density(
    const std::vector<double> &r_mesh,
    double h,
    std::size_t end,
    const std::vector<double> &potential,
    int core_charge,
    const std::vector<int> &n,
    const std::vector<int> &l,
    const std::vector<int> &spin,
    const std::vector<int> &kappa,
    const std::vector<double> &occupation,
    std::vector<double> &e_eig,
    Matrix<Real> &density,
    int max_eig_iter,
    double eig_tol
);

double total_energy(const std::vector<double> &r_mesh,
                    Matrix<Real> &rho,
                    std::vector<double> &v_hartree,
                    Matrix<Real> &e_xc,
                    std::vector<double> &v_pot,
                    double E_band,
                    double Z,
                    int N
);

} // lsms

#endif //LSMS_SRC_CORE_ATOMIC_DFT_HPP_
