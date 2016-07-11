// Interface to the Embedded Conic Solver (ECOS)
//
// See: https://github.com/embotech/ecos
//
// TODO(mwytock): It is expected that this interface will be implemented
// directly in the SCS code base at some point and called via a plugin
// architecture, see comment about pure C interface in README.md.

#ifndef CVXCANON_SOLVER_CONE_ECOS_CONE_SOLVER_H
#define CVXCANON_SOLVER_CONE_ECOS_CONE_SOLVER_H

#include <memory>
#include <vector>
#include "cvxcanon/solver/cone/ConeSolver.hpp"

class EcosConeSolver : public ConeSolver {
public:
   EcosConeSolver();
   ~EcosConeSolver();

   ConeSolution solve(const ConeProblem& problem) override;

//(fabioftv): The following should be private:

   SolverStatus get_ecos_status(int exitflag);

   int get_constr_size(const std::vector<ConeConstraint>& constraints);

   void build_ecos_problem(const ConeProblem& problem, ConeSolution* solution);

   void build_ecos_constraint(
      const Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
      const DenseVector& b,
      const std::vector<ConeConstraint>& constraints,
      int offset);

   std::vector<Triplet> A_coeffs_;
   std::vector<Triplet> G_coeffs_;

private:
   struct EcosData;

   // ECOS Data Structures
   std::unique_ptr<EcosData> ecos_data_;

   // ECOS Supporting Data Structures
   std::vector<long> q_;

   // Constraints Ordered
   Eigen::SparseMatrix<double, Eigen::RowMajor, long> A_;
   Eigen::SparseMatrix<double, Eigen::RowMajor, long> G_;
   DenseVector b_;
   DenseVector h_;

   // Extra Attributes for Constraints
   int num_eq_constrs_;
   int num_leq_constrs_;
   int num_seco_constrs_;
   int num_exp_constrs_;
};

#endif  // CVXCANON_SOLVER_CONE_ECOS_CONE_SOLVER_H
