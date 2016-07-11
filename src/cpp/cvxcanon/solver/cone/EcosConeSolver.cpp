
#include "cvxcanon/solver/cone/EcosConeSolver.hpp"
#include "cvxcanon/util/MatrixUtil.hpp"
#include <unordered_map>
#include <vector>
#include <assert.h>

#define LEN_EXP 3

// ECOS Environment 
namespace ecos {
#include <ecos/external/SuiteSparse_config/SuiteSparse_config.h>
#include <ecos/include/ecos.h>
#include <ecos/include/spla.h>
#include <ecos/include/cone.h>
typedef double pfloat;
}  // namespace ecos

struct EcosConeSolver::EcosData {
   // ECOS Data Structures
   ecos::settings settings_;
   ecos::stats stats_;
   ecos::pwork* pwork_;
};

EcosConeSolver::EcosConeSolver() : ecos_data_(new EcosData()) {}
EcosConeSolver::~EcosConeSolver() {}


//(fabioftv): Determine the size of constraints
int EcosConeSolver::get_constr_size(
   const std::vector<ConeConstraint>& constraints) {
   int size_constraint = 0;
   for (const ConeConstraint& constr : constraints) {
      size_constraint += constr.size;
   }
   return size_constraint;
}

//(fabioftv): Build the constraints. If equality constraints, then the coefficients are assigned to A_coeffs_. If not, the coefficients are assigned to G_coeffs_.
void EcosConeSolver::build_ecos_constraint(
   const Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
   const DenseVector& b,
   const std::vector<ConeConstraint>& constraints,
   int offset) {
   for (const ConeConstraint& constr : constraints) {
      if (constr.cone != ConeConstraint::ZERO){
         append_block_triplets(
            A.middleRows(constr.offset, constr.size), offset, 0, 
               &G_coeffs_);
         h_.segment(offset, constr.size) = 
            b.segment(constr.offset, constr.size);
         offset += constr.size;
      }
      else {
         append_block_triplets(
               A.middleRows(constr.offset, constr.size), offset, 0, 
               &A_coeffs_);
         b_.segment(offset, constr.size) = 
            b.segment(constr.offset, constr.size);
         offset += constr.size;
      }
   }
}

//(fabioftv): Build the problem
void EcosConeSolver::build_ecos_problem(
   const ConeProblem& problem,
   ConeSolution* solution) {
   Eigen::SparseMatrix<double, Eigen::RowMajor> A = problem.A;
   const DenseVector& b = problem.b;

   std::unordered_map<int, std::vector<ConeConstraint>> constr_map;
   for (const ConeConstraint& constr : problem.constraints) {
      constr_map[constr.cone].push_back(constr);
   }
   //(fabioftv): Define the size of each set of constraints
   num_eq_constrs_ = get_constr_size(constr_map[ConeConstraint::ZERO]);
   num_leq_constrs_ = get_constr_size(constr_map[ConeConstraint::NON_NEGATIVE]);
   std::vector<ConeConstraint> seco_constrs = constr_map[ConeConstraint::SECOND_ORDER];
   int len_seco_constrs_ = get_constr_size(seco_constrs);
   num_seco_constrs_ = seco_constrs.size();
   num_exp_constrs_ = get_constr_size(constr_map[ConeConstraint::PRIMAL_EXPO])/LEN_EXP;

   const int n = problem.A.cols();
   const int m = num_leq_constrs_ + len_seco_constrs_ + num_exp_constrs_*LEN_EXP;
   const int p = num_eq_constrs_;
   // Count the number of second order cones.
   q_.clear();
   for (const ConeConstraint& constr : seco_constrs) {
     q_.push_back(constr.size);
   }

   //(fabioftv): Initialize variables and determine size of b_, h_, and s_
   A_coeffs_.clear();
   G_coeffs_.clear();
   b_ = DenseVector(p);
   h_ = DenseVector(m);
   //(fabioftv): Build the constraints
   build_ecos_constraint(A, b, constr_map[ConeConstraint::ZERO], 0);
   build_ecos_constraint(A, b, constr_map[ConeConstraint::NON_NEGATIVE], 0);
   build_ecos_constraint(A, b, constr_map[ConeConstraint::SECOND_ORDER],
                         num_leq_constrs_);
   build_ecos_constraint(A, b, constr_map[ConeConstraint::PRIMAL_EXPO],
                         num_leq_constrs_ + num_seco_constrs_);

   //(fabioftv): Assign coefficients to A_ and G_
   A_ = sparse_matrix(p, n, A_coeffs_);
   G_ = sparse_matrix(m, n, G_coeffs_);

   VLOG(1) << "ECOS Constraints:\n"
           << "A:\n" << matrix_debug_string(A_)
           << "b:\n" << vector_debug_string(b_)
           << "G:\n" << matrix_debug_string(G_)
           << "h:\n" << vector_debug_string(h_);

//(fabioftv): Fill out the information needed by ECOS.
// See more at: https://github.com/embotech/ecos/blob/develop/include/ecos.h
   // TODO fill out these arguments.
   ecos_data_->pwork_ = ecos::ECOS_setup(n, m, p,
                                  num_leq_constrs_,
                                  num_seco_constrs_, q_.data(), num_exp_constrs_,
                                  G_.valuePtr(), G_.outerIndexPtr(),
                                  G_.innerIndexPtr(),
                                  A_.valuePtr(), A_.outerIndexPtr(),
                                  A_.innerIndexPtr(),
                                        const_cast<double *>(problem.c.data()),
                                        const_cast<double *>(h_.data()),
                                        const_cast<double *>(b_.data()));
   //(fabioftv): Dimensions
   ecos_data_->settings_.verbose = 0;
}

//(fabioftv): Determine the status of the problem
SolverStatus EcosConeSolver::get_ecos_status(int exitflag) {
   if (exitflag == 0) {
       return OPTIMAL;
   } else if (exitflag == 1) {
      return INFEASIBLE;
   } else if (exitflag == 2) {
      return UNBOUNDED;
   } else if (exitflag == -1) {
      return USER_LIMIT;
   } else {
      return ERROR;
   }
}

//(fabioftv): Determine the solution of the problem.
ConeSolution EcosConeSolver::solve(const ConeProblem& problem) {
   ConeSolution solution;
   build_ecos_problem(problem, &solution);
   int exitflag = ECOS_solve(ecos_data_->pwork_);
   solution.p_objective_value = ecos_data_->pwork_->info->pcost;
   solution.d_objective_value = ecos_data_->pwork_->info->dcost;
   solution.status = get_ecos_status(exitflag);
   return solution;
}
