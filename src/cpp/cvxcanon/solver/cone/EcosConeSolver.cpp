
#include "cvxcanon/solver/cone/EcosConeSolver.hpp"
#include "cvxcanon/util/MatrixUtil.hpp"
#include <unordered_map>
#include <vector>

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
   ecos::pwork pwork_;

   // ECOS Supporting Data Structures
   ecos::spmat A_matrix_;
   ecos::spmat G_matrix_;
   ecos::cone cone_;
};

EcosConeSolver::EcosConeSolver() : ecos_data_(new EcosData()) {}
EcosConeSolver::~EcosConeSolver() {}


//(fabioftv): Determine the size of constraints
void EcosConeSolver::define_size_ecos_constraint(
   const std::vector<ConeConstraint>& constraints,
   int size_constraint) {
   size_constraint = 0;
   for (const ConeConstraint& constr : constraints) {
      size_constraint += constr.size;
   }
}

//(fabioftv): Build the constraints. If equality constraints, then the coefficients are assigned to A_coeffs_. If not, the coefficients are assigned to G_coeffs_.
void EcosConeSolver::build_ecos_constraint(
   const Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
   const DenseVector& b,
   const std::vector<ConeConstraint>& constraints) {
   for (const ConeConstraint& constr : constraints) {
      if (constr.cone != ConeConstraint::ZERO){
         append_block_triplets(
            A.middleRows(constr.offset, constr.size), num_constrs_, 0, 
               &G_coeffs_);
         h_.segment(num_constrs_, constr.size) = 
            b.segment(constr.offset, constr.size);
      }
      else {
         append_block_triplets(
            A.middleRows(constr.offset, constr.size), num_constrs_, 0, 
               &A_coeffs_);
         b_.segment(num_constrs_, constr.size) = 
            b.segment(constr.offset, constr.size);
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
   define_size_ecos_constraint(constr_map[ConeConstraint::ZERO],  
      num_eq_constrs_);
   define_size_ecos_constraint(constr_map[ConeConstraint::NON_NEGATIVE],
      num_leq_constrs_);
   define_size_ecos_constraint(constr_map[ConeConstraint::SECOND_ORDER],
      num_seco_constrs_);
   define_size_ecos_constraint(constr_map[ConeConstraint::PRIMAL_EXPO],
      num_exp_constrs_);

   const int n = problem.A.cols();
   const int m = num_leq_constrs_ + num_seco_constrs_ + num_exp_constrs_;
   const int p = num_eq_constrs_;

   //(fabioftv): Initialize variables and determine size of b_, h_, and s_
   A_coeffs_.clear();
   G_coeffs_.clear();
   b_ = DenseVector(p);
   h_ = DenseVector(m);
   s_ = DenseVector(m);

   //(fabioftv): Build the constraints
   build_ecos_constraint(A, b, constr_map[ConeConstraint::ZERO]);
   build_ecos_constraint(A, b, constr_map[ConeConstraint::NON_NEGATIVE]);
   build_ecos_constraint(A, b, constr_map[ConeConstraint::SECOND_ORDER]);
   build_ecos_constraint(A, b, constr_map[ConeConstraint::PRIMAL_EXPO]);

   //(fabioftv): Assign coefficients to A_ and G_
   A_ = sparse_matrix(p, n, A_coeffs_);
   G_ = sparse_matrix(m, n, G_coeffs_);

   VLOG(1) << "ECOS Constraints:\n"
           << "A:\n" << matrix_debug_string(A_)
           << "b:\n" << vector_debug_string(b_)
           << "G:\n" << matrix_debug_string(G_)
           << "h:\n" << vector_debug_string(h_);

//(fabioftv): Fill ou the information needed by ECOS.
// See more at: https://github.com/embotech/ecos/blob/develop/include/ecos.h
   
   //(fabioftv): Dimensions
   ecos_data_->pwork_.n = n;
   ecos_data_->pwork_.m = num_leq_constrs_;
   ecos_data_->pwork_.p = p;
   
   //(fabioftv): Variables
   ecos_data_->pwork_.x = const_cast<double*>(solution->x.data());
   ecos_data_->pwork_.z = const_cast<double*>(solution->y.data());
   ecos_data_->pwork_.y = const_cast<double*>(solution->y.data());
   ecos_data_->pwork_.s = const_cast<double*>(s_.data());

   //(fabioftv): Problem Data
   ecos_data_->A_matrix_.m = p;
   ecos_data_->A_matrix_.n = n;
   ecos_data_->A_matrix_.ir = (long int*) A_.innerIndexPtr();
   ecos_data_->A_matrix_.jc = (long int*) A_.outerIndexPtr();
   ecos_data_->A_matrix_.pr = A_.valuePtr();
   ecos_data_->A_matrix_.nnz = A_.nonZeros();

   ecos_data_->G_matrix_.m = m;
   ecos_data_->G_matrix_.n = n;
   ecos_data_->G_matrix_.ir = (long int*) G_.innerIndexPtr();
   ecos_data_->G_matrix_.jc = (long int*) G_.outerIndexPtr();
   ecos_data_->G_matrix_.pr = G_.valuePtr();
   ecos_data_->G_matrix_.nnz = G_.nonZeros();

   ecos_data_->pwork_.A = &ecos_data_->A_matrix_;
   ecos_data_->pwork_.G = &ecos_data_->G_matrix_;
   ecos_data_->pwork_.c = const_cast<double*>(problem.c.data());
   ecos_data_->pwork_.b = const_cast<double*>(b_.data());
   ecos_data_->pwork_.h = const_cast<double*>(h_.data());

   ecos_data_->cone_.nsoc = num_seco_constrs_;
   ecos_data_->cone_.nexc = num_exp_constrs_;

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

//(fabioftv): Determine the solution of the problem  
ConeSolution EcosConeSolver::solve(const ConeProblem& problem) {
   ConeSolution solution;
   build_ecos_problem(problem, &solution);
   int exitflag = ECOS_solve(&ecos_data_->pwork_);
   solution.p_objective_value = ecos_data_->stats_.pcost;
   solution.d_objective_value = ecos_data_->stats_.dcost;
   solution.status = get_ecos_status(exitflag);
   return solution;
}
