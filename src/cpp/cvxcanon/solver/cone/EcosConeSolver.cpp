
#include "cvxcanon/solver/cone/EcosConeSolver.hpp"
#include "cvxcanon/util/MatrixUtil.hpp"
#include <unordered_map>
#include <vector>

// ECOS Environment 
namespace ecos {
#include <ecos/external/SuiteSparse_config/SuiteSparse_config.h>
#include <ecos/include/ecos.h>
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
};

EcosConeSolver::EcosConeSolver() : ecos_data_(new EcosData()) {}
EcosConeSolver::~EcosConeSolver() {}
/*
void EcosConeSolver::build_ecos_constraint(
   const Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
   const DenseVector& b,
   const std::vector<ConeConstraint>& constraints,
   int* total_size,
   int* sizes) {
   for (const ConeConstraint& constr : constraints) {
      if (constr.cone != ZERO){
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
      if (total_size != nullptr) *total_size += constr.size;
      num_constrs_ += constr.size;        
  }
}

void EcosConeSolver::build_ecos_problem(
   const ConeProblem& problem,
   ConeSolution* solution) {
   
   const int m = problem.A.rows();
   const int n = problem.A.cols();

   {
      Eigen::SparseMatrix<double, Eigen::RowMajor> A = problem.A;
      const DenseVector& b = problem.b;

      A_coeffs_.clear();
      G_coeffs_.clear();
      // (fabioftv): Need to correctly define size of b_ and h_
      b_ = DenseVector(m);
      h_ = DenseVector(m);
      num_constrs_ = 0;

      std::unordered_map<int, std::vector<ConeConstraint>> constr_map;
      for (const ConeConstraint& constr : problem.constraints) {
         constr_map[constr.cone].push_back(constr);
      }

      // (fabioftv): Need to correctly define pwork_.m for each cone
      build_ecos_constraint(
         A, b, constr_map[ConeConstraint::ZERO], &ecos_data_->pwork_.p, nullptr);
      build_ecos_constraint(
         A, b, constr_map[ConeConstraint::NON_NEGATIVE], &ecos_data_->pwork_.m,
            nullptr);
      build_ecos_constraint(
         A, b, constr_map[ConeConstraint::SECOND_ORDER], &ecos_data_->pwork_.m,
            nullptr);
      build_ecos_constraint(
         A, b, constr_map[ConeConstraint::PRIMAL_EXPO], &ecos_data_->pwork_.m,
            nullptr);

      // (fabioftv): Need to correctly define size of A_ and G_
      A_ = sparse_matrix(m, n, A_coeffs_);
      G_ = sparse_matrix(m, n, G_coeffs_);

      VLOG(1) << "ECOS Constraints:\n"
              << "A:\n" << matrix_debug_string(A_)
              << "b:\n" << vector_debug_string(b_)
              << "G:\n" << matrix_debug_string(G_)
              << "h:\n" << vector_debug_string(h_);
   }
}*/
/*
scs_data_->A_matrix_.m = m;
  scs_data_->A_matrix_.n = n;
  scs_data_->A_matrix_.p = A_.outerIndexPtr();
  scs_data_->A_matrix_.i = A_.innerIndexPtr();
  scs_data_->A_matrix_.x = A_.valuePtr();

  scs_data_->data_.m = m;
  scs_data_->data_.n = n;
  scs_data_->data_.A = &scs_data_->A_matrix_;
  scs_data_->data_.b = const_cast<double*>(b_.data());
  scs_data_->data_.c = const_cast<double*>(problem.c.data());
  scs_data_->data_.stgs = &scs_data_->settings_;
  scs::setDefaultSettings(&scs_data_->data_);
  scs_data_->settings_.verbose = 0;

  s_ = DenseVector(A_.rows());
  solution->x = DenseVector(A_.cols());
  solution->y = DenseVector(A_.rows());
  scs_data_->sol_.x = const_cast<double*>(solution->x.data());
  scs_data_->sol_.y = const_cast<double*>(solution->y.data());
  scs_data_->sol_.s = const_cast<double*>(s_.data());



   // Dimension
   pwork_.n = std::max(n_eq, n_leq); // Can n_eq and n_leq be different?
   pwork_.m = m_leq;
   pwork_.p = m_eq;

// Variables
	solution->x = DenseVector(A.cols()); // Or DenseVector(G.cols()) => If n_eq and n_leq cannot differ, then I think it doesn't really matter
	pwork_.x = const_cast<double*>(solution->x.data());
	s_ = DenseVector(G.rows());

// Cone
	pwork_.cone const_cast<cone_*>(constr.cone());

// Problem Data
	pwork_.A = &A_matrix_;
	pwork_.G = &G_matrix_;
	pwork_.c = const_cast<double*>(problem.c.data());
	pwork_.b = const_cast<double*>(b_.data());
	pwork_.h = const_cast<double*>(h_.data());

// Info Struct
	pwork_.info = &stats_;

// Settings Struct
	pwork_.stgs = &settings_;

// TODO (fabioftv)
// 1) Check:
//		pwork_.D => degree of the cone
//		pwork_.y => multipliers for equality constaints
//		pwork_.z => multipliers for conic inequalities
}

// TODO(fabioftv): Check Description of the Function
SolverStatus EmbeddedConicSolver::get_ecos_status() {
	if (strcmp(stats_.info, "Solved") == 0) {
		return OPTIMAL;
	} else if (strcmp(stats_.info, "Infeasible") == 0) {
		return INFEASIBLE;
	} else if (strcmp(stats_.info, "Unbounded") == 0) {	
		return UNBOUNDED;
	} else if (strcmp(stats_.info, "User Limit") == 0) {
		return USER_LIMIT;
	else {
		return ERROR;
	}
}


*/







ConeSolution EcosConeSolver::solve(const ConeProblem& problem) {
  ConeSolution solution;


  // TODO(mwytock): solve using ecos, e.g. ECOS_setup()...

  ECOS_solve(&ecos_data_->pwork_);

  return solution;
}
