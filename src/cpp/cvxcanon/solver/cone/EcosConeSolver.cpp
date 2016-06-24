
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

void EcosConeSolver::build_ecos_constraint(
   const Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
   const Eigen::SparseMatrix<double, Eigen::RowMajor>& G,
   const DenseVector& b,
   const DenseVector& h,
   const std::vector<ConeConstraint>& constraints,
   int* total_size,
   int* sizes) {
   int j = 0;
   for (const ConeConstraint& constr : constraints) {
      append_block_triplets(
         A.middleRows(constr.offset, constr.size), num_constrs_, 0, &A_coeffs_);
      b_.segment(num_constrs_, constr.size) = 
         b.segment(constr.offset, constr.size);
    if (total_size != nullptr) *total_size += constr.size;
    if (sizes != nullptr) sizes[j++] = constr.size;
    num_constrs_ += constr.size;
  }
}

/*
void EmbeddedConicSolver::build_ecos_eq_constraint(
	const Eigen::SparseMatrix<double, Eigen::RowMajor>&A,
	const DenseVector& b,
	const std::vector<ConeConstraint>& constraints_eq,
	int* total_size_eq,
	int* sizes_eq){
		int j = 0;
		for (const ConeConstraint& constr : constraints_eq){
			append_block_triplets(A.middleRows(constr.offset_eq, constr.size_eq), num_eq_constrs_, 0, &A_coeffs_);
			b_.segment(num_eq_constrs_, constr.size_eq) = b.segment(constr.offset_eq, constr.size_eq);
			if (total_size_eq != nullptr) *total_size_eq += constr.size_eq;
			if (sizes_eq != nullptr) sizes_eq[j++] = constr.size_eq;
			num_eq_constrs_ += constr.size_eq;
		}
}
*/









ConeSolution EcosConeSolver::solve(const ConeProblem& problem) {
  ConeSolution solution;

  ecos::settings settings;
  ecos::pwork pwork;
  ecos::stats stats;

  // TODO(mwytock): solve using ecos, e.g. ECOS_setup()...

  ECOS_solve(&pwork);

  return solution;
}


