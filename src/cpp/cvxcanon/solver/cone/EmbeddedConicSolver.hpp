// Interface to the Embedded Conic Solver (ECOS)
//
// See: https://github.com/embotech/ecos

//(fabioftv) Copy from TODO(mwytock) in SCS: 
// It is expected that this interface will be implemented
// directly in the SCS code base at some point and called via a plugin
// architecture, see comment about pure C interface in README.md.

#ifndef CVXCANON_SOLVER_CONE_EMBEDDED_CONIC_SOLVER_H
#define CVXCANON_SOLVER_CONE_EMBEDDED_CONIC_SOLVER_H

#include <memory>
#include <vector>
#include "cvxcanon/solver/cone/ConeSolver.hpp"

// ECOS Environment 
namespace ecos {
	typedef double pfloat;
	typedef int idxint;
#include <ecos/include/ecos.h>
}

class EmbeddedConicSolver : public ConeSolver {
public:
	ConeSolution solve(const ConeProblem& problem) override;
private:
	void build_ecos_problem(const ConeProblem& problem, ConeSolution* solution);
	void build_ecos_eq_constraint(
		const Eigen::SparseMatrix<double, Eigen::RowMajor>&A,
		const DenseVector& b,
		const std::vector<ConeConstraint>& constraints_eq,
		int* total_size_eq,
		int* sizes_eq);
	void build_ecos_leq_constraint(
		const Eigen::SparseMatrix<double, Eigen::RowMajor>&G,
		const DenseVector& h,
		const std::vector<ConeConstraint>& constraints_leq,
		int* total_size_leq,
		int* sizes_leq);

	SolverStatus get_ecos_status();

// ECOS Data Structures
	ecos::settings settings_;
	ecos::stats stats_;
	ecos::pwork pwork_;

// ECOS Supporting Data Structures
	ecos::spmat A_matrix_;
	ecos::spmat G_matrix_;

// Constraints
	SparseMatrix A_;
	SparseMatrix G_;
	DenseVector b_;
	DenseVector h_;

// Extra Attributes for Constraints
	int num_eq_constrs_;
	int num_leq_constrs_;
	std::vector<Triplet> A_coeffs_;
	std::vector<Triplet> G_coeffs_;

// Cone
	ecos::cone cone_;

#endif  // CVXCANON_SOLVER_CONE_EMBEDDED_CONIC_SOLVER_H
