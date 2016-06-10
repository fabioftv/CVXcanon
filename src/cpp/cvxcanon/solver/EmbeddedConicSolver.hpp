// Interface to the Embedded Conic Solver (ECOS)
//
// See: https://github.com/embotech/ecos

#ifndef EMBEDDED_CONIC_SOLVER_H
#define EMBEDDED_CONIC_SOLVER_H

#include <memory>

#include "cvxcanon/solver/ConeSolver.hpp"

namespace ecos {
typedef double pfloat;
typedef int idxint;
#include <ecos/include/ecos.h>
}

class EmbeddedConicSolver : public ConeSolver {
public:
	ConeSolution solve(const ConeProblem& problem) override;
private:
	void build_ecos_problem(
		const ConeProblem& problem,
		ConeSolution* solution);
	void build_ecos_equa_constraint(
		const Eigen::SparseMatrix<double, Eigen::RowMajor>&A,
		const DenseVector& b,
		int* total_size_equa,
		int* sizes_equa);
	void build_ecos_leq_constraint(
		const Eigen::SparseMatrix<double, Eigen::RowMajor>&G,
		const DenseVector& h,
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
	int num_equa_constrs_;
	std::vector<Triplet> A_coeffs_;
	std::vector<Triplet> G_coeffs_;

#endif  // EMBEDDED_CONIC_SOLVER_H
