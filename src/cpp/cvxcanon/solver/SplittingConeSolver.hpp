// Interface to the Splitting Conic Solver (SCS)
//
// See: https://github.com/cvxgrp/scs

#ifndef SPLITTING_CONIC_SOLVER_H
#define SPLITTING_CONIC_SOLVER_H

#include <memory>

#include "cvxcanon/solver/ConeSolver.hpp"

// SCS Environment
namespace scs {
	typedef double scs_float;
	typedef int scs_int;
#include <scs/linsys/amatrix.h>
#include <scs/include/scs.h>
}

class SplittingConeSolver : public ConeSolver {
public:
	ConeSolution solve(const ConeProblem& problem) override;
private:
	void build_scs_problem(const ConeProblem& problem, ConeSolution* solution);
	void build_scs_constraint(
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
		const DenseVector& b,
		const std::vector<ConeConstraint>& constraints,
		int* total_size,
		int* sizes);

	SolverStatus get_scs_status();

// SCS Data Structures
	scs::Data data_;
	scs::Settings settings_;
	scs::Sol sol_;
	scs::Info info_;
	scs::Scaling scaling_;
	scs::Work work;
	scs::Cone cone_;

// SCS Supporting Data Structures
	scs::AMatrix A_matrix_;
	std::unique_ptr<int[]> q_;
	std::unique_ptr<int[]> sd_;
	DenseVector s_;

// Constraints Ordered the way SCS needs them
	SparseMatrix A_;
	DenseVector b_;

// Used for Building Constraints
	int num_constrs_;
	std::vector<Triplet> A_coeffs_;

};

#endif  // SPLITTING_CONIC_SOLVER_H
