// Interface for Cone Solvers
//
// TODO(mwytock): This should likely be a pure C interface to ease
// implementation for existing cone solvers (ECOS, SCS)

// See: http://mathprogbasejl.readthedocs.io/en/latest/conic.html

#ifndef CONE_SOLVER_H
#define CONE_SOLVER_H

#include <Eigen/Sparse>

#include "cvxcanon/solver/SolverStatus.hpp"
#include "cvxcanon/util/Utils.hpp"

// TODO(fabioftv): Need another version of SymbolicConeSolver to account for other Cones
class ConeConstraint {
public:
	enum Cone {
		FREE,						// No Restrictions
		ZERO,						// All Components = 0
		NON_NEGATIVE,				// Nonnegative Orthant
		NON_POSITIVE,				// Nonpositive Orthant
		SECOND_ORDER,				// Second Order Cone
		ROTATED_SECOND_ORDER,		// Rotated Second Order Cone
		SYM_POS_SEMI,				// Symmetric Positive Semidefinite Matrices
		PRIMAL_EXPO,				// Primal Exponential Cone
		DUAL_EXPO					// Dual Exponential Cone
	};

	Cone cone;
	int offset_eq, offset_leq, size_eq, size_leq;
};

class ConeProblem {

// SCS
// Minimize		c'x
// Subject to:	Ax + s = b
//				s in K

// ECOS
// Minimize		c'x
// Subject to:	Ax = b
//				h - Gx in K

public:
	SparseMatrix A;
	SparseMatrix G;
	DenseVector b, h, c;
	std::vector<ConeConstraint> constraints_eq;
	std::vector<ConeConstraint> constraints_leq;
};

class ConeSolution {
public:
	SolverStatus status;
	DenseVector x, y;
	double objective_value;
};

class ConeSolver {
public:
	virtual ConeSolution solve(const ConeProblem& problem) = 0;
	virtual ~ConeSolver() {}
};

#endif  // CONE_SOLVER_H
