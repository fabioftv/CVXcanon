// Interface for Cone Solvers
//
// TODO(mwytock): This should likely be a pure C interface to ease
// implementation for existing cone solvers (ECOS, SCS)

// See: http://mathprogbasejl.readthedocs.io/en/latest/conic.html

// Specifical of cone problems and interface for cone solvers.

#ifndef CVXCANON_SOLVER_CONE_CONE_SOLVER_H
#define CVXCANON_SOLVER_CONE_CONE_SOLVER_H

#include <Eigen/Sparse>
#include <vector>
#include "cvxcanon/solver/SolverStatus.hpp"
#include "cvxcanon/util/Utils.hpp"

// TODO(fabioftv): Need another version of SymbolicConeSolver to account for other Cones
// A cone constraint
class ConeConstraint {
public:
	enum Cone {
		FREE,						// No Restrictions
		ZERO,						// All Components = 0
		NON_NEGATIVE,				// Nonnegative Orthant
		NON_POSITIVE,				// Nonpositive Orthant
		SECOND_ORDER,				// Second Order Cone
		ROTATED_SECOND_ORDER,		// Rotated Second Order Cone
		SEMIDEFINITE,				// Symmetric Positive Semidefinite Matrices
		PRIMAL_EXPO,				// Primal Exponential Cone
		DUAL_EXPO					// Dual Exponential Cone
	};

	Cone cone;
	int offset, size;
};

class ConeProblem {

// SCS
// Minimize		c'x
// Subject to:	Ax + s = b
//				s in K


public:
	SparseMatrix A;
	DenseVector b, h, c;
	std::vector<ConeConstraint> constraints;
};

// The solution to a cone problem
class ConeSolution {
 public:
  SolverStatus status;

  // Primal and dual variables
  DenseVector x, y;

  // Primal objective value
  double objective_value;
};

// The cone solver interface.
class ConeSolver {
public:
	virtual ConeSolution solve(const ConeProblem& problem) = 0;
	virtual ~ConeSolver() {}
};

#endif  // CVXCANON_SOLVER_CONE_CONE_SOLVER_H
