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

// Cone Constraint
class ConeConstraint {
public:
	enum Cone {
		FREE,			// No Restrictions
		ZERO,			// All Components = 0
		NON_NEGATIVE,		// Nonnegative Orthant
		NON_POSITIVE,		// Nonpositive Orthant
		SECOND_ORDER,		// Second Order Cone
		ROTATED_SECOND_ORDER,	// Rotated Second Order Cone
		SEMIDEFINITE,		// Semidefinite Matrices
		PRIMAL_EXPO,		// Primal Exponential Cone
		DUAL_EXPO,		// Dual Exponential Cone
	};

	Cone cone;
	int offset, size;
};

// Standard Form Problem (Primal)
// Minimize	c'x
// Subject to:	b - Ax in K1
//		x in K2

class ConeProblem {
public:
	SparseMatrix A;
	DenseVector b, c;
	std::vector<ConeConstraint> constraints;
};

// Solution to Cone Problem
class ConeSolution {
 public:
  SolverStatus status;

  // Primal and Dual Variables
  DenseVector x, y;

  // Primal/Dual Objective Value
  double p_objective_value;
  double d_objective_value;
};

// Cone Solver Interface
class ConeSolver {
public:
	virtual ConeSolution solve(const ConeProblem& problem) = 0;
	virtual ~ConeSolver() {}
};

#endif  // CVXCANON_SOLVER_CONE_CONE_SOLVER_H
