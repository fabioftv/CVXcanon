// Interface for All Solvers

#ifndef SOLVER_H
#define SOLVER_H

#include <map>

#include <Eigen/Sparse>

#include "cvxcanon/expression/Expression.hpp"
#include "cvxcanon/solver/SolverStatus.hpp"

//TODO(fabioftv): List All Solvers Here in the Future

class SolverOptions {
public:
	enum SolverOptions {
	ECOS,
	SCS,
	};
};

class Solution {
 public:
  SolverStatus status;
  std::map<int, DenseVector> variable_values;
  double objective_value;
};

class Solver {
 public:
  virtual ~Solver() {}
  virtual Solution solve(const Problem& problem) = 0;
};

Solution solve(const Problem& problem, const SolverOptions& solver_options);

#endif  // SOLVER_H
