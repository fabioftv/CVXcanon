
#include "Solver.hpp"

#include <glog/logging.h>

#include "cvxcanon/expression/TextFormat.hpp"
#include "cvxcanon/solver/SymbolicConeSolver.hpp"
#include "cvxcanon/solver/cone/EmbeddedConicSolver.hpp"
#include "cvxcanon/solver/cone/SplittingConeSolver.hpp"
#include "cvxcanon/transform/LinearConeTransform.hpp"

bool validate(const Problem& problem, const SolverOptions& solver_options) {
  LinearConeTransform transform;
  return transform.accepts(problem);
}

Solution solve(const Problem& problem, const SolverOptions& solver_options) {

// TODO(mwytock): Allow for different transforms/solvers as per SolveOptions

	VLOG(1) << "input problem:\n" << format_problem(problem);
	LinearConeTransform transform;
	Problem cone_problem = transform.apply(problem);

	VLOG(1) << "cone problem:\n" << format_problem(cone_problem);

  std::unique_ptr<ConeSolver> cone_solver;
// (fabioftv): Added Different Solvers
  if (solver_options.name == SCS) {
    SymbolicConeSolver solver(std::make_unique<SplittingConeSolver>());
    return solver.solve(cone_problem);
  } else if (solver_options.name == ECOS) {
    SymbolicConeSolver solver(std::make_unique<EmbeddedConicSolver>());
    return solver.solve(cone_problem);
  }
}
