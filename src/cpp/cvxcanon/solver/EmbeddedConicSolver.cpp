
#include "EmbeddedConicSolver.hpp"

#include <unordered_map>

#include "cvxcanon/util/MatrixUtil.hpp"

// ECOS Environment
namespace ecos {
#include <ecos/src/ecos.c>
}
// TODO (fabioftv)
//		Can I include <ecos/src/ecos.c>? Reason: checkExitConditions(pwork_, mode)
//		Include Mixed-Integer SOCP Module in the Future

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

void EmbeddedConicSolver::build_ecos_leq_constraint(
	const Eigen::SparseMatrix<double, Eigen::RowMajor>&G,
	const DenseVector& h,
	const std::vector<ConeConstraint>& constraints_leq,
	int* total_size_leq,
	int* sizes_leq){
		int j = 0;
		for (const ConeConstraint& constr : constraints_leq){
			append_block_triplets(G.middleRows(constr.offset_leq, constr.size_leq), num_leq_constrs_, 0, &G_coeffs_);
			h_.segment(num_leq_constrs_, constr.size_leq) = h.segment(constr.offset_leq, constr.size_leq);
			if (total_size_leq != nullptr) *total_size_leq += constr.size_leq;
			if (sizes_leq != nullptr) sizes_leq[j++] = constr.size_leq;
			num_leq_constrs_ += constr.size_leq;
		}
}

void EmbeddedConicSolver::build_ecos_problem(const ConeProblem& problem, ConeSolution* solution){
	const int m_eq = problem.A.rows();
	const int n_eq = problem.A.cols();
	const int m_leq = problem.G.rows();
	const int n_leq = problem.G.rows();
	{
		Eigen::SparseMatrix<double, Eigen::RowMajor> A = problem.A;
		Eigen::SparseMatrix<double, Eigen::RowMajor> G = problem.G;
		const DenseVector& b = problem.b;
		const DenseVector& h = problem.h;

		A_coeffs_.clear();
		G_coeffs_.clear();
		b_ = DenseVector(m_eq);
		h_ = DenseVector(m_leq);
		num_eq_constrs = num_leq_constrs = 0;.

// TODO (fabioftv)
// Equality Constraint Assignment
//		std::unordered_map<int, std::vector<ConeConstraint>> constr_eq_map;
//		build_ecos_eq_constraint(A, b, constraints_eq, &pwork.p, nullptr);		

			std::unordered_map<int, std::vector<ConeConstraint>> constr_leq_map;
		for (const ConeConstraint& constr : problem.constraints_leq){
			constr_leq_map[constr.cone].push_back(constr);
		}

		build_ecos_leq_constraint(G, h, constr_leq_map[ConeConstraint::NON_NEGATIVE], &pwork.m, nullptr);

		A_ = sparse_matrix(m_eq, n_eq, A_coeffs_);
		G_ = sparse_matrix(m_leq, n_leq, G_coeffs_);
	}

// Dimension
	pwork_.n = std::max(n_eq, n_leq); // Can n_eq and n_leq be different?
	pwork_.m = m_leq;
	pwork_.p = m_eq;

// Variables
	solution->x = DenseVector(A.cols()); // Or DenseVector(G.cols()) => If n_eq and n_leq cannot differe, then I think it doesn't really matter
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
// 1) Set Parameters from <ecos/include/ecos.h>
//		pwork_.D => degree of the cone
//		pwork_.y => multipliers for equality constaints
//		pwork_.z => multipliers for conic inequalities
//		pwork_.lambda => scaled variable
//		pwork_.kap => kappa (homogeneous embedding)
//		pwork_.tau => tau (homogeneous embedding)
// 2) Check Parameters from "best iterate seen so far" and "temporary stuff holding search direction"
// 3) Check Othee Parameters: "indices that map entries of A and G to the KKT matrix", "equilibration vector", "scalings of problem data", "residuals", "norm iterates", and "KKT System"

}

// TODO(fabioftv): Consider Infeasible, Unbounded, and User Limit
SolverStatus EmbeddedConicSolver::get_ecos_status() {
	if (ecos::idxint checkExitConditions(pwork_, mode) == 0) {
		return OPTIMAL;
	} else {
		return ERROR;
	}
}

ConeSolution EmbeddedConicSolver::solve(const ConeProblem& problem) {
	ConeSolution solution;
	build_ecos_problem(problem, &solution);
// TODO(fabioftv): Check if ecos::ECOS_solve(&pwork_) is correct
	ecos::ECOS_solve(&pwork_);
	solution.objective_value = pwork_.cx;
	solution.x = pwork_.x;
	solution.status = get_ecos_status();
	return solution;
}